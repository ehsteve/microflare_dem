pro aia_teem_map,fileset = fileset,FOV = fov,wave_,npix,teem_table,teem_map, FILElist = filelist, VERBOSE = verbose,u

;+
; Project     : AIA/SDO
;
; Name        : AIA_TEMPMAP 
;
; Category    : Data analysis   
;
; Explanation : calculates EM and Te temperature maps
;		based on single-Gaussian fits in each macropixel
;
; Syntax      : IDL>aia_teem_map,fileset,fov,wave_,npix,teem_table,teem_map
;
; Inputs      : fileset  = common filename part of 6 wavelength FITS images
;		wave_	 = strarr(6) with wavelengths in Angstroem
;		fov(4)	 = [i1,j1,i2,j2] pixel ranges of field-of-view
;		npix	 = size of macropixel (spatial resolution)
;		teem_table = filename of DEM lookup table
;                            (calculated previously with AIA_TEEM_TABLE.PRO)
;		teem_map = savefile containing EM and Te maps 
;
; Outputs     : postscript file <plotname>_col.ps (if io=2)
;
; History     :  3-Mar-2011, Version 1 written by Markus J. Aschwanden
;				 16-May-2011, made fov optional. If not given, uses the
;								whole images. 
;								made fileset optional. Added filelist keyword to just give a list of files.
;		: 5-Oct-2011 - Added search for 1 sigma contour in chi^2 space to calculate errors in telog and tsig - A. Inglis
;
; Contact     : aschwanden@lmsal.com
;-

;_________________________________________________________________________
t1 = systime(0,/seconds)
nwave = n_elements(wave_)

IF keyword_set(fileset) THEN BEGIN
	files = strarr(6)
	FOR iw = 0, nwave-1 DO BEGIN
		 ;searchstring=fileset+'*'+wave_(iw)+'.fits'
		 searchstring=fileset+'*'+wave_[iw]+'_.fts'
	 	file_iw = file_search(searchstring, count=nfiles)
	 	files[iw]=file_iw[u]
	ENDFOR
ENDIF ELSE BEGIN 
	;need to sort the files in the proper order
	files = filelist
	s = fltarr(nwave)
	FOR i = 0, n_elements(wave_)-1 DO s[i] = where(strmatch(files,'*' + num2str(wave_[i]) + '*') EQ 1)
	files = files[s]
ENDELSE
;_____________________REBINNING IMAGES___________________________________

texp_	=fltarr(nwave)
for iw=0,nwave-1 do begin
	read_sdo,files(iw),index,data
	
	; schriste (16-may-2011)
	; If Keyword FOV not given then determine the FOV
	; of the whole image automatically
	IF keyword_set(fov) THEN BEGIN
		i1=fov[0] & j1=fov[1] & i2=fov[2] & j2=fov[3]
	ENDIF ELSE BEGIN
		s = size(data)
		i1=0 & j1=0 & i2=s[1]-1 & j2=s[2]-1
	ENDELSE
	
	dim	=size(data)
	nx0	=dim[1]
	ny0	=dim[2]
	if (i1 ge nx0) or (i2 ge nx0) or (j1 ge ny0) or (j2 ge ny0) then begin
		print,'FOV=',i1,i2,j1,j2
		print,'image size=',nx0,ny0
		stop,'ERROR in subimage range i1,i2,j1,j2 for image with size nx0,ny0
	endif
	image	=data[i1:i2,j1:j2]
	texp   =index.exptime
	dateobs=index.date_obs
	if (iw eq 0) then begin
		dim	=size(image)
		nx	=dim[1]
		ny	=dim[2]
		nxx	=(nx/npix)
		nyy	=(ny/npix)
		i3	=nxx*npix-1
		j3	=nyy*npix-1
		x	=i1+(npix*findgen(nxx)+0.5)
		y	=j1+(npix*findgen(nyy)+0.5)
		images=fltarr(nxx,nyy,nwave)
	endif
	if (npix eq 1) then images(*,*,iw)=float(image)/texp
	if (npix gt 1) then images(*,*,iw)=rebin(float(image[0:i3,0:j3]),nxx,nyy)/texp
	texp_(iw)=texp

	;PLOT the image as a MAP
	id = strmid(index.TELESCOP,0,3) + '/' + strmid(index.instrume,0,3)
	aia_lct, rr, gg, bb, wavelnth=index.WAVELNTH, /load
	a = index.CROTA2
	xcen   = index.CRVAL1 + index.CDELT1*cos(a)*(0.5*(index.NAXIS1+1)- index.CRPIX1)-index.CDELT2*sin(a)*((index.NAXIS2+1)/2-index.CRPIX2)
	ycen = index.CRVAL2 + index.CDELT1*sin(a)*((index.NAXIS1+1)*0.5-index.CRPIX1) + index.CDELT2*cos(a)*((index.NAXIS2+1)*0.5-index.CRPIX2)
	map = make_map(data, time = index.t_obs, id = id, dur = index.exptime, xc = xcen, yc = ycen, dx = index.cdelt1, dy = index.cdelt2)
	plot_map, map, /limb
	
endfor

;________________________TEMPERATURE MAP_________________________________ 
restore,teem_table	;-->wave_,q94,area,resp_corr,telog,dte,tsig,flux
dim	=size(flux)
nte	=dim[1]
nsig	=dim[2]
nwave	=dim[3]
ntot	=nte*nsig
te_map	=fltarr(nxx,nyy)
em_map	=fltarr(nxx,nyy)
sig_map	=fltarr(nxx,nyy)
chi_map	=fltarr(nxx,nyy)
telog_errmap=fltarr(nxx,nyy,2)
tsig_errmap=fltarr(nxx,nyy,2)
telog_errmap_symmetric=fltarr(nxx,nyy)
tsig_errmap_symmetric=fltarr(nxx,nyy)
chi_map6=fltarr(nxx,nyy,6)
te_best	=0.
em_best	=0.
sig_best=0.
r0	=0.95*(4096/2)
x0	=4096/2
y0	=4096/2
nfree	=3

FOR j = 0, nyy-1 DO BEGIN
 if (nx eq 4096) then ind=where(sqrt((x-x0)^2+(y(j)-y0)^2) le r0,nind)
 if (nx lt 4096) then ind=findgen(nxx)
 i1	= min(ind) > 0
 i2	= max(ind) < (nxx-1)
 for i = i1, i2 do begin
  flux_obs=reform(images[i,j,*])
  counts=flux_obs*texp_
  noise=sqrt(counts)/texp_
  chimin=9999.
  chi6min=9999.
  chi2d = fltarr(nte,nsig)
  for k=0,nte-1 do begin
   for l=0,nsig-1 do begin
    flux_dem1=reform(flux(k,l,*))
    em1	=total(flux_obs)/total(flux_dem1)
    flux_dem=flux_dem1*em1
    chi	=sqrt(total((flux_obs-flux_dem)^2/noise^2)/(nwave-nfree))
    chi6=abs(flux_obs-flux_dem)/noise
    chi2d[k,l] = chi
    if (chi le chimin) then begin
     chimin	=chi
     chi6min	=chi6
     em_best	=alog10(em1)
     te_best	=telog(k)
     sig_best	=tsig(l)
    endif
   endfor 
  endfor

; find errors in tsig and telog by finding the 1 sigma contour in chi^2 space for each pixel - A. Inglis 5-Oct-2011
IF (chimin ne 9999.) THEN BEGIN
	;chimin defaults to 9999. if bad pixel (i.e. no data), so skip those
	;find the 1 sigma contour in chi^2-space to get errors in tsig and telog
	contour, chi2d, telog, tsig, levels = [chimin +1.0],path_xy=path_xy,path_info=path_info,/path_data_coords,closed=0
	testxmax=MAX(path_xy(0,*))
	testxmin=MIN(path_xy(0,*))
	testymax=MAX(path_xy(1,*))
	testymin=MIN(path_xy(1,*))

	telog_err=[testxmin,testxmax] - te_best
	;error bars are asymmetric - symmetrise!
	telog_errsymmetric=sqrt(telog_err(0)^2 + telog_err(1)^2)

	tsig_err=[testymin,testymax] - sig_best
	;error bars are asymmetric - symmetrise!
	tsig_errsymmetric=sqrt(tsig_err(0)^2 + tsig_err(1)^2)
	;print,telog_errsymmetric
	;print,tsig_errsymmetric
	;print,telog_err
	;print,tsig_err
	ENDIF ELSE BEGIN
	;contingency statements in case bad pixel
	telog_err=[0,0]
	tsig_err=[0,0]
	tsig_errsymmetric=0.
	telog_errsymmetric=0.
ENDELSE 

	IF keyword_set(VERBOSE) THEN BEGIN
		nlevels = 20
		levels = chimin * (findgen(nlevels)*0.1 + 1.0)
		anot = strarr(20)
		FOR y = 0, nlevels-1 DO anot[y] = num2str(levels[y])
		contour, chi2d, telog, tsig, levels = levels, c_annotation = anot, xtitle = 'log[Temp]', ytitle = 't_sig'
		oplot, [te_best], [sig_best], psym = 4
		leg = num2str(chimin) + '[' + num2str(te_best) + ',' + num2str(sig_best) + ']'
		legend, leg, psym = 4
		contour,chi2d,telog,tsig,levels=[chimin + 1.],/over,thick=2
		stop
	
	ENDIF
  em_map(i,j)=em_best
  te_map(i,j)=te_best
  sig_map(i,j)=sig_best
  chi_map(i,j)=chimin
  chi_map6(i,j,*)=chi6min
  telog_errmap(i,j,0:1)=telog_err[0:1]
  tsig_errmap(i,j,0:1)=tsig_err[0:1]
  telog_errmap_symmetric(i,j)=telog_errsymmetric
  tsig_errmap_symmetric(i,j)=tsig_errsymmetric
  ;telog_errmap(i,j,*)=telog_err
  ;tsig_errmap(i,j,*)=tsig_err
 endfor
 if (j mod 10) eq 0 then print,j,nyy
endfor

help,telog_errmap_symmetric
help,tsig_errmap_symmetric

print,telog_errmap_symmetric(50:70,50:70)
print,tsig_errmap_symmetric(50:70,50:70)

aia_map_cube = make_map( images[*,*,0], xc = xcen, yc = ycen, dx = index.cdelt1*npix, dy = index.cdelt2*npix, id = '', time = dateobs)
aia_map_cube = replicate(aia_map_cube, 6)
	
FOR i = 0, n_elements(wave_)-1 DO BEGIN
	aia_map_cube[i].data = images[*,*,i]
	aia_map_cube[i].id = 'SDO/AIA ' + num2str(wave_[i]) + ' A'
ENDFOR

temperature_map = make_map( te_map, xc = xcen, yc = ycen, dx = index.cdelt1*npix, dy = index.cdelt2*npix, id = 'log(Temperature [MK])', time = dateobs)

emission_map = make_map( em_map, xc = xcen, yc = ycen, dx = index.cdelt1*npix, dy = index.cdelt2*npix, id = 'log(Emission Measure [cm!U-3!N]', time = dateobs)

sigma_map = make_map( sig_map, xc = xcen, yc = ycen, dx = index.cdelt1*npix, dy = index.cdelt2*npix, id = 'Sigma', time = dateobs)

chisq_map = make_map( chi_map, xc = xcen, yc = ycen, dx = index.cdelt1*npix, dy = index.cdelt2*npix, id = textoidl('\chi^2'), time = dateobs)

tsig_error_map_minus = make_map( tsig_errmap(*,*,0), xc = xcen, yc = ycen, dx = index.cdelt1*npix, dy = index.cdelt2*npix, id = textoidl('\chi^2'), time = dateobs)

telog_error_map_minus = make_map( telog_errmap(*,*,0), xc = xcen, yc = ycen, dx = index.cdelt1*npix, dy = index.cdelt2*npix, id = textoidl('\chi^2'), time = dateobs)

tsig_error_map_plus = make_map( tsig_errmap(*,*,1), xc = xcen, yc = ycen, dx = index.cdelt1*npix, dy = index.cdelt2*npix, id = textoidl('\chi^2'), time = dateobs)

telog_error_map_plus = make_map( telog_errmap(*,*,1), xc = xcen, yc = ycen, dx = index.cdelt1*npix, dy = index.cdelt2*npix, id = textoidl('\chi^2'), time = dateobs)

;________________________STATISTICS OF CHI-2 FITS_________________
print,'Statistics of chi2='
statistic,chi_map
print,'Statistics of chi2 (94 A)'
statistic,chi_map6(*,*,5)
print,'log(EM)-range = ',minmax(em_map)

;________________________SAVE MAPS________________________________
;save,filename=teem_map,te_map,em_map,sig_map,chi_map,dateobs 
save,filename=teem_map,te_map,em_map,sig_map,chi_map,tsig_errmap_symmetric,telog_errmap_symmetric,temperature_map,emission_map,sigma_map,chisq_map,aia_map_cube, dateobs 
save,temperature_map,emission_map,sigma_map,chisq_map,tsig_error_map_minus,tsig_error_map_plus,telog_error_map_minus,telog_error_map_plus,filename='andy_tempmaptesting' + strtrim(string(u),2) + '.sav',/verbose
print,'TE+EM maps saved in file : ',teem_map
t2	=systime(0,/seconds)
cpu	=t2-t1
print,'Computation time = ',cpu,' s'
end
