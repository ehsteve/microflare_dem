pro aia_teem_map,fileset = fileset,wave_,npix,teem_table,teem_map, FILElist = filelist, VERBOSE = verbose, filename_extra = FILENAME_extra, SAVE_DIR = save_dir, DEBUG = debug, XRANGE = xrange, YRANGE = yrange

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
; Syntax      : IDL> aia_teem_map,fileset,fov,wave_,npix,teem_table,teem_map
;
; Inputs      : fileset  = common filename part of 6 wavelength FITS images
;		wave_	 = strarr(6) with wavelengths in Angstroem
;		/// fov(4)	 = [i1,j1,i2,j2] pixel ranges of field-of-view - deprecated
;       xrange = fltarr(2) the xrange in arcsec for the subimage
;       yrange = fltarr(2) the yrange in arcsec for the subimage
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
;		: 1-Feb-2011 - Changed definition of files (see lines 51-55) so that files = filelist directly. Solves issues with AIA2011xxx file formats.
;
; Contact     : aschwanden@lmsal.com
;-

default, save_dir, ''

;_________________________________________________________________________
t1 = systime(0,/seconds)
nwave = n_elements(wave_)
default, filename_extra, ''

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
	;s = fltarr(nwave)
	;FOR i = 0, n_elements(wave_)-1 DO s[i] = where(strmatch(files,'*' + num2str(wave_[i]) + '_.fts') EQ 1)
	;files = files[s]
ENDELSE
;_____________________REBINNING IMAGES___________________________________

texp_ = fltarr(nwave)

FOR iw = 0, nwave-1 DO BEGIN
	read_sdo,files[iw],index,data
    index2map,index,data,map
		
	IF keyword_set(xrange) AND keyword_set(yrange) THEN BEGIN
	    sub_map, map, smap, xrange = xrange, yrange = yrange
	    map = smap
	    data = smap.data
    ENDIF
    
	s = size(data)
	i1=0 & j1=0 & i2=s[1]-1 & j2=s[2]-1

	dim	=size(data)
	nx0	=dim[1]
	ny0	=dim[2]

	image = map.data
	texp = map.dur
	dateobs = map.time
	
	IF (iw EQ 0) then begin
		dim	= size(image)
		nx	= dim[1]
		ny	= dim[2]
		nxx	= (nx/npix)
		nyy	= (ny/npix)
		i3	= nxx*npix-1
		j3	= nyy*npix-1
		;x	= i1+(npix*findgen(nxx)+0.5)
		;y	= j1+(npix*findgen(nyy)+0.5)
		images = fltarr(nxx,nyy,nwave)
		x = get_map_xp(map, /oned)
		y = get_map_yp(map, /oned)
	ENDIF
	
	IF (npix EQ 1) THEN images[*,*,iw] = float(image)/texp
	IF (npix GT 1) THEN BEGIN
		images[*,*,iw] = rebin(float(image[0:i3,0:j3]),nxx,nyy)/texp
		sub_map, map, smap, xrange = [0,i3], yrange = [0,j3], /pixel
		rmap = rebin_map(smap, nxx, nyy)
		x = get_map_xp(rmap, /oned)
		y = get_map_yp(rmap, /oned)
		map = rmap
	ENDIF
	
	texp_[iw]=texp 
    
    IF iw EQ 0 THEN aia_map_cube = replicate(map, nwave) ELSE aia_map_cube[iw] = map
        
	IF keyword_set(DEBUG) THEN BEGIN
		aia_lct, rr, gg, bb, wavelnth=wave_[iw], /load
		plot_map, map, /limb, /log
	ENDIF
ENDFOR

;________________________TEMPERATURE MAP_________________________________ 
restore, save_dir + teem_table	;-->wave_,q94,area,resp_corr,telog,dte,tsig,flux
dim	= size(flux)
nte	= dim[1]
nsig = dim[2]
nwave = dim[3]
ntot = nte*nsig
te_map = fltarr(nxx,nyy)
em_map = fltarr(nxx,nyy)
sig_map	= fltarr(nxx,nyy)
flux_dem_map = fltarr(nxx, nyy, nwave)
chi_map = fltarr(nxx,nyy)
em_errmap = fltarr(nxx,nyy)
telog_errmap = fltarr(nxx,nyy,2)
tsig_errmap = fltarr(nxx,nyy,2)
telog_errmap_symmetric = fltarr(nxx,nyy)
tsig_errmap_symmetric = fltarr(nxx,nyy)
chi_map6 = fltarr(nxx,nyy,6)
te_best	= 0.
em_best	= 0.
sig_best = 0.
;r0 = 0.95*(4096/2)
;x0 = 4096/2
;y0 = 4096/2
nfree = 3
;stop
temp = get_sun(dateobs)
solar_radius = temp[1]

FOR j = 0, nyy-1 DO BEGIN
	; the following code makes sure that regions with radial distance beyond r0
	; are not considered.
	ry = sqrt(x^2 + y[j]^2)
	ind = where(ry LE solar_radius, count)
	;if (nx eq 4096) then ind=where(sqrt((x-x0)^2+(y(j)-y0)^2) le r0,nind)
	;if (nx lt 4096) then ind=findgen(nxx)
	;i1	= min(ind) > 0
	;i2	= max(ind) < (nxx-1)
	if count NE 0 THEN BEGIN
		i1 = min(ind)
		i2 = max(ind) 
	ENDIF ELSE continue
	
	FOR i = i1, i2-1 DO BEGIN
		flux_obs=reform(images[i,j,*])
		counts=flux_obs*texp_
		noise=sqrt(counts)/texp_
		chimin = 9999.
		chi6min = 9999.
		flux_dem_best = 9999.
		chi2d = fltarr(nte,nsig)
		FOR k=0, nte-1 DO BEGIN
   			FOR l = 0, nsig-1 DO BEGIN
				flux_dem1=reform(flux(k,l,*))
				em1	=total(flux_obs)/total(flux_dem1)
				em1_err = sqrt(total(noise^2))/total(flux_dem1)
				flux_dem=flux_dem1*em1
				chi	=sqrt(total((flux_obs-flux_dem)^2/noise^2)/(nwave-nfree))
				chi6=abs(flux_obs-flux_dem)/noise
				chi2d[k,l] = chi				
				IF (chi le chimin) THEN BEGIN
					chimin	=chi
					chi6min	=chi6
					em_best	=alog10(em1)
					em_best_err = em1_err/em_best
					te_best	=telog(k)
					sig_best	=tsig(l)
					flux_dem_best = flux_dem
				ENDIF
			ENDFOR 
		ENDFOR

	; find errors in tsig and telog by finding the 1 sigma contour in chi^2 space for each pixel - A. Inglis 5-Oct-2011
	IF (chimin NE 9999.) THEN BEGIN
		;chimin defaults to 9999. if bad pixel (i.e. no data), so skip those
		;find the 1 sigma contour in chi^2-space to get errors in tsig and telog
		contour, chi2d, telog, tsig, levels = [chimin +2.3],path_xy=path_xy,path_info=path_info,/path_data_coords,closed=0
		IF exist(path_xy) THEN BEGIN
			testxmax=MAX(path_xy(0,*))
			testxmin=MIN(path_xy(0,*))
			testymax=MAX(path_xy(1,*))
			testymin=MIN(path_xy(1,*))
		
			telog_err=[testxmin,testxmax] - te_best
			;error bars are asymmetric - symmetrise!
			telog_errsymmetric=sqrt(telog_err(0)^2 + telog_err(1)^2)
		
			tsig_err=[sig_best - testymin,testymax - sig_best]
			;error bars are asymmetric - symmetrise!
			;tsig_errsymmetric=sqrt(tsig_err(0)^2 + tsig_err(1)^2)
			tsig_errsymmetric = max(tsig_err)
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
	ENDIF ELSE BEGIN
		;contingency statements in case bad pixel
		telog_err=[0,0]
		tsig_err=[0,0]
		tsig_errsymmetric=0.
		telog_errsymmetric=0.
	ENDELSE 

	IF keyword_set(VERBOSE) THEN BEGIN
		print, 'j = ', j
		print, 'chimin = ', chimin
		print, 'Wavelength:', wave_
		print, 'Observed flux:', flux_obs
		print, 'Fitted flux:', flux_dem_best
		print, 'Fitted/Obs flux:', flux_dem_best/flux_obs
		
		IF chimin LT 9999. THEN BEGIN
			loadct, 0
			hsi_linecolors
			nlevels = 20
			levels = chimin * (findgen(nlevels)*0.1 + 1.0)
			anot = strarr(20)
			FOR k = 0, nlevels-1 DO anot[k] = num2str(levels[k])
			contour, chi2d, telog, tsig, levels = levels, c_annotation = anot, xtitle = 'log[Temp]', ytitle = 't_sig'
			oplot, [te_best], [sig_best], psym = symcat(16), color = 6
			leg = num2str(chimin) + '[' + num2str(te_best) + ',' + num2str(sig_best) + ']'
			legend, leg, psym = 4
			contour,chi2d,telog,tsig,levels=[chimin + 2.3],/over,thick=2, color = 6
		ENDIF
		stop
	ENDIF
	;stop
	flux_dem_map[i,j,*] = flux_dem_best
	em_map[i,j]=em_best
	te_map[i,j]=te_best
	sig_map[i,j]=sig_best
	chi_map[i,j]=chimin
	chi_map6[i,j,*]=chi6min
	telog_errmap[i,j,0:1]=telog_err[0:1]
	tsig_errmap[i,j,0:1]=tsig_err[0:1]
	telog_errmap_symmetric[i,j]=telog_errsymmetric
	tsig_errmap_symmetric[i,j]=tsig_errsymmetric
	;telog_errmap(i,j,*)=telog_err
	;tsig_errmap(i,j,*)=tsig_err
	
	ENDFOR
 	IF (j mod 10) eq 0 THEN print,j,nyy
ENDFOR

help,telog_errmap_symmetric
help,tsig_errmap_symmetric

;print,telog_errmap_symmetric(50:70,50:70)
;print,tsig_errmap_symmetric(50:70,50:70)

aia_simul_map_cube = aia_map_cube

FOR i = 0, nwave-1 DO BEGIN
	aia_simul_map_cube[i].data = flux_dem_map[*,*,i]
	aia_simul_map_cube[i].id = 'Simulated SDO/AIA ' + num2str(wave_[i]) + ' A'
ENDFOR

temperature_map = map
temperature_map.data = te_map
temperature_map.id = 'log(Temperature [MK])'
temperature_map.time = dateobs

emission_map = map
emission_map.data = em_map
emission_map.id = 'log(Emission Measure [cm!U-3!N]'
emission_map.time = dateobs

sigma_map = map
sigma_map.data = sig_map
sigma_map.id = 'Sigma'
sigma_map.time = dateobs

chisq_map = map
chisq_map.data = chi_map
chisq_map.id = textoidl('\chi^2')
chisq_map.time = dateobs

tsig_error_map_minus = map
tsig_error_map_minus.data = tsig_errmap[*,*,0]
tsig_error_map_minus.id = 'Temp sigma err -'
tsig_error_map_minus.time = dateobs

telog_error_map_minus = map
telog_error_map_minus.data = telog_errmap[*,*,0]
telog_error_map_minus.id = 'log(Temp)  err -'
telog_error_map_minus.time = dateobs

tsig_error_map_plus = map
tsig_error_map_plus.data = tsig_errmap[*,*,1]
tsig_error_map_plus.id = 'Temp sigma err +'
tsig_error_map_plus.time = dateobs

tsig_error_map_plus = map
tsig_error_map_plus.data = telog_errmap[*,*,1]
tsig_error_map_plus.id = 'log(Temp)  err +'
tsig_error_map_plus.time = dateobs

;________________________STATISTICS OF CHI-2 FITS_________________
print,'Statistics of chi2='
statistic,chi_map
print,'Statistics of chi2 (94 A)'
statistic,chi_map6(*,*,5)
print,'log(EM)-range = ',minmax(em_map)

;________________________SAVE MAPS________________________________
;save,filename=teem_map,te_map,em_map,sig_map,chi_map,dateobs 
save,filename=save_dir + teem_map,te_map,em_map,sig_map,chi_map,tsig_errmap_symmetric,telog_errmap_symmetric,temperature_map,emission_map,sigma_map,chisq_map,aia_map_cube, dateobs , aia_simul_map_cube
save,temperature_map,emission_map,sigma_map,chisq_map,tsig_error_map_minus,tsig_error_map_plus,telog_error_map_minus,telog_error_map_plus,aia_simul_map_cube, filename=save_dir + 'tempmap' + filename_extra + '.sav'

print,'TE+EM maps saved in file : ',teem_map
t2	=systime(0,/seconds)

print,'Computation time = ',(t2-t1)/60.0,' min'

IF keyword_set(DEBUG) THEN stop

END
