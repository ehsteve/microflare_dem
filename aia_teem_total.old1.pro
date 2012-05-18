pro aia_teem_total,fileset = fileset,fov = fov,npix,wave_,q94,teem_table,teem_map,teem_tot, mask = mask, filelist = filelist, PLOT = plot, SAVE_DIR = save_dir
;+
; Project     : AIA/SDO
;
; Name        : AIA_TEMPMAP 
;
; Category    : Data analysis   
;
; Explanation : calculates temperature map
;		based on single-Gaussian fit
;
; Syntax      : IDL>aia_teem_total,fileset,fov,wave_,teem_map
;
; Inputs      : fileset  = strarr(6) with filenames of 6 wavelength FITS images
;		fov(4)   = [i1,j1,i2,j2] pixel ranges of field-of-view
;               npix     = macropixel size 
;		wave_	 = strarr(6) with wavelengths in Angstroem
;       q94      = empirical correction factor of 94 A response
;       teem_table = savefile containing DEM lookup table
;       teem_map = savefile containing EM and Te maps
;		mask 	= given a mask, return the dem only from within it.
;
; Outputs     : teem_tot = savefile containing total DEM and fluxes
;
; History     :  9-Mar-2011, Version 1 written by Markus J. Aschwanden
; 			  :  10-May-2011, Version 2 added mask keyword by S. Christe
;
; Contact     : aschwanden@lmsal.com
;-

default, save_dir, ''
;_____________________TOTAL FLUX_________________________________________
nwave = n_elements(wave_)
flux_ = fltarr(nwave)

IF keyword_set(fileset) THEN BEGIN
	searchstring=fileset+'*'+wave_(0)+'*'
	file_iw=file_search(searchstring,count=nfiles)
ENDIF ELSE BEGIN 
	;need to sort the files in the proper order
	;s = fltarr(nwave)
	;FOR i = 0, n_elements(wave_)-1 DO s[i] = where(strmatch(filelist,'*' + num2str(wave_[i]) + '_.fts') EQ 1)
	file_iw = filelist
ENDELSE

FOR iw = 0, nwave-1 DO BEGIN
	 
	fits2map, file_iw[iw], map
	
	IF keyword_set(FOV) THEN BEGIN
		sub_map, map, smap, xrange = fov[0:1], yrange = fov[2:3]
		map = smap
	ENDIF
	
	dim	= size(map.data)

	IF keyword_set(mask) THEN BEGIN
		print,total(map.data)
		FOR k=0, dim[0] do begin
			FOR l=0, dim[1] do begin
				data[k,l]=data[k,l]*mask[k,l]
			ENDFOR
		ENDFOR
		mask_frac=total(mask)/(dim[0]*dim[1])
	ENDIF ELSE mask_frac = 1
	
	texp = map.dur
	flux_[iw] = total(map.data)/texp
	;flux_(iw)=total(data)/texp
	print, wave_[iw], flux_[iw]
ENDFOR

;_____________________AIA RESPONSE FUNCTION________________________
restore,save_dir + teem_table      ;-->wave_,q94,area,resp_corr,telog,dte,tsig,flux
nte	=n_elements(telog)

;_____________________READ DEM PER PIXEL__________________________ 
restore,save_dir + teem_map	;-->te_map,em_map,sig_map,chi_map,dateobs 
dim	=size(em_map)
nx	=dim(1)
ny	=dim(2)
IF NOT keyword_set(mask) THEN tempmask = replicate(1,nx,ny) ELSE tempmask = congrid(mask,nx,ny)
em_tot	=fltarr(nte)
em	=fltarr(nte)
te	=fltarr(nte)
sig	=fltarr(nte)
em_errtot=fltarr(nte)
error_array=fltarr(nx,ny,nte)
nmacro = float(nx)*float(ny)
FOR j=0,ny-1 do begin
 FOR i=0,nx-1 do begin
  em[*]  = tempmask[i,j]*10.^em_map[i,j]	;log(EM)-->EM
  te[*]  =te_map[i,j]*tempmask[i,j]	        ;log(te)
  sig[*] =sig_map[i,j]*tempmask[i,j]
  em_tot =em_tot+em*exp(-(telog-te)^2/(2.*sig^2))/(nmacro * mask_frac)
  IF (sig[0] gt 0.) THEN begin
		error1=double(tsig_errmap_symmetric[i,j])
  		error2=double(telog_errmap_symmetric[i,j])
		error3=(2.*(error2/telog-te)*(-(telog-te)^2.)) ; error in -(telog-te)^2 - careful with minus sign.
		error4=2.*(error1/sig)*(sig^2.) ; error in sig^2
		error5=2.*error4 ; error in 2sig^2
		error6=sqrt( (error3/(telog-te)^2.) + (error5/(2.*sig^2.)) ) * ((telog-te)^2. / (2.*sig^2.)) ;error in exponent: {(telog-te)^2/(2*sig^2)} 
		error7=error6 * exp(-(telog-te)^2./(2.*sig^2.)) ;error once the exp is done.
		error8=error7*(em/nmacro) ; final error in em at i,j
		;check FOR -NaN values
		p=finite(error8)
		index = WHERE (p gt 0,count,complement=index_c)
		comp_count = n_elements(p) - count
		;replace NaN values with 0.
		IF comp_count NE 0 THEN BEGIN 
			error8[index_c] = 0.
			;store the error FOR this i,j in an array.
			error_array[i,j,*]=error8
		ENDIF
	ENDIF ELSE begin
		error_array[i,j,*]=0.
	;em_errtot = em_errtot + 0.
	ENDELSE
	
 ENDFOR
ENDFOR
emlog = alog10(em_tot)


;get sum of the squares FOR each element of em_errtot = fltarr(nte). This gives the final error in em_tot summed over all the pixels.
FOR j = 0,nte-1 DO em_errtot[j] = reform(sqrt(total(error_array[*,*,j]^2.)))

;get error in em_log rather than em_tot
emlog_err = (em_errtot/em_tot)

clearplot
IF keyword_set(PLOT) THEN BEGIN
	plot,telog,emlog,yrange=minmax(emlog),xtitle='Temperature  log(T)',$
   		ytitle='Emission measure  log(EM [cm!U-5!N K!U-1!N])'
   	oploterr, telog, emlog, emlog_err
ENDIF

;______________________TOTAL FLUX FROM DEM________________________
flux_dem=fltarr(nwave)
qflux	=fltarr(nwave)
FOR iw=0,nwave-1 do begin
 flux_dem[iw]=total(resp_corr[*,iw] * em_tot*dte) * nmacro*npix^2 * mask_frac
 qflux[iw] = flux_dem[iw]/flux_[iw]
 print,flux_(iw),flux_dem(iw),qflux(iw)
ENDFOR

;______________________SAVE RESULTS_______________________________
save,filename=save_dir + teem_tot,telog,emlog,flux_,flux_dem,qflux, emlog_err
;save,filename='teem_tot_test'+strtrim(string(u),2)+'.sav',telog,emlog,emlog_err,flux_,flux_dem,qflux

END
