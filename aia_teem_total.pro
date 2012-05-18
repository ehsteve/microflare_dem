pro aia_teem_total,fileset = fileset,npix,wave_,q94,teem_table,teem_map,teem_tot, mask_map = mask_map, filelist = filelist, PLOT = plot, SAVE_DIR = save_dir, XRANGE = xrange, YRANGE = yrange, macro_dem = macro_dem
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
;		mask_map 	= given a mask as a map, return the dem only where it equals 1.
;
; Outputs     : teem_tot = savefile containing total DEM and fluxes
;
; History     :  9-Mar-2011, Version 1 written by Markus J. Aschwanden
; 			  :  10-May-2011, Version 2 added mask keyword by S. Christe
;				:1-Feb-2012 Version 2 - fixed outstanding bugs in code - A. Inglis
;				a) Lines 44-50 did not handle the AIA2011xxxx file format and have been replaced by setting 
;				file_iw=filelist directly (filelist is now in correct handling order already).
;				b)Error calculation was missing mask_frac on the denominator when error 8 is calculated (line 113)
;				c)mask_frac is now calculated correctly! See line 79.
;		:3-Feb-2012  fixed energy calculation, and now generates an energy_map with the energy in each pixel. Saved in teem_tot.....sav file
;
; Contact     : aschwanden@lmsal.com
;-
;
default, save_dir, ''
;_____________________TOTAL FLUX_________________________________________
nwave = n_elements(wave_)
flux_ = fltarr(nwave)

file_iw = filelist

FOR iw = 0, nwave-1 DO BEGIN
	 
	read_sdo,file_iw[iw],index,data
    index2map,index,float(data),map
			
	IF keyword_set(xrange) AND keyword_set(yrange) THEN BEGIN
	    sub_map, map, smap, xrange = xrange, yrange = yrange
	    map = smap
	    data = smap.data
    ENDIF
    
	s = size(data)
	i1=0 & j1=0 & i2=s[1]-1 & j2=s[2]-1

	image = data
	texp = map.dur
	dateobs = map.time	
	
	IF keyword_set(mask_map) then begin
		mask = mask_map.data
		; zero out everything that is not in the mask
		FOR k = 0, i2 DO BEGIN
			FOR l = 0, j2 DO BEGIN
				data[k,l] = data[k,l]*mask[k,l]
			ENDFOR
		ENDFOR	
	ENDIF
	
	flux_[iw] = total(data[i1:i2,j1:j2])/texp
	print,'Total flux in ', wave_(iw),flux_(iw)
ENDFOR

IF keyword_set(mask) THEN mask_frac=total(mask)/n_elements(data) ELSE mask_frac = 1.

;_____________________AIA RESPONSE FUNCTION________________________
restore,save_dir + teem_table      ;-->wave_,q94,area,resp_corr,telog,dte,tsig,flux
nte	=n_elements(telog)

;_____________________READ DEM PER PIXEL__________________________ 
restore,save_dir + teem_map	;-->te_map,em_map,sig_map,chi_map,dateobs 
dim	=size(em_map)
nx	=dim(1)
ny	=dim(2)
IF NOT keyword_set(mask) THEN tempmask = replicate(1,nx,ny) ELSE tempmask = congrid(mask,nx,ny)
em_tot	=dblarr(nte)
em	=dblarr(nte)
te	=dblarr(nte)
sig	=dblarr(nte)
em_errtot=dblarr(nte)
error_array=dblarr(nx,ny,nte)
energy=fltarr(nx,ny)
kb=1.3806503e-16
nmacro = float(nx)*float(ny)

FOR j=0,ny-1 do begin
 FOR i=0,nx-1 do begin
  em[*]  = tempmask[i,j]*10.^em_map[i,j]	;log(EM)-->EM
  te[*]  =telog_map[i,j]*tempmask[i,j]	        ;log(te)
  sig[*] =sig_map[i,j]*tempmask[i,j]
  em_tot =em_tot+em*exp(-(telog-te)^2/(2.*sig^2))/(nmacro * mask_frac)
  IF keyword_set(macro_dem) THEN BEGIN
		if ((i EQ 0) AND (j EQ 0)) THEN cnt = 0
		if cnt mod macro_dem THEN stop
	ENDIF
  	
  IF (sig[0] gt 0.) THEN begin
		error1=double(sig_errmap_symmetric[i,j])
  		error2=double(telog_errmap_symmetric[i,j])
		error3=(2.*(error2/telog-te)*(-(telog-te)^2.)) ; error in -(telog-te)^2 - careful with minus sign.
		error4=2.*(error1/sig)*(sig^2.) ; error in sig^2
		error5=2.*error4 ; error in 2sig^2
		error6=sqrt( (error3/(telog-te)^2.) + (error5/(2.*sig^2.)) ) * ((telog-te)^2. / (2.*sig^2.)) ;error in exponent: {(telog-te)^2/(2*sig^2)} 
		error7=error6 * exp(-(telog-te)^2./(2.*sig^2.)) ;error once the exp is done.
		error8=error7*(em/(nmacro*mask_frac)) ; final error in em at i,j 01/27/12 - added mask_frac - ARI
		;check FOR -NaN values
		p=finite(error8)
		index = WHERE (p gt 0,count,complement=index_c)
		;replace NaN values with 0.
		error8(index_c)=0.
		;store the error FOR this i,j in an array.
		error_array[i,j,*]=error8

		;calculate the energy in each pixel by integrating the emission measure. ARI 2011/11/18
		;1.8805014e11 converts from pixels to cm^2
		;Assume V=A^1.5
  	
		;Energy = 1.5k * V^1/2 * Integral(EM ^ 1/2 dT)

		;integrate the square root of emission measure. Need to convert to cm^-5 hence the /1.8805013e15
		area=1.8805013e15*npix^2
		emgaussian=em*exp(-(telog-te)^2/(2.*sig^2))/(nmacro * mask_frac)
		em_integral=int_tabulated((10^telog),sqrt(emgaussian*area),/double)	
		energy(i,j)=1.5 * em_integral * sqrt((area)^1.5) * kb 
		;print,'test numbers:  ',total(sqrt(em)),em_integral, sqrt((nmacro*npix^2*mask_frac*1.8805014e15)^1.5), kb, energy(i,j)
		
	ENDIF ELSE begin
		error_array[i,j,*]=0.
	;em_errtot = em_errtot + 0.
	ENDELSE
	;print,em,error1,error2,error3,error4,error5,error6,error7,error8
	;stop
 ENDFOR
ENDFOR
emlog = alog10(em_tot)


;get sum of the squares FOR each element of em_errtot = fltarr(nte). This gives the final error in em_tot summed over all the pixels.
FOR j = 0,nte-1 DO em_errtot[j] = reform(sqrt(total(error_array[*,*,j]^2.)))

;get error in em_log rather than em_tot
emlog_err = (em_errtot/em_tot)
energy_map=temperature_map
energy_map.data = energy
energy_map.id = 'Energy [erg]'
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
save,filename=save_dir + teem_tot,telog,emlog,flux_,flux_dem,qflux, emlog_err, error_array, energy_map
;save,filename='teem_tot_test'+strtrim(string(u),2)+'.sav',telog,emlog,emlog_err,flux_,flux_dem,qflux

END
