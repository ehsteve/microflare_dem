pro aia_teem_total,fileset,fov = fov,npix,wave_,q94,teem_table,teem_map,teem_tot, mask = mask, u
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
;			  :  5-Oct-2011 Added error calculations for em_tot based on errors in telog and sig - A. Inglis
;
; Contact     : aschwanden@lmsal.com
;-

;_____________________TOTAL FLUX_________________________________________
nwave	=n_elements(wave_)
flux_	=fltarr(nwave)

for iw=0,nwave-1 do begin
 ;searchstring=fileset+'*'+wave_(iw)+'.fits'
 ;file_iw=findfile(searchstring,count=nfiles)
 searchstring=fileset+'*'+wave_(iw)+'_.fts'
 file_iw=file_search(searchstring,count=nfiles)
 
 read_sdo,file_iw(u),index,data
 	IF keyword_set(fov) THEN BEGIN
		i1=fov[0] & j1=fov[1] & i2=fov[2] & j2=fov[3]
	ENDIF ELSE BEGIN
		s = size(data)
		i1=0 & j1=0 & i2=s[1]-1 & j2=s[2]-1
	ENDELSE
masktotal=0.
IF keyword_set(mask) then begin
tempmask=congrid(mask,i2+1,j2+1)

print,total(data)
for k=0,i2 do begin
for l=0,j2 do begin
data(k,l)=data(k,l)*tempmask(k,l)
IF (tempmask(k,l) eq 1) then masktotal=masktotal+1
endfor
endfor	
endif
 texp   =index.exptime
 flux_(iw)=total(data(i1:i2,j1:j2))/texp
 ;flux_(iw)=total(data)/texp
print,wave_(iw),flux_(iw)
endfor
mask_frac=masktotal/(i2*j2)


;_____________________AIA RESPONSE FUNCTION________________________
restore,teem_table      ;-->wave_,q94,area,resp_corr,telog,dte,tsig,flux
nte	=n_elements(telog)

;_____________________READ DEM PER PIXEL__________________________ 
restore,teem_map	;-->te_map,em_map,sig_map,chi_map,dateobs 
dim	=size(em_map)
nx	=dim(1)
ny	=dim(2)
print,telog_errmap_symmetric(50:70,50:70)
print,tsig_errmap_symmetric(50:70,50:70)
IF NOT keyword_set(mask) THEN mask = replicate(1,nx,ny)
em_tot	=fltarr(nte)
em_errtot=fltarr(nte)
error_array=fltarr(nx,ny,nte)
em	=fltarr(nte)
te	=fltarr(nte)
sig	=fltarr(nte)
nmacro	=float(nx)*float(ny)
for j=0,ny-1 do begin
 for i=0,nx-1 do begin
  em(*)  =mask[i,j]*10.^em_map[i,j]	;log(EM)-->EM
  te(*)  =te_map[i,j]*mask[i,j]	        ;log(te)
  sig(*) =sig_map[i,j]*mask[i,j]
  em_tot =em_tot+em*exp(-(telog-te)^2/(2.*sig^2))/nmacro

  	;estimate the error in em_tot based on the error in telog and tsig - ARI 2011/09/30
	;IF statement makes sure we don't count pixels where mask = 0.
	IF (sig(0) gt 0.) THEN begin
		error1=double(tsig_errmap_symmetric(i,j))
  		error2=double(telog_errmap_symmetric(i,j))
		error3=(2.*(error2/telog-te)*(-(telog-te)^2.)) ; error in -(telog-te)^2 - careful with minus sign.
		error4=2.*(error1/sig)*(sig^2.) ; error in sig^2
		error5=2.*error4 ; error in 2sig^2
		error6=sqrt( (error3/(telog-te)^2.) + (error5/(2.*sig^2.)) ) * ((telog-te)^2. / (2.*sig^2.)) ;error in exponent: {(telog-te)^2/(2*sig^2)} 
		error7=error6 * exp(-(telog-te)^2./(2.*sig^2.)) ;error once the exp is done.
		error8=error7*(em/nmacro) ; final error in em at i,j
		;check for -NaN values
		p=finite(error8)
		index = WHERE (p gt 0,count,complement=index_c)
		;replace NaN values with 0.
		error8(index_c)=0.
		;store the error for this i,j in an array.
		error_array(i,j,*)=error8
	
	ENDIF ELSE begin
		error_array(i,j,*)=0.
	;em_errtot = em_errtot + 0.
	ENDELSE
        ;DEBUGGING STATEMENTS
	;help,error1,error2,error3,error4,error5,error6,error7,error8
	;print,error1,error2,error3,error4,error5,error6,error7,error8
	;print,error8
	;print,em_errtot
	;IF (total(sig) gt 0.) THEN print,sig
	;IF (total(sig) gt 0.) THEN print,error3(0),error4(0),error5(0),error6(0),error7(0),error8(0),em_errtot(0)
 endfor
endfor
emlog	=alog10(em_tot)

for v=0,32 do begin 
	;get sum of the squares for each element of em_errtot = fltarr(nte). This gives the final error in em_tot summed over all the pixels.
	em_errtot(v)=reform(sqrt(total(error_array(*,*,v)^2.)))
endfor
;get error in em_log rather than em_tot
emlog_err = (em_errtot/em_tot)
;MORE DEBUGGING STATEMENTS
;help,emlog
;help,emlog_err
;print,emlog
;print,emlog_err
;stop
clearplot
window,0,xsize=756,ysize=756
plot,telog,emlog,yrange=[min(emlog)-max(emlog_err),max(emlog)+max(emlog_err)],xtitle='Temperature  log(T)',$
   ytitle='Emission measure  log(EM [cm!U-5!N K!U-1!N])'
;overplot the error bars
oploterr, telog, emlog, emlog_err
;______________________TOTAL FLUX FROM DEM________________________
flux_dem=fltarr(nwave)
qflux	=fltarr(nwave)
for iw=0,nwave-1 do begin
 flux_dem(iw)=total(resp_corr(*,iw)*em_tot*dte)*nmacro*npix^2; * mask_frac
 qflux(iw)   =flux_dem(iw)/flux_(iw)
 print,flux_(iw),flux_dem(iw),qflux(iw)
endfor

;______________________SAVE RESULTS_______________________________
save,filename='teem_tot_test'+strtrim(string(u),2)+'.sav',telog,emlog,emlog_err,flux_,flux_dem,qflux
end
