pro aia_teem_table, wave_, tsig, TE_RANGE = te_range, TELOG = telog, Q94 = q94, teem_table, SAVE_DIR = save_dir, area = area
;+
; Project     : AIA/SDO
;
; Name        : AIA_DEM_LOOKUP 
;
; Category    : Data analysis   
;
; Explanation : calculates AIA fluxes for 6 wavelengths (WAVE_)
;		for single-Gaussian DEM distributions with
;		TEMIN < TE < TEMAX, DTE1 < DTE < DTE2
;
; Syntax      : IDL>aia_dem_table,wave_,tsig,te_range,q94,tem_table
;
; Inputs      : wave_ = strarr(6) with wavelengths in Angstroem
;
; Outputs     : postscript file <plotname>_col.ps (if io=2)
;
; History     :  3-Mar-2011, Version 1 written by Markus J. Aschwanden
;
; Contact     : aschwanden@lmsal.com
;-

default, save_dir, ''
default, wave_, ['131','171','193','211','335','94']

;_____________________AIA RESPONSE FUNCTION________________________
; old response
;tresp = aia_get_response(/temp,/full, /dn, version = 1)
tresp = aia_get_response(/temp, /full, /dn, version = 2, /chiantifix, /evenorm)

telog_  = tresp.logte
nwave = n_elements(wave_)
ichan_ = fltarr(nwave)

; put the response in the right order based on wave_
FOR iw = 0, nwave-1 DO BEGIN
	IF iw EQ 0 THEN BEGIN
		nte_ = n_elements(tresp.tresp[*,0])
		nwave_ = n_elements(tresp.tresp[0,*])
		resp_ = fltarr(nte_, nwave)
	ENDIF
	filter ='A'+wave_(iw)
	IF (wave_[iw] eq '094') THEN filter='A94'
	ichan = where(tresp.channels EQ filter)
	resp_[*,iw] = tresp.tresp[*,ichan]
ENDFOR

; if te_range is set, limit the response to that range
IF keyword_set(te_range) THEN BEGIN
	telog1	= alog10(te_range[0])
	telog2	= alog10(te_range[1])
	ind_te	= where((telog_ GE telog1) AND (telog_ LE telog2), nte)
	resp    = fltarr(nte,nwave)
	telog	= telog_[ind_te]
	FOR iw = 0, nwave-1 DO resp[*,iw] = resp_[ind_te, iw]
ENDIF ELSE IF NOT keyword_set(telog) THEN BEGIN
	nte = nte_
	resp    = fltarr(nte,nwave)
	telog = telog_
ENDIF

; if telog is set than need to interpolate from the tabulated values
; to the temperatures defined in telog
IF keyword_set(telog) THEN BEGIN
	nte = n_elements(telog)
	resp = fltarr(nte,nwave)
	FOR iw = 0, nwave-1 DO resp[*,iw] = interpol(resp_[*,iw], telog_, telog, /quadratic)
ENDIF

;_____________________EMPIRICAL CORRECTION 94 A____________________
resp_corr = resp
IF keyword_set(q94) THEN BEGIN
	ind1 = where(telog le 6.3)
	resp_corr[ind1,5] = resp[ind1,5] * q94
ENDIF ELSE q94 = 1

;_____________________CALCULATES LOOPUP TABLES_____________________
dte1 = 10.^telog(1:nte-1)-10.^telog(0:nte-2)
dte	= [dte1[0], dte1]
em1	= 1.
nsig = n_elements(tsig)
flux = fltarr(nte, nsig, nwave)

FOR i = 0, nte-1 DO BEGIN
	FOR j = 0, nsig-1 DO BEGIN
		em_te = em1 * exp(-(telog - telog[i])^2 / (2. * tsig[j]^2))
  		FOR iw = 0, nwave-1 DO flux[i,j,iw] = total(resp_corr[*,iw] * em_te * dte)
	ENDFOR
ENDFOR

save, filename = save_dir + teem_table, wave_, q94, area, resp_corr, telog, dte, tsig, flux
print, 'Lookup table created : ', teem_table

END