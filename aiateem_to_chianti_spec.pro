PRO aiateem_to_chianti_spec, x = x, y = y, filename = filename

restore, '/Users/schriste/Dropbox/idl/aia_deem/teem_data_00020110215_000009.sav', /verbose
restore, '/Users/schriste/Dropbox/idl/aia_deem/teem_table.sav'

;IF (NOT keyword_set(i) OR NOT keyword_set(j)) THEN BEGIN
	s = size(telog_map)
	dimx = s[1]
	dimy = s[2]
;ENDIF ELSE BEGIN
;	dimx = 1
;	dimy = 1
;ENDELSE

FOR i = 0, dimx-1 DO BEGIN
	FOR j = 0, dimy-1 DO BEGIN
	
	em_best = em_map[i,j]
	telog_best = telog_map[i,j]
	sig_best = sig_map[i,j]

	emlog =em_best*exp(-(telog-telog_best)^2/(2.*sig_best^2))

	plot, telog, emlog, yrange=minmax(emlog),xtitle='Temperature  log(T)',$
   		ytitle='Emission measure  log(EM [cm!U-5!N K!U-1!N])'
   
	r = chianti_spec_from_dem(telog, emlog, /plot)

	save, emlog, telog, r, filename = 'xray_spectrum_i' + num2str(i) + '_j' + num2str(j) + '.sav'
ENDFOR
	ENDFOR
END