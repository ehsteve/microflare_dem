PRO make_teem_table,epstein=epstein,n=n,force_table=force_table

;Utility to remake the AIA lookup table (teem_table.sav or teem_table_epstein.sav) if needed.
;Andrew Inglis - 2012/10/23

wave_ =['131','171','193','211','335','94'] 
nwave =n_elements(wave_) 
nfiles = fltarr(nwave)

IF keyword_set(epstein) THEN BEGIN
	teem_table='teem_table_epstein.sav'
ENDIF ELSE BEGIN
	teem_table='teem_table.sav'
ENDELSE

default,n,10


t_min = 6.0
t_max = 7.5
t_d = 0.05
telog = t_d * findgen((t_max - t_min)/t_d) + t_min

tsig_min = 0.01
tsig_max = 0.8;0.40
tsig_d = 0.01 ;0.005
tsig = tsig_d * findgen((tsig_max - tsig_min)/tsig_d) + tsig_min

f = file_search(teem_table)

IF (f[0] NE '') AND NOT keyword_set(force_table) THEN BEGIN
	print,''
	print,'---------------------------'
	print,'Warning: ',teem_table,' already exists! Aborting.'
	print,'To OVERWRITE the lookup table, call make_teem_table with the /FORCE_TABLE keyword'
	print,'---------------------------'
	print,''
	return

ENDIF ELSE BEGIN
	aia_teem_table, wave_, tsig, telog = telog, q94 = q94, teem_table, save_dir = save_dir, area = area, n=n,epstein=epstein
	IF (f[0] NE '') THEN BEGIN
		print,''
		print,'---------------------'
		print,'Notice: ',teem_table,' created. Previous existing file was overwritten.'
		print,'---------------------'
		print,' '
	ENDIF ELSE BEGIN
		print,''
		print,'---------------------'
		print,'Notice: ',teem_table,' created.'
		print,'---------------------'
		print,' '
	ENDELSE
ENDELSE


END