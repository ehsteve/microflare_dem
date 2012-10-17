PRO make_teem_table,epstein=epstein,n=n,force_table=force_table


wave_ =['131','171','193','211','335','94'] 
nwave =n_elements(wave_) 
nfiles = fltarr(nwave)

teem_table='teem_table.sav'

default,n,1

;file_list = get_aia_file_list(dir, fileset = fileset)
;FOR i = 0, nwave-1 DO nfiles[i] = n_elements(file_list[*,i])

;hardcoded this for 21 June 2011 flare
marker=120;120

;area of 1 AIA pixel in cm^2
aia_pixel_area = 1.85589e+15


t_min = 5.5
t_max = 8.0

t_min = 5.5 ; 6.0
t_max = 7.5
t_d = 0.05
telog = t_d * findgen((t_max - t_min)/t_d) + t_min

tsig_min = 0.01
tsig_max = 1.0;0.40
tsig_d = 0.002 ;0.005
tsig = tsig_d * findgen((tsig_max - tsig_min)/tsig_d) + tsig_min

f = file_search(teem_table)

;area = aia_teem_pixel_area(file_list[0,0])

IF keyword_set(epstein) THEN BEGIN
	teem_table='teem_table_epstein.sav'
	IF f[0] EQ '' OR keyword_set(FORCE_TABLE) THEN aia_teem_table_epstein, wave_, tsig, telog = telog, q94 = q94, teem_table, save_dir = save_dir, n=n
ENDIF ELSE BEGIN
	IF f[0] EQ '' OR keyword_set(FORCE_TABLE) THEN aia_teem_table, wave_, tsig, telog = telog, q94 = q94, teem_table, save_dir = save_dir, area = area
ENDELSE


END