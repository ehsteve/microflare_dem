FUNCTION get_aia_file_list, dir, WAVE_LIST = wave_list ,fileset = fileset

; PURPOSE: Given a directory get all of the aia cutout files in the directory.
;
; KEYWORDS: 
;           wave_list - 
;
; WRITTEN: Steven Christe (8-Oct-2011)

default, wave_list, ['131','171','193','211','335','94'] 
default, fileset, 'ssw_cutout'

nwave =n_elements(wave_list) 
nfiles = fltarr(nwave)

; search for all of the files
FOR i = 0, nwave-1 DO BEGIN
	IF fileset EQ 'ssw_cutout' THEN files = file_search(dir + fileset + '*' + wave_list[i] + '_.fts')
	IF fileset EQ 'AIA' THEN files = file_search(dir + fileset + '*0' + wave_list[i] + '.fits')
	nfiles[i] = n_elements(files)
ENDFOR

; what to do if number of files is not the same?
; right now just get the minimum number of files

file_list = strarr(min(nfiles), nwave)

FOR i = 0, nwave-1 DO BEGIN
	IF fileset EQ 'ssw_cutout' THEN BEGIN
		files = file_search(dir + fileset + '*' + wave_list[i] + '_.fts')
		nfiles[i] = n_elements(files)
		times = anytim(aiacutout_to_time(files))
	ENDIF
	
	IF fileset EQ 'AIA' THEN BEGIN
		files = file_search(dir + fileset + '*0' + wave_list[i] + '.fits')
		nfiles[i] = n_elements(files)
		times = anytim(aiaprep_to_time(files))
	ENDIF
	
	s = sort(times)
	files = files[s]
	file_list[*,i] = files[0:min(nfiles)-1]
ENDFOR

RETURN, file_list

END