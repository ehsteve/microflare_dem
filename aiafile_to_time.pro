FUNCTION aiafile_to_time, f, fileset = fileset

;Convert a aia cutout filename to a time from the filename
;
;Written: Steven Christe (7-Feb-2011)

dim = n_elements(f)
result = strarr(dim)
default, fileset, 'aia'

FOR i = 0, dim-1 DO BEGIN

	filename = file_basename(f[i])
	
	IF fileset EQ 'aia' THEN BEGIN
		year = strmid(filename,3,4)
		month = strmid(filename,7,2)
		day = strmid(filename,9,2)
		time = strmid(filename,12,2) + ':' + strmid(filename,14,2) + ':' + strmid(filename,16,2)
	ENDIF ELSE IF fileset EQ 'ssw_cutout' THEN BEGIN
		year = strmid(filename,11,4)
		month = strmid(filename,15,2)
		day = strmid(filename,17,2)
		time = strmid(filename,20,2) + ':' + strmid(filename,22,2) + ':' + strmid(filename,24,2)
	ENDIF

	result[i] = anytim(year + '/' + month + '/' + day + ' ' + time,/ecs)
ENDFOR

RETURN, result

END