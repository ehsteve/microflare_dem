FUNCTION aiacutout_to_time, f

;Convert a aia cutout filename to a time from the filename
;
;Written: Steven Christe (7-Feb-2011)

dim = n_elements(f)
result = strarr(dim)

FOR i = 0, dim-1 DO BEGIN

	filename = file_basename(f[i])
	
	year = strmid(filename,11,4)
	month = strmid(filename,15,2)
	day = strmid(filename,17,2)
	time = strmid(filename,20,2) + ':' + strmid(filename,22,2) + ':' + strmid(filename,24,2)

	result[i] = anytim(year + '/' + month + '/' + day + ' ' + time,/ecs)
ENDFOR

RETURN, result

END