FUNCTION aiaprep_to_time, f

;Convert a aia cutout filename to a time from the filename
;
;Written: Steven Christe (7-Feb-2011)

dim = n_elements(f)
result = strarr(dim)

FOR i = 0, dim-1 DO BEGIN

	filename = file_basename(f[i])
	
	year = strmid(filename,3,4)
	month = strmid(filename,7,2)
	day = strmid(filename,9,2)
	time = strmid(filename,12,2) + ':' + strmid(filename,14,2) + ':' + strmid(filename,16,2)

	result[i] = anytim(year + '/' + month + '/' + day + ' ' + time,/ecs)
ENDFOR

RETURN, result

END