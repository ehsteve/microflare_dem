FUNCTION break_time, t_in, qdebug = qdebug

; IDL Version 5.4 (sunos sparc)
; Journal File for hhudson@delsol
; Working directory: /home/hhudson/fe_line_ii
; Date: Sun Jul 20 11:45:38 2003
;
;MODIFIED: Steven Christe (10-apr-09), made compatible with arrays of t_in 

if n_elements(t_in) eq 0 then t_in = !stime

dim = n_elements(t_in)
filname = strarr(dim)

FOR i = 0, dim-1 DO BEGIN

    extime = anytim(t_in[i],/ext)
    year = strtrim(string(extime(6)),1)
    month = strtrim(string(extime(5)),1)

    if strlen(month) eq 1 then month = '0'+month
    day = strtrim(string(extime(4)),1)
    if strlen(day) eq 1 then day = '0'+day
    hour = strtrim(string(extime(0)),1)
    if strlen(hour) eq 1 then hour = '0'+hour
    minute = strtrim(string(extime(1)),1)
    if strlen(minute) eq 1 then minute = '0'+minute
    second = strtrim(string(extime(2)),1)
    if strlen(second) eq 1 then second = '0'+second
    filname[i] = year+month+day+'_'+hour+minute+second
    if keyword_set(qdebug) then stop

ENDFOR

RETURN, filname

end