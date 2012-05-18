PRO fig_rhessi_lc_overview, OUTPS = outps

IF keyword_set(OUTPS) THEN BEGIN
	charsize = 1.0
	popen, filename, xsize = 8, ysize = 5
ENDIF

!P.MULTI = [0,2,2]

;flare 3
;AIA files time ranges
;20110603_064601, 20110603_074559
t3 = ['2011/06/03 07:00:01', '2011/06/03 07:30:59']
;Flare 3: 11060303, B3.69, location [19.04, -362.25]
;3-Jun-2011, 07:15, 07:16, 07:19

lc = get_rhessi_lc(t3, /corr, /plot)

;flare 2
;AIA files time ranges
;20110621_175201, 20110621_185159
t2 = ['2011/06/21 18:10:01', '2011/06/21 18:41:59']
;Flare 3: 11062114, B2.83, location [284, 263]
;21-Jun-2011, 18:19, 18:22, 18:24
;istart = 
;

lc = get_rhessi_lc(t2, /corr, /plot)

;key
;flare 0
;AIA files time ranges
;20110716_163301, 20110716_173259
t0 = ['2011/07/16 16:50:01', '2011/07/16 17:32:59']
;Flare 0: 11071604, B6.39, location = [-586, -437]
;16-jul-2011, 17:01, 17:03, 17:13
;istart = 117, 

lc = get_rhessi_lc(t0, /corr, /plot)

;flare 1 - the limb
;AIA files time ranges
;20110826_202202, 20110826_212159
t1 = ['2011/08/26 20:35:02', '2011/08/26 21:21:59']
;Flare 1: 11082631, B4.16, location = [-908, -303]
;26-aug-2011, 20:51, 20:53, 20:58
;istart = 130

lc = get_rhessi_lc(t1, /corr, /plot)



IF keyword_set(OUTPS) THEN pclose


END