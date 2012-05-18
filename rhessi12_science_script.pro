
;key
;flare 0
;AIA files time ranges
;20110716_163301, 20110716_173259
t0 = ['2011/07/16 16:33:01', '2011/07/16 17:32:59']
;Flare 0: 11071604, B6.39, location = [-586, -437]
;16-jul-2011, 17:01, 17:03, 17:13
;istart = 117, 
;
;flare 1 - the limb
;AIA files time ranges
;20110826_202202, 20110826_212159
t1 = ['2011/08/26 20:22:02', '2011/08/26 21:21:59']
;Flare 1: 11082631, B4.16, location = [-908, -303]
;26-aug-2011, 20:51, 20:53, 20:58
;istart = 130
;
;flare 2
;AIA files time ranges
;20110621_175201, 20110621_185159
t2 = ['2011/06/21 17:52:01', '2011/06/21 18:51:59']
;Flare 3: 11062114, B2.83, location [284, 263]
;21-Jun-2011, 18:19, 18:22, 18:24
;istart = 113
;

;flare 3
;AIA files time ranges
;20110603_064601, 20110603_074559
t3 = ['2011/06/03 06:46:01', '2011/06/03 07:45:59']
;Flare 3: 11060303, B3.69, location [19.04, -362.25]
;3-Jun-2011, 07:15, 07:16, 07:19

;flare 4
;

; TO DO, now do ospex fit
; 12 second integration
; do one time interval during impulsive phase of flare to compare
; multi_therm and a therm and bpow

;need to create the dem for flare 3
;make image for flare 3
;plot light curves of slope for aia dem to compare with nonthermal HXR emission
;create light curves of "high temperature" filters of AIA to compare to rhessi

;how does the area changes as a function of temperature?

;to get goes lightcurves
lc = get_goes_lc(time_range)

;to create plot overlays
plot_aia_dem_hist, flare_num = 1, index = 117, wave = 5