pro aia_deem_plots, aia_deem_map_save = aia_deem_map_save, outps = outps,filename=filename

;+
; Project     : AIA/SDO
;
; Name        : AIA_TEMPMAP 
;
; Category    : Data analysis   
;
; Explanation : calculates temperature map
;		based on single-Gaussian fit
;
; Syntax      : IDL>aia_teem_plots,fileset,fov,wave_,teem_map
;
; Inputs      : fileset  = strarr(6) with filenames of 6 wavelength FITS images
;		fov(4)   = [i1,j1,i2,j2] pixel ranges of field-of-view
;               npix     = macropixel size 
;		wave_	 = strarr(6) with wavelengths in Angstroem
;       q94      = empirical correction factor of 94 A response
;       teem_table = savefile containing DEM lookup table
;       teem_map = savefile containing EM and Te maps
;		mask_map 	= given a mask as a map, return the dem only where it equals 1.
;
; Outputs     : teem_tot = savefile containing total DEM and fluxes
;
; History     :  9-Mar-2011, Version 1 written by Markus J. Aschwanden
; 			  :  10-May-2011, Version 2 added mask keyword by S. Christe
;				:1-Feb-2012 Version 2 - fixed outstanding bugs in code - A. Inglis
;				a) Lines 44-50 did not handle the AIA2011xxxx file format and have been replaced by setting 
;				file_iw=filelist directly (filelist is now in correct handling order already).
;				b)Error calculation was missing mask_frac on the denominator when error 8 is calculated (line 113)
;				c)mask_frac is now calculated correctly! See line 79.
;		:3-Feb-2012  fixed energy calculation, and now generates an energy_map with the energy in each pixel. Saved in teem_tot.....sav file
;
; Contact     : aschwanden@lmsal.com
;-
;

;restore, '/Users/schriste/Dropbox/idl/aia_deem/teem_data_00020110215_000009_npix10.sav'
IF keyword_set(filename) THEN BEGIN
restore,filename
ENDIF ELSE BEGIN
restore,'teem_data_00020110215_000009_55_80_q10_chiantifix.sav'
ENDELSE

IF keyword_set(OUTPS) THEN BEGIN
    loadct, 0
    hsi_linecolors
    popen, 'aia_teem_plots', xsize = 7, ysize = 6
ENDIF

nwave = n_elements(aia_map_cube)
wave =['131','171','193','211','335','94'] 

hist_color = 4

temp_title = 'log[temperature (K)]'
em_title = 'log[emission measure]'
sig_title = 'sigma'
chisq_title = 'reduced chisq'
qflux_title = 'qflux'

chisq = chisq_map.data
temp = temperature_map.data
em = EMISSION_MAP.data
sig = sigma_map.data

temp_range = [minnotzero(temp), max(temp)]
em_range = [minnotzero(em), max(em)]
sig_range = [minnotzero(sig), max(sig)]
chisq_range = [0, 5]
em_bins = 0.01
temp_bins = 0.05
sig_bins = 0.01
chisq_bins = 0.01

aia_qflux_cube = aia_map_cube
FOR i = 0, nwave-1 DO aia_qflux_cube[i].data = aia_simul_map_cube[i].data/float(aia_map_cube[i].data)

!P.MULTI = [0,1,2, 0, 0]
!X.MARGIN=[2,2];[0,0]
!Y.MARGIN=[2,2];[0,0]
FOR i = 0, nwave-1 DO BEGIN
    aia_lct, rr, gg, bb, wavelnth=wave[i], /load
    plot_map, aia_map_cube[i], /limb, /log, grid = 10, gcolor = 255,charsize=0.7
    legend, wave[i] + ' data', textcolor = 255, /right
    plot_map, aia_simul_map_cube[i], /limb, /log, grid = 10, title = '', gcolor = 255,charsize=0.7
    legend, wave[i] + ' simul', textcolor = 255, /right
ENDFOR
!P.MULTI = 0

loadct, 5
plot_map, temperature_map, dmin = 6, /limb, /isotropic,cbar=1
loadct, 1
plot_map, emission_map, dmin = 20, /limb, /isotropic,cbar=1
plot_map, sigma_map, /isotropic,cbar=1
plot_map, chisq_map, dmin = 0, /isotropic,cbar=1

!P.MULTI = [0,3,2]
!X.MARGIN=[8,3];[0,0]
!Y.MARGIN=[1,1];[0,0]
FOR i = 0, nwave-1 DO BEGIN
    aia_lct, rr, gg, bb, wavelnth=wave[i], /load
    plot_map, aia_map_cube[i], /limb, /log, title = ''
    legend, wave[i] + ' temp', textcolor = 255, /right
    plot_map, temperature_map, /over, levels = findgen(10)*0.1 + 6, color = 255
ENDFOR

!P.MULTI = [0,3,2]
FOR i = 0, nwave-1 DO BEGIN
    aia_lct, rr, gg, bb, wavelnth=wave[i], /load
    plot_map, aia_map_cube[i], /limb, /log, title = ''
    legend, wave[i] + ' em', textcolor = 255, box = 0, /right
    plot_map, emission_map, /over, levels = findgen(10)*0.5 + 20, color = 255
ENDFOR

!P.MULTI = [0,3,2]
FOR i = 0, nwave-1 DO BEGIN
    aia_lct, rr, gg, bb, wavelnth=wave[i], /load
    plot_map, aia_map_cube[i], /limb, /log
    legend, wave[i] + ' sig', textcolor = 255, box = 0, /right
    plot_map, sigma_map, /over, levels = findgen(30)*0.1 , color = 255
ENDFOR

!P.MULTI = 0
!X.MARGIN=[6,3];[0,0]
!Y.MARGIN=[4,2];[0,0]
loadct, 0
my_linecolors

plothist, temp, bin = temp_bins, /fill, xtitle = temp_title, xrange = temp_range, fcolor = hist_color
statistic, temp, xav, sigma, median, sig1, sig2
index = where(temp EQ 0, zero_count)
leg = ['#0s = ' + num2str(zero_count), 'mean = ' + num2str(xav), 'sig = ' + num2str(sigma)]
legend, leg , box = 0

plothist, em, bin = em_bins, /fill, xtitle = em_title, xrange = em_range, fcolor = hist_color
legend, '# of zeros = ' + num2str(zero_count), box = 0
statistic, em, xav, sigma, median, sig1, sig2
index = where(em EQ 0, zero_count)
leg = ['#0s = ' + num2str(zero_count), 'mean = ' + num2str(xav), 'sig = ' + num2str(sigma)]
legend, leg , box = 0

plothist, sig, bin = sig_bins, /fill, xtitle = sig_title, xrange = sig_range, fcolor = hist_color
statistic, sig, xav, sigma, median, sig1, sig2
index = where(sig EQ 0, zero_count)
leg = ['#0s = ' + num2str(zero_count), 'mean = ' + num2str(xav), 'sig = ' + num2str(sigma)]
legend, leg , box = 0

plothist, chisq, bin = chisq_bins, /fill, xtitle = chisq_title, fcolor = hist_color
statistic, chisq, xav, sigma, median, sig1, sig2
index = where(chisq EQ 0, zero_count)
leg = ['#0s = ' + num2str(zero_count), 'mean = ' + num2str(xav), 'sig = ' + num2str(sigma)]
legend, leg , box = 0

;FOR i = 0, nwave-1 DO plothist, aia_qflux_cube[i].data, bin = 0.1, /fill, xtitle = aia_qflux_cube[i].id


plot, chisq, temp, psym = 4, xtitle = chisq_title, ytitle = temp_title, $
	yrange = temp_range

plot, chisq, em, psym = 4, xtitle = chisq_title, ytitle = em_title, $
	yrange = em_range
	
plot, chisq, sig, psym = 4, xtitle = chisq_title, ytitle = sig_title, $
	yrange = sig_range

plot, em, temp, psym = 4, xtitle = em_title, ytitle = temp_title, $
	yrange = temp_range, xrange = em_range

plot, sig, em, psym = 4, xtitle = sig_title, ytitle = em_title, $
	yrange = em_range

plot, sig, temp, psym = 4, xtitle = sig_title, ytitle = temp_title, $
	yrange = temp_range

loadct, 1

img = hist_2d( chisq, temp, min1 = chisq_range[0], max1 = chisq_range[1], min2 = temp_range[0], max2 = temp_range[1], bin1 = chisq_bins, bin2 = temp_bins)
plot_image, img, origin = [chisq_range[0], temp_range[0]], scale = [chisq_bins,temp_bins], xtitle = chisq_title, ytitle = temp_title

img = hist_2d( em, temp, min1 = em_range[0], max1 = em_range[1], min2 = temp_range[0], max2 = temp_range[1], bin1 = em_bins, bin2 = temp_bins)
plot_image, img, origin = [em_range[0], temp_range[0]], scale = [em_bins,temp_bins], xtitle = em_title, ytitle = temp_title

img = hist_2d( sig, em, min1 = sig_range[0], max1 = sig_range[1], min2 = em_range[0], max2 = em_range[1], bin1 = sig_bins, bin2 = em_bins)
plot_image, img, origin = [sig_range[0], em_range[0]], scale = [em_bins,temp_bins], xtitle = em_title, ytitle = temp_title

img = hist_2d( sig, temp, min1 = sig_range[0], max1 = sig_range[1], min2 = temp_range[0], max2 = temp_range[1], bin1 = sig_bins, bin2 = temp_bins)
plot_image, img, origin = [sig_range[0], temp_range[0]], scale = [sig_bins,temp_bins], xtitle = sig_title, ytitle = temp_title


!P.multi = [0,2,nwave/2]

FOR i = 0, nwave-1 DO BEGIN 
    plothist, aia_qflux_cube[i].data < 10 > 0, bin = 0.1, xrange = [0, 3], /fill, xtitle = 'simul/data', fcolor = hist_color
    legend, wave[i], /right
ENDFOR

FOR i = 0, nwave-1 DO BEGIN 
	plot, aia_qflux_cube[i].data, temp, ytitle = temp_title, yrange = temp_range, psym = 3, $
	xtitle = qflux_title
	oplot, [1.0, 1.0], [-1d19, 1e19], linestyle = 2
	legend, wave[i], /right
ENDFOR

FOR i = 0, nwave-1 DO BEGIN
	plot, aia_qflux_cube[i].data, em, ytitle = em_title, yrange = em_range, psym = 3, $
	xtitle = qflux_title
	oplot, [1.0, 1.0], [-1d19, 1e19], linestyle = 2
	legend, wave[i], /right
ENDFOR
FOR i = 0, nwave-1 DO BEGIN
	plot, aia_qflux_cube[i].data, sig, ytitle = sig_title, yrange = sig_range, psym = 3, $
	xtitle = qflux_title
	oplot, [1.0, 1.0], [-1d19, 1e19], linestyle = 2
	legend, wave[i], /right
ENDFOR
IF keyword_set(OUTPS) THEN pclose

END