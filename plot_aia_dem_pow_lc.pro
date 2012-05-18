PRO plot_aia_dem_pow_lc, file, flare_num = flare_num



IF keyword_set(flare_num) THEN BEGIN
	dir = '~/Desktop/' + ['flare0_deem/', $
	'flare1_deem/', $
	'flare2_deem/', $
	'flare3_deem/']
	dir = dir[flare_num-1]
	f = file_search(dir + 'aia_dem_slope_*.sav')
	file = f[0]
	
	hsi_fits_dir = '/Users/schriste/Desktop/flare_' + num2str(flare_num-1) + '_rhessi/'
	f = file_search(hsi_fits_dir + 'hsi_spectrum*.fits')
	hsi_spectrum = f[0]
	
	f = file_search(hsi_fits_dir + 'hsi_srm*.fits')
	hsi_srm = f[0]
	
	
	;['hsi_image_20110716_170350.fits', $
	;'hsi_image_20110826_205258.fits', $
	;'hsi_image_20110603_071626.fits', $
	;'hsi_image_20110621_181902.fits']
	if flare_num eq 1 then timerange = ['16-Jul-11 17:00:11.000','16-Jul-11 17:20:00.000']
	;dir_aia = dir_aia[flare_num-1]
ENDIF
yrange = [50,20000]
print, file
restore, file, /verbose
!P.multi = [0,1,2]

if not exist(timerange) then timerange = anytim(minmax(fit_result[*,0]))
fit_yrange = minmax([fit_result[*,2],fit_result[*,1]] )

hsi_linecolors
colors = [7,8]
charsize = 1.5
utplot, anytim(fit_result[*,0],/yoh), fit_result[*,1], /nodata, charsize = charsize, timerange = timerange, ytitle = textoidl('\alpha_{fit}'), yrange = fit_yrange
outplot, anytim(fit_result[*,0],/yoh), fit_result[*,1], thick = 2.0, color = colors[0]
outplot, anytim(fit_result[*,0],/yoh), fit_result[*,2], thick = 2.0, color = colors[1]
legend, ['flare', 'background'], color = colors, charsize = charsize, linestyle = [0,0], thick = [2,2], /right

o = ospex(/no_gui)
o->set, spex_specfile = hsi_spectrum
o->set, spex_drmfile= hsi_srm
o->set, spex_eband=get_edge_products([4,15,25,50,100],/edges_2)

o -> plot_time, /data, /no_plotman, this_band = [0,1], charsize = charsize, timerange = timerange, yrange = yrange, /ystyle

stop

END