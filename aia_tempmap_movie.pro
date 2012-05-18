PRO aia_tempmap_movie, dir, FLARE_NUM = flare_num, index = index

IF keyword_set(flare_num) THEN BEGIN
	dir = '~/Desktop/' + ['flare0_deem/', $
	'flare1_deem/', $
	'flare2_deem/', $
	'flare3_deem/']
	dir = dir[flare_num-1]
	
	dir_aia = '~/Desktop/' + ['hsi_flare_20110716_170350/', $
	'hsi_flare_20110826_205258/', $
	'hsi_flare_20110603_071626/', $
	'hsi_flare_20110621_182202/']
	
	
	hsi_fits_dir = '/Users/schriste/Desktop/flare_' + num2str(flare_num-1) + '_rhessi/'
	f = file_search(hsi_fits_dir + 'ospex_*.fits')
	hsi_fits = f[0]
		
	f = file_search(hsi_fits_dir + 'hsi_image*.fits')
	hsi_image = f[0]
	
	;['hsi_image_20110716_170350.fits', $
	;'hsi_image_20110826_205258.fits', $
	;'hsi_image_20110603_071626.fits', $
	;'hsi_image_20110621_181902.fits']
	dir_aia = dir_aia[flare_num-1]
	
	print, dir
	print, hsi_fits
	print, hsi_image
ENDIF

f = file_search(dir + 'teem_data_*.sav')

ilist = strmid(filename(f), 10, 3)
time_list = unbreak_time(strmid(filename(f),13,15))

restore, f[index], /verbose
fits2map, hsi_image, rmap

min_logtemperature = 5.5
max_logtemperature = 7.0
nbins = 20

!P.MULTI = [0,2,2]

plot_map, temperature_map, /cbar, dmin = min_logtemperature, dmax = max_logtemperature
plot_map, rmap, /over, /percent, levels = [50,60,70,80,90]

imap = inter_map(rmap,temperature_map)
m = max(imap.data)
; set the mask at everything above 50% contour
l = where(imap.data LE m*0.5, complement = complement)
temperature_map.data[l] = 0
yrange_t = minmax(temperature_map.data[complement])
EMISSION_MAP.data[l] = 0
;yrange_em = minmax(EMISSION_MAP.data[complement])
yrange_em = [21,22.5]
SIGMA_MAP.data[l] = 0
;yrange_sig = minmax(SIGMA_MAP.data[complement])
yrange_sig = [0,1]

hist = histogram(temperature_map.data, min = min_logtemperature, max = max_logtemperature, locations = loc, nbins = nbins)

plot, loc, hist, psym = 10, ytitle = 'Number of pixels', xtitle = 'log(Temperature [K])', /nodata, title = time_list[index] 
oplot, loc, hist, thick = 2.0, psym = 10

hist = histogram(EMISSION_MAP.data, min = yrange_em[0], max = yrange_em[1], locations = loc, nbins = nbins)

plot, loc, hist, psym = 10, ytitle = 'Number of pixels', xtitle = 'log(Emission measure [cm!U-5!N])', /nodata, title = time_list[index] 
oplot, loc, hist, thick = 2.0, psym = 10

hist = histogram(SIGMA_MAP.data, min = yrange_sig[0], max = yrange_sig[1], locations = loc, nbins = nbins)

plot, loc, hist, psym = 10, ytitle = 'Number of pixels', xtitle = 'log(sigma [K])', /nodata, title = time_list[index], /ylog, yrange = [1,max(hist)]
oplot, loc, hist, thick = 2.0, psym = 10


stop

END