PRO temp_flareteem_info

dir = ['hsi_flare_20110716_170350', $
	'hsi_flare_20110826_205258', $
	'hsi_flare_20110603_071626', $
	'hsi_flare_20110621_182202']
	
hsi_fits = ['hsi_image_20110716_170350.fits', $
	'hsi_image_20110826_205258.fits', $
	'hsi_image_20110603_071626.fits', $
	'hsi_image_20110621_181902.fits']

FOR i = 0, n_elements(dir) DO BEGIN
	f = file_search('~/Desktop/' + dir[i] + '/*.fts')
	times = anytim(aiacutout_to_time(f))
	
	time_range = minmax(times)
	flares = hsi_whichflare(time_range)
	flare_str = hsi_whichflare(time_range, /only_one, /structure)
	flare_str = hsi_whichflare(time_range, /closest, /structure)
		
	print, 'start_time:', anytim(flare_str.start_time,/yoh)
	print, 'peak_time:', anytim(flare_str.peak_time,/yoh)
	print, 'end_time:', anytim(flare_str.end_time, /yoh)
	print, 'flare_id:', flare_str.id_number
	
	print, get_goes_class([flare_str.start_time, flare_str.end_time])
	
	stop
ENDFOR


END
