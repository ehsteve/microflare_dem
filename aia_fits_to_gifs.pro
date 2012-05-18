PRO aia_fits_to_gifs, dir = dir, OUTPLOT = outplot

dir = '/Volumes/sdo-aia/hsi_flare_20110716_/'
;fileset = file_search(dir + 'ssw_cutout_20110716_*')

fileset = 'ssw_cutout'
wave_ =['131','171','193','211','335','94'] 
nwave =n_elements(wave_) 
nfiles = fltarr(nwave)

FOR i = 0, nwave-1 DO BEGIN
	files = file_search(dir + fileset + '*' + wave_[i] + '_.fts')
	nfiles[i] = n_elements(files)
ENDFOR

;the number of files need to be the same!
print, nfiles
file_list = strarr(nfiles[0], nwave)

FOR i = 0, nwave-1 DO BEGIN
	files = file_search(dir + fileset + '*' + wave_[i] + '_.fts')
	nfiles[i] = n_elements(files)
	times = anytim(aiacutout_to_time(files))
	s = sort(times)
	files = files[s]
	file_list[*,i] = files
ENDFOR

FOR j = 0, nwave-1 DO BEGIN

	aia_lct, rr, gg, bb, wavelnth=wave_[j], /load 

	FOR i = 0, nfiles[0]-1 DO BEGIN

	curfile = file_list[i,j]

	read_sdo, curfile, header, data
	wcs = fitshead2wcs( header )
	wcs2map, data, wcs, map, id = header.wavelnth
	
	IF keyword_set(OUTPLOT) THEN BEGIN
        set_plot, 'z'  
		device, set_resolution = [600, 600]
	ENDIF
	
	cur_time = anytim(map.time)
	plot_map, map, charsize = 1.5
	legend, wave_[j], box = 0, charsize = 1.5

	IF keyword_set(OUTPLOT) THEN BEGIN
		tvlct, r, g, b, /get
        outfile = 'dem_' + wave_[j] + '_' + break_time(cur_time) + '_' + num2str(i, padchar = '0', length = 3) + '.png'
        write_png, outfile, tvrd(), r,g,b
        set_plot, 'x'
	ENDIF

ENDFOR
ENDFOR

END
