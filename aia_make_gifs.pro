PRO aia_make_gifs, dir = dir

dir = '/Volumes/sdo-aia/hsi_flare_20110716_/'
;fileset = file_search(dir + 'ssw_cutout_20110716_*')

hsi_fits = '/Users/schriste/idlpro/schriste/aia_deem/hsi_image_20110716_170350.fits'

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

FOR i = 0, nfiles DO BEGIN

	FOR j = 0, nwave-1 DO BEGIN
		

	fits2map, curfile, map
	
	IF keyword_set(OUTPLOT) THEN BEGIN
        set_plot, 'z'
        loadct,0
        hsi_linecolors
        device, set_resolution = [800, 600]
	ENDIF
	
	plot_map, map, charsize = 1.5, grid = 10

	IF keyword_set(OUTPLOT) THEN BEGIN
		tvlct, r, g, b, /get
        outfile = 'dem_' + break_time(cur_time) + '_' + num2str(i, padchar = '0', length = 3) + '.png'
        write_png, outfile, tvrd(), r,g,b
        set_plot, 'x'
	ENDIF


ENDFOR

END
