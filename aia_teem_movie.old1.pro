PRO aia_teem_movie, DIR = dir, istart = istart, HSI_FITS = hsi_fits, DEBUG = debug, OUTPLOT = outplot, XRAYSPEC = xrayspec, SAVE_DIR = save_dir, flare_num = FLARE_NUM, VERBOSE = verbose, ISTOP = istop, fileset = fileset, fov = fov, npix = npix, FORCE_TABLE = force_table, q94 = q94

; 
;
;KEYWORDS: 
;			SAVE_DIR - 
;			FOV - must be given in arcsec [xmin,xmax,ymin,ymax]
;
;WRITTEN: Steven Christe (11-Oct-2011)
;REVISION: Steven Christe (6-Jan-2011)

default, istart, 0
default, dir, '~/Desktop/hsi_flare_20110716_170350/'
default, hsi_fits, '~/Dropbox/idl/aia_deem/hsi_image_20110716_170350.fits'
default, save_fit, ''
default, hsi_image, ''
default, fileset, 'ssw_cutout'
default, npix, 6 ;(macropixel size=4x4 pixels, yields 512x512 map) 
default, q94, 6.7 ;(correction factor for low-temperature 94 A response)

IF keyword_set(flare_num) THEN BEGIN
	;dir = '~/Desktop/' + ['hsi_flare_20110716_170350/', $
	;'hsi_flare_20110826_205258/', $
	;'hsi_flare_20110603_071626/', $
	;'hsi_flare_20110621_182202/']
	
	dir = '~/Desktop/' + [ $
	'hsi_flare_20110603_071626/', $
	'hsi_flare_20110621_182202/', $
	'hsi_flare_20110716_170350/', $
	'hsi_flare_20110826_205258/', $
	'hsi_flare_20101016/']
	
	;key 0->2 1->3, 2->0, 3->1
	
	save_dir = '~/Desktop/' + ['flare0_deem/', $
	'flare1_deem/', $
	'flare2_deem/', $
	'flare3_deem/', $
	'flare4_deem/']
	
	hsi_fits_dir = '/Users/schriste/Desktop/flare_' + num2str(flare_num-1) + '_rhessi/'
	f = file_search(hsi_fits_dir + 'ospex_*.fits')
	hsi_fits = f[0]
		
	f = file_search(hsi_fits_dir + 'hsi_image*.fits')
	hsi_image = f[0]
	
	dir = dir[flare_num-1]
	save_dir = save_dir[flare_num-1]
ENDIF

;note: flare happens at i = 120
loadct,0
hsi_linecolors

wave_ =['131','171','193','211','335','94'] 
nwave =n_elements(wave_) 
nfiles = fltarr(nwave)

file_list = get_aia_file_list(dir, fileset = fileset)

FOR i = 0, nwave-1 DO nfiles[i] = n_elements(file_list[*,i])

;t_min = 5.0
;t_max = 6.5
;t_d = 0.1
;telog1 = t_d * findgen((t_max - t_min)/t_d) + t_min

t_min = 5.5
t_max = 8.0

t_min = 5.5
t_max = 7.0
t_d = 0.05
telog = t_d * findgen((t_max - t_min)/t_d) + t_min

;telog = [telog1,telog2]

;te_range=[0.1,100]*1.e6 ;   ([K], valid temperature range for DEM solutions)

;tsig=0.1*(1+1*findgen(10)) ;   (values of Gaussian logarithmic temperature widths) 
tsig_min = 0.01
tsig_max = 1.0
tsig_d = 0.01
tsig = tsig_d * findgen((tsig_max - tsig_min)/tsig_d) + tsig_min

;not sure what the following does
vers='a' ;   (version number of label in filenames used for full images) 

; Savefile that contains DEM loopup table
teem_table='teem_table.sav' 
; Savefile that contains EM and Te maps
teem_map =fileset+vers+'_teem_map.sav' 
;  jpg-file that shows EM and Te maps
teem_jpeg=fileset+vers+'_teem_map.jpg' 
; only need to do this once

f = file_search(teem_table)
IF f[0] EQ '' OR keyword_set(FORCE_TABLE) THEN aia_teem_table,wave_,tsig,telog=telog,q94,filelist=file_list[0,*], teem_table, save_dir = save_dir

;aia_teem_table,wave_,tsig,te_range=te_range,q94,filelist=file_list[0,*], teem_table, save_dir = save_dir

IF NOT keyword_set(istop) THEN istop = nfiles[0]

IF hsi_image[0] NE '' THEN BEGIN
	fits2map, file_list[0,0], aiamap
	fits2map, hsi_image, map
	; interpolate the rhessi map to the aia map
	imap = inter_map(map,aiamap)
	m = max(imap.data)
	; set the mask at everything above 50% contour
	index = where(imap.data LE m*0.5, complement = complement)
	imap.data[index] = 0
	imap.data[complement] = 1
	mask = imap.data
	invmask = NOT mask
ENDIF


FOR i = istart, nfiles[0]-1 DO BEGIN
	;window,0
	IF fileset EQ 'ssw_cutout' THEN cur_time = anytim(aiacutout_to_time(file_list[i,0]),/yoh)
	IF fileset EQ 'AIA' THEN cur_time = anytim(aiaprep_to_time(file_list[i,0]), /yoh)
	print, cur_time	
	; Savefile that contains EM and Te maps
	teem_map ='teem_data_' + num2str(i, padchar = '0', length = 3) + break_time(cur_time) + '.sav' 
	
	teem_tot ='teem_tot_' + num2str(i, padchar = '0', length = 3) + break_time(cur_time) + '_q94_' + num2str(10*q94, length = 2) + '.sav'
			
	aia_teem_map,wave_,npix,teem_table,teem_map,filelist = file_list[i,*], filename_extra = '_' + num2str(i, padchar = '0', length = 3), save_dir = save_dir, debug = debug, verbose = verbose, fov = fov
	
	IF hsi_image[0] NE '' THEN BEGIN
		tx ='teem_tot_' + num2str(i, padchar = '0', length = 3) + break_time(cur_time)
		teem_tot = [tx + '_mask_.sav',tx + '_invmask_.sav']
		aia_teem_total,npix,wave_,q94,teem_table,teem_map,teem_tot[0],filelist = file_list[i,*], mask = mask, /plot, save_dir = save_dir, fov = fov
		aia_teem_total,npix,wave_,q94,teem_table,teem_map,teem_tot[1],filelist = file_list[i,*], mask = invmask, save_dir = save_dir, fov = fov
	ENDIF ELSE BEGIN
		aia_teem_total,npix,wave_,q94,teem_table,teem_map,teem_tot,filelist = file_list[i,*], fov = fov
	ENDELSE
	
	stop
	
	IF keyword_set(OUTPLOT) THEN BEGIN
        set_plot, 'z'
        loadct,0
        hsi_linecolors
        device, set_resolution = [800, 600]
	ENDIF

	IF keyword_set(HSI_FITS) THEN BEGIN

	restore, save_dir + teem_tot[0]
	telog_mask = telog
	emlog_mask = emlog
	emlog_err_mask = emlog_err
	restore, save_dir + teem_tot[1]

	hsi_linecolors
	yrange = [min([emlog_mask, emlog]), max([emlog_mask, emlog])]
	plot, telog, emlog, xtitle = 'log(Temperature [MK])', ytitle = 'log(Emmission Measure [cm!U-5!N])', /nodata, yrange = yrange, /ystyle,title = cur_time, charsize = 2.0
	
	oplot, telog, emlog, thick = 2.0, linestyle = 1, psym = symcat(16), color = 6
	oploterr, telog, emlog, emlog_err, /nohat, color = 6, errcolor = 6
	
	oplot, telog_mask, emlog_mask, thick = 2.0, psym = symcat(16), color = 7
	oploterr, telog_mask, emlog_mask, emlog_err_mask, /nohat, color = 7, errcolor = 7

	;get the maxima
	m = max(emlog, mindex)
	max_t = telog[mindex]
	
	m = max(emlog_mask, mindex)
	max_t_mask = telog_mask[mindex]
	
	;text = ['Flare T!Lmax!N = ' + mynum2str(10^max_t_mask),'Background T!Lmax!N = ' + mynum2str(10^max_t)]
	
	; now fit a line/power-law to the data above 5 MK
	index = where(telog GE alog10(5e6))
	x = telog_mask[index]
	y = emlog_mask[index]
	eq_mask = linfit(x,y)
		
	oplot, [alog10(5e6), alog10(1e8)], eq_mask[1]*[alog10(5e6), alog10(1e8)] + eq_mask[0], color = 5

	x = telog[index]
	y = emlog[index]
	eq_bkg = linfit(x,y)
	
	oplot, [alog10(5e6), alog10(1e8)], eq_bkg[1]*[alog10(5e6), alog10(1e8)] + eq_bkg[0], color = 6

	IF keyword_set(HSI_FITS) THEN BEGIN
		print, hsi_image
		print, hsi_fits
		fits2map, hsi_image, map
		c = 0.5
		frac = n_elements(where(map.data/max(map.data) GE c))/float(n_elements(map.data))
		area = n_elements(map.data)*frac*map.dx*map.dy * 712e5^2.0
	
		result = spex_read_fit_results(hsi_fits)
		p = result.spex_summ_params
		
		;now find the fit that is closest in time to time t
		hsi_time = result.spex_summ_time_interval
		;con1 = anytim(t) LE anytim(hsi_time[1,*])
		;con2 = anytim(t) GE anytim(hsi_time[0,*])
		
		params = p[*,0]
	
		;multi_therm_pow
		;DEM(T) = a[0] * (2/T)^a[3]
		;a[0] diff emission measure at T = 2 keV, 10^49 cm^(-3) keV^-1
		;a[1] min plasma temperature (keV)
		;a[2] max plasma temperature (keV)
		;a[3] power law index
		;a[4] relative abundance
		
		hsi_telog = alog10([params[1], params[2]]*11.6d6)
		hsi_emlog = alog10(1/11.6d6*1d49*params[0]*([2.0/params[1], 2.0/params[2]])^(params[3]))
				
		hsi_emlog = hsi_emlog - alog10(area)
		oplot, hsi_telog, hsi_emlog, thick = 2.0, color = 4   
		
		;legend, ['rhessi (' + textoidl('\alpha = ') + num2str(-params[3], length = 4) + ')', 'aia', 'aia (bkg)'], color = [4,1,1], linestyle = [0,1,2], pspacing = 1, charsize = 1.3
		;legend, ['rhessi (' + textoidl('\alpha = ') + num2str(-params[3], length = 4) + ')', color = [4,1,1], linestyle = [0,1,2], pspacing = 1, charsize = 1.3
		
		hsi_time_str = anytim(hsi_time, /yoh)
		t = strmid(hsi_time_str[0], 0,10) + strmid(hsi_time_str[0], 10,8) + ' to ' + strmid(hsi_time_str[1], 10,8)
		text = ['aia', 'aia (bkg)', 'rhessi (' + t + ')']
	
		legend, text, linestyle = [0,1, 0], /left, box = 1, /bottom, charsize = 1.5, pspacing = 1, color = [7,6,4], textcolor = [7,6,4]
	
		text = textoidl('\alpha_{fit} = ') + [num2str(eq_bkg[1], length = 4), num2str(eq_mask[1], length = 4)]
		text = [text, textoidl('\alpha_{hsi} = ') + num2str(-params[3], length = 4)]
		legend, text, linestyle = [0,0,0], color = [5,6,4], charsize = 1.5, /top, /right, pspacing = 1, thick = [1,1,2]
	
	ENDIF ELSE BEGIN
		text = ['aia', 'aia (bkg)']
	
		legend, text, linestyle = [0,1], /left, box = 1, /bottom, charsize = 1.5, pspacing = 1, color = [7,6], textcolor = [7,6]
		
		text = textoidl('\alpha_{fit} = ') + [num2str(eq_bkg[1], length = 4), num2str(eq_mask[1], length = 4)]
		legend, text, linestyle = [0,0], color = [5,6], charsize = 1.5, /top, /right, pspacing = 1
	
			;legend, ['rhessi', 'aia', 'aia (bkg)'], color = [4,1,1], linestyle = [0,1,2], pspacing = 1, charsize = 1.3
	ENDELSE 
			
	IF keyword_set(OUTPLOT) THEN BEGIN
		tvlct, r, g, b, /get
        outfile = save_dir + 'dem_' + num2str(i, length = 3, padchar = '0') + '_' + break_time(cur_time) + '_' + num2str(i, padchar = '0', length = 3) + '.png'
        write_png, outfile, tvrd(), r,g,b
        set_plot, 'x'
	ENDIF
	
ENDIF

	IF keyword_set(XRAYSPEC) THEN BEGIN
		; now simulate the x-ray spectrum for this area
			dem = [transpose(telog), transpose(emlog)]
			result = chianti_spec_from_dem(dem, /plot)
		IF keyword_set(hsi_fits) THEN BEGIN	
			dem_mask = [transpose(telog_mask), transpose(emlog_mask)]
			result_mask = chianti_spec_from_dem(dem_mask)
			oplot, result_mask[0,*], result_mask[1,*]
		ENDIF
	ENDIF
		
	IF keyword_set(debug) THEN stop
	
ENDFOR

END