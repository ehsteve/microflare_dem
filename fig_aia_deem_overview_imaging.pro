PRO fig_aia_deem_overview_imaging, FORCE = force, OUTPS = outps, filename = filename, image_algo = image_algo, image_wave = image_wave, energy_range = energy_range, _EXTRA = _EXTRA

default, filename, 'fig_aia_deem_overview_imaging'
default, charsize, 1.5
default, thick, 2.0
default, image_algo, 'clean'
default, image_wave, 5
default, energy_range, [4,10]

;example fig_aia_deem_overview_imaging, dmin = 6, /log

;vis fwdfit
;mem njit
;uv_smooth

;energy_range = [[4.0000000D, 10.000000D], [15.000000D, 30.000000D]]
wave = ['131','171','193','211','335','94']

!P.MULTI = [0,2,2]
;flare 0

image_wave = 5

tr0 = [['03-jun-2011 07:16', '03-jun-2011 07:17'], $
	   ['21-jun-2011 18:22', '21-jun-2011 18:23'], $
	   ['16-Jul-2011 17:02:01.424', '16-Jul-2011 17:02:56.116'], $
	   ['26-aug-2011 20:52', '26-aug-2011 20:53']]

dir_aia = '/Users/schriste/Desktop/' + ['hsi_flare_20110603_071626/', 'hsi_flare_20110621_182202/', 'hsi_flare_20110716_170350/', 'hsi_flare_20110826_205258/']

aia_image_index = [123,121,120,120]

IF keyword_set(OUTPS) THEN BEGIN
	charsize = 1.0
	popen, filename, xsize = 6, ysize = 6
ENDIF

FOR i = 0, 4-1 DO BEGIN
	
	tr0_str = break_time(tr0[0,i])
	hsi_image_fname = 'hsi_image_' + tr0_str + '_' + image_algo + '.fits'
	
	file_list = get_aia_file_list(dir_aia[i])
	aia_times = aiacutout_to_time(file_list[*,image_wave])
	
	aia_file = file_list[aia_image_index[i], image_wave]
	aia_bkg_file = file_list[0,image_wave]
	
	f = file_search(hsi_image_fname)
		
	IF f[0] EQ '' OR keyword_set(force) THEN BEGIN
		obj = hsi_image()                                                                         
		obj-> set, det_index_mask= [0B, 1B, 1B, 1B, 1B, 1B, 0B, 0B, 0B]                           
		obj-> set, im_energy_binning= energy_range[*,0]        
		obj-> set, im_time_interval= tr0[*,i] 
		;obj-> set, image_algorithm= 'VIS FWDFIT'  
		obj-> set, image_algorithm = image_algo
		obj-> set, clean_show_maps= 0                                                             
obj-> set, clean_progress_bar= 0   
		obj-> set, image_dim= [128, 128]                                                          
		obj-> set, pixel_size= [1.00000, 1.00000]                                                 
		obj-> set, time_bin_def= [1.00000, 2.00000, 4.00000, 4.00000, 8.00000, 16.0000, 32.0000, $
		 64.0000, 64.0000]                                                                        
		obj-> set, time_bin_min= 256L                                                             
		obj-> set, use_phz_stacker= 1L 
		obj-> set, vf_loop= 1
		obj-> set, vis_plotfit= 0B
		obj-> set, profile_show_plot= 0                                                           
		obj-> set, profile_plot_rate= 0                                                           
		obj-> set, profile_plot_resid= 0   
		data = obj-> getdata()
		obj->fitswrite, this_out_filename = hsi_image_fname
	ENDIF
	
	fits2map, hsi_image_fname, hsi_map
	fits2map, aia_file, aia_map
	fits2map, aia_bkg_file, aia_bkg_map
	
	aia_bkg_drot_map = drot_map(aia_bkg_map, ref_map = aia_map)
	
	r = get_map_center(hsi_map)
	dr = 40
	xrange = r[0] + dr*[-1,1]
	yrange = r[1] + dr*[-1,1]
	
	aia_lct,r,g,b,wavelnth=wave[5],/load
	tmp = [r[n_elements(r)-1], g[n_elements(r)-1], b[n_elements(r)-1]]
	r = reverse(r)
	g = reverse(g)
	b = reverse(b)
	r[n_elements(r)-1] = tmp[0]
	g[n_elements(r)-1] = tmp[1]
	b[n_elements(r)-1] = tmp[2]
	tvlct,r,g,b
	hsi_linecolors

	dmap = diff_map(aia_map, aia_bkg_drot_map)
	
	title = strmid(dmap.time, 0,20)
	
	plot_map, dmap, xrange = xrange, yrange = yrange, charsize = charsize, /limb, grid_spacing = 2, title = title, bottom = 13, _EXTRA = _EXTRA, dmin = 30
	
	text = '' + image_algo + ' ' + num2str(energy_range[0]) + '-' + num2str(energy_range[1]) + ' keV'
	
	plot_map, hsi_map, /over, /percent, levels = [50,60,70,80,90], color = 6, thick = thick
	legend, text, linestyle = 0, charsize = charsize, color = 6, thick = thick, pspacing = 1
	legend, wave[image_wave] + ' A ', box = 0, /bottom, /right, charsize = charsize
	
ENDFOR

IF keyword_set(OUTPS) THEN pclose

stop

END