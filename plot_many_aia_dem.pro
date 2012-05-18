PRO plot_many_aia_dem, dir, istart = istart, hsi_fits = hsi_fits, hsi_image = hsi_image, OUTPLOT = outplot, DEBUG = debug, FLARE_NUM = flare_num

default, istart, 0
default, dir, '/Users/schriste/Desktop/flare2_deem/'
default, hsi_fits, '~/Dropbox/idl/aia_deem/ospex_results_8_oct_2011.fits'
default, hsi_image, '~/Dropbox/idl/aia_deem/hsi_image_20110826_205258.fits'

IF keyword_set(flare_num) THEN BEGIN
	dir = '~/Desktop/' + ['flare0_deem/', $
	'flare1_deem/', $
	'flare2_deem/', $
	'flare3_deem/']
	
	hsi_fits_dir = '/Users/schriste/Desktop/flare_' + num2str(flare_num-1) + '_rhessi/'
	f = file_search(hsi_fits_dir + 'ospex_*.fits')
	hsi_fits = f[0]
		
	f = file_search(hsi_fits_dir + 'hsi_image*.fits')
	hsi_image = f[0]
	
	;['hsi_image_20110716_170350.fits', $
	;'hsi_image_20110826_205258.fits', $
	;'hsi_image_20110603_071626.fits', $
	;'hsi_image_20110621_181902.fits']
	
	dir = dir[flare_num-1]
	print, dir
	print, hsi_fits
	print, hsi_image
ENDIF

f = file_search(dir + 'teem_tot*_mask_.sav')

dim = n_elements(f)

fit_result = dblarr(dim,3)

FOR i = istart, dim-1 DO BEGIN

    IF keyword_set(OUTPLOT) THEN BEGIN
        set_plot, 'z'
        loadct,0
        hsi_linecolors
        device, set_resolution = [800, 600]
	ENDIF

    plot_aia_dem, f[i], hsi_fits = hsi_fits, hsi_image = hsi_image, /bkg_file, fit = fit
    fit_result[i,*] = fit
    IF keyword_set(OUTPLOT) THEN BEGIN
		tvlct, r, g, b, /get
        outfile = 'dem_' + num2str(i, padchar = '0', length = 4) + '_' + break_time(anytim(fit[0],/yoh)) + '.png'
        write_png, outfile, tvrd(), r,g,b
        set_plot, 'x'
	ENDIF
    IF keyword_set(DEBUG) THEN stop
    
ENDFOR

save, fit_result, filename = 'aia_dem_slope_' + num2str(flare_num) + '.sav'

END