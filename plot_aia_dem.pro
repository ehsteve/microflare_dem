PRO plot_aia_dem, file, HSI_FITS = hsi_fits, HSI_IMAGE = hsi_image, BKG_FILE = bkg_file, FIT = fit

default, file, '/Users/schriste/Desktop/flare2_deem/teem_tot_00020110603_064611_mask_.sav'
default, hsi_fits, '~/Dropbox/idl/aia_deem/ospex_results_8_oct_2011.fits'
default, hsi_image, '~/Dropbox/idl/aia_deem/hsi_image_20110826_205258.fits'
    
hsi_linecolors

;parse the time out of the filename
f = filename(file)
t = anytim(unbreak_time(strmid(f,12,15)),/yoh)

IF keyword_set(BKG_FILE) THEN BEGIN
    bkg_file = str_replace(file, 'mask', 'invmask')
    f = file_search(bkg_file)
    IF f[0] EQ '' THEN stop
    restore, bkg_file
    
    emlog_bkg = emlog
    telog_bkg = telog
ENDIF

f = file_search(bkg_file)
IF f[0] EQ '' THEN stop
restore, file
    
yrange = minmax(emlog)

title = t

plot, telog, emlog, xtitle = 'log(Temperature [K])', ytitle = 'Emission measure  log(EM [cm!U-5!N K!U-1!N])', /nodata, yrange = yrange, /ystyle, charsize = 1.5, title = title
oplot, telog, emlog, thick = 2.0, linestyle = 1

IF keyword_set(BKG_FILE) THEN oplot, telog_bkg, emlog_bkg, thick = 2.0, linestyle = 2

;text = ['Flare T!Lmax!N = ' + mynum2str(10^max_t_mask),'Background T!Lmax!N = ' + mynum2str(10^max_t)]

;legend, text, linestyle = [0,1], /right, box = 0, charsize = 1.5, /bottom

; now fit a line/power-law to the data above 5 MK
index = where(telog GE alog10(5e6))

x = telog[index]
y = emlog[index]
res0 = linfit(x,y)

oplot, [alog10(5e6), alog10(1e8)], res0[1]*[alog10(5e6), alog10(1e8)] + res0[0], color = 6

x = telog_bkg[index]
y = emlog_bkg[index]
res1 = linfit(x,y)

fit = [anytim(t), res0[1], res1[1]]

oplot, [alog10(5e6), alog10(1e8)], res1[1]*[alog10(5e6), alog10(1e8)] + res1[0], color = 7

text = [textoidl('\alpha_{fit} = ') + num2str(res0[1], length = 4), textoidl('\alpha_{fit} = ') + num2str(res1[1], length = 4)]
legend, text, linestyle = 0, color = [6,7], charsize = 1.5, /bottom, pspacing = 1

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
    con1 = anytim(t) LE anytim(hsi_time[1,*])
    con2 = anytim(t) GE anytim(hsi_time[0,*])
    
    IF con1 AND con2 THEN BEGIN
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
        
        legend, ['rhessi (' + textoidl('\alpha = ') + num2str(-params[3], length = 4) + ')', 'aia', 'aia (bkg)'], color = [4,1,1], linestyle = [0,1,2], pspacing = 1, charsize = 1.3
    ENDIF ELSE BEGIN
        legend, ['rhessi', 'aia', 'aia (bkg)'], color = [4,1,1], linestyle = [0,1,2], pspacing = 1, charsize = 1.3
    ENDELSE 
    
ENDIF

END