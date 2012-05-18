PRO plot_aia_lightcurve, dir, FLARE_NUM = flare_num, index = index

IF keyword_set(flare_num) THEN BEGIN
	dir = '~/Desktop/' + ['flare0_deem/', $
	'flare1_deem/', $
	'flare2_deem/', $
	'flare3_deem/']
	dir = dir[flare_num-1]
	
	hsi_fits_dir = '/Users/schriste/Desktop/flare_' + num2str(flare_num-1) + '_rhessi/'
	f = file_search(hsi_fits_dir + 'ospex_*.fits')
	hsi_fits = f[0]
		
	f = file_search(hsi_fits_dir + 'hsi_image*.fits')
	hsi_image = f[0]
	
	hsi_fits_dir = '/Users/schriste/Desktop/flare_' + num2str(flare_num-1) + '_rhessi/'
	f = file_search(hsi_fits_dir + 'hsi_spectrum*.fits')
	hsi_spectrum = f[0]
	
	f = file_search(hsi_fits_dir + 'hsi_srm*.fits')
	hsi_srm = f[0]
	
	;['hsi_image_20110716_170350.fits', $
	;'hsi_image_20110826_205258.fits', $
	;'hsi_image_20110603_071626.fits', $
	;'hsi_image_20110621_181902.fits']
	
	print, dir
	print, hsi_fits
	print, hsi_image
ENDIF

f = file_search(dir + 'teem_data_*.sav')

ilist = strmid(filename(f), 10, 3)
time_list = unbreak_time(strmid(filename(f),13,15))

dim = n_elements(f)
!P.MULTI = [0,1,2]

t = strarr(dim)
aia_lc = fltarr(5,dim)
dem_lc = fltarr(4,dim)

restore, f[0], /verbose
fits2map, hsi_image, rmap
imap = inter_map(rmap,temperature_map)
m = max(imap.data)
index = where(imap.data LE m*0.5, complement = complement)

FOR i = 0, dim-1 DO BEGIN
    restore, f[i], /verbose
    temperature_map.data[index] = 0

    t[i] = temperature_map.time
    dem_lc[0,i] = average(temperature_map.data[complement])
    dem_lc[1,i] = average(emission_map.data[complement])
    dem_lc[2,i] = average(sigma_map.data[complement])
    dem_lc[3,i] = average(chisq_map.data[complement])
    
    FOR j = 0, 5-1 DO aia_lc[j,i] = average(aia_map_cube[j].data[complement])
    
ENDFOR

charsize = 2.0
!P.multi = [0,1,5]

o = ospex(/no_gui)
o->set, spex_specfile = hsi_spectrum
o->set, spex_drmfile= hsi_srm
o->set, spex_eband=get_edge_products([4,15,25,50,100],/edges_2)

o -> plot_time, /data, /no_plotman, this_band = [0,1], charsize = charsize, timerange = timerange, yrange = yrange, /ystyle


utplot, t, dem_lc[0,*], yrange = minmax(dem_lc[0,*]), ytitle = "log(Temperature [K])", /nodata, charsize = charsize
outplot, t, dem_lc[0,*], thick = 2
utplot, t, dem_lc[1,*], yrange = minmax(dem_lc[1,*]), ytitle = "log(Emission Measure [cm!U-5!N])", charsize = charsize, /nodata
outplot, t, dem_lc[1,*], thick = 2
utplot, t, dem_lc[2,*], yrange = minmax(dem_lc[2,*]), ytitle = "log(sigma [K])", charsize = charsize, /nodata
outplot, t, dem_lc[2,*], thick = 2
utplot, t, dem_lc[3,*], yrange = minmax(dem_lc[3,*]), ytitle = textoidl('\chi^2'), charsize = charsize, /nodata
outplot, t, dem_lc[3,*], thick = 2
stop

!P.multi = [0,1,5]
FOR i = 0, 5-1 DO BEGIN
    utplot, t, aia_lc[i,*], yrange = minmax(aia_lc[i,*]), ytitle = aia_map_cube[i].id, /nodata, charsize = charsize
    outplot, t, aia_lc[i,*], thick = 2
ENDFOR

stop

END