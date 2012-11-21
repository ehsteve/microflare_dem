PRO 

aia_dir = '/Users/schriste/Desktop/foxsi_aia_data/'
dir = '/Users/schriste/Dropbox/idl/schriste/foxsi_obs/'

file_list = file_search(aia_dir + 'ssw_cutout_20121102_18034*')

hsi_image = dir + 'hsi_image_20121102_175924.fits

fits2map, hsi_image, rmap

fileset = 'ssw_cutout'
fit_time = ['2012/11/02 18:00:00', '2012/11/02 18:05:00']
bkg_time = fit_time

aia_hsi_dem_analysis, fit_time = fit_time, bkg_time = bkg_time, dir = aia_dir, FILESET = fileset, /force_table, /aia_only, /epstein, hsi_image = hsi_image

filename = 'aia_fit_results_epstein.sav'
aia_hsi_dem_results, filename







aia_hsi_dem_analysis, dir = aia_dir, HSI_IMAGE = hsi_image, FILESET = fileset, /force_table, /aia_only, /epstein

!P.multi = [0,2,3]

xrange = [850,1050]
yrange = [-300,-100]

ps, 'aia_rhessi.ps', /color, /portrait

wavelength = ['94', '131', '171', '193', '211', '304', '335']
dmin = [50,50,50,50,50,50,50]
dmax = [400,400,3000,5000,4000,400]

FOR i = 0, n_elements(file_list)-1 DO BEGIN

	fits2map, file_list[i], aia_map
	fname = file_basename(file_list[i])
	w = str_replace(strmid(fname,7,3, /reverse_offset), '_', '')
	
	index = where(wavelength EQ w)
	print, w[index]
	aia_lct, r, g, b, wavelnth = w, /load
	hsi_linecolors

	plot_map, aia_map, bottom = 13, xrange = xrange, yrange = yrange, title = strmid(aia_map.time,0,20), dmin = dmin[index], dmax = dmax[index], grid = 10, /limb, charsize = 1.5
	plot_map, rmap, /over, levels = [50,60,70,80,90], /percent, color = 6
	legend, w, textcolor = 1
ENDFOR

psclose

END