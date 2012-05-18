PRO aia_teem_stats, filename, OUTPLOT = outplot, outfile = outfile, totfile = totfile

default, filename, '/Users/schriste/Dropbox/idl/aia_deem/tempmap_000.sav'
default, wave_, ['131','171','193','211','335','94'] 
default, charsize, 1.5
default, outfile, 'aia_teem_stats'
default, totfile, '/Users/schriste/Dropbox/idl/aia_deem/teem_tot_00020110803_182309.sav'
nwave =n_elements(wave_) 

;dir = '/Users/schriste/data/sdo/aia/ar11263/'
;fileset = 'AIA'
;fov = [-198.857, 200.000, 50.0, 350.000]

;aia_teem_movie, dir = dir, fileset = fileset, fov = fov, npix = 4

restore, filename

IF keyword_set(OUTPLOT) THEN popen, outfile, xsize = 4, ysize = 6 

!P.MULTI = [0,2,6]
FOR i = 0, nwave-1 DO BEGIN
	aia_lct, rr, gg, bb, wavelnth=wave_[i], /load
	plot_map, aia_map_cube[i]
	dmap = diff_map(aia_map_cube[i], aia_simul_map_cube[i])
	plot_map, dmap, title = 'obs - simul'
ENDFOR

!P.MULTI = 0

loadct, 5
plot_map, temperature_map, charsize = charsize, /cbar
loadct, 0
plot_map, emission_map, charsize = charsize, /cbar
plot_map, sigma_map, charsize = charsize, /cbar

loadct, 0
plot_map, chisq_map, dmin = 0, dmax = 10, /cbar

loadct, 0
plot_map, TSIG_ERROR_MAP_MINUS, charsize = 1.0
plot_map, TSIG_ERROR_MAP_PLUS, charsize = 1.0
plot_map, TELOG_ERROR_MAP_MINUS, charsize = 1.0
plot_map, TELOG_ERROR_MAP_PLUS, charsize = 1.0

!P.MULTI = [0,1,2]

;look at the statistics of the fits

ratio_map = aia_map_cube
FOR i = 0, nwave-1 DO BEGIN
	newdata = (aia_map_cube[i].data - aia_simul_map_cube[i].data)/(aia_map_cube[i].data)
	ratio_map[i].data = newdata
ENDFOR

;plot_map, ratio_map[0], dmin = -2.0, dmax = 2.0,/cbar
hist = histogram(ratio_map[0].data, binsize = 0.1, min = -2.0, max =2.0, locations = loc)
plot, loc, hist, psym = 10, /nodata, xtitle = 'ratio (obs-fit)/obs', charsize = 1.5

loadct, 0
hsi_linecolors
colors = findgen(nwave)+2
FOR i = 0, nwave-1 DO BEGIN
	hist = histogram(ratio_map[i].data, binsize = 0.1, min = -2.0, max =2.0, locations = loc)
	oplot, loc, hist, psym = 10, linestyle = 0, color = colors[i], thick = 2.0
ENDFOR

legend, wave_, color = colors, linestyle = 0, pspacing = 0.5

hist = histogram(aia_map_cube[0].data, binsize = 1, min = 0, max = 150, locations = loc)

plot, loc, hist, psym = 10, /nodata, xtitle = 'obs flux', charsize = charsize

loadct, 0
hsi_linecolors
colors = findgen(nwave)+1
FOR i = 0, nwave-1 DO BEGIN
	hist = histogram(aia_map_cube[i].data, binsize = 1, min = 0, max = 150, locations = loc)
	oplot, loc, hist, psym = 10, linestyle = 0, color = colors[i]
ENDFOR

legend, wave_, color = colors, linestyle = 0, pspacing = 0.5

hist = histogram(temperature_map[0].data, binsize = 0.1, min = 5.5, max = 7.0 , locations = loc)
plot, loc, hist, psym = 10, /nodata, xtitle = 'log(temperature)', charsize = charsize
oplot, loc, hist, psym = 10, linestyle = 0

hist = histogram(emission_map[0].data, binsize = 0.1, min = min(emission_map[0].data), max = max(emission_map[0].data) , locations = loc)
plot, loc, hist, psym = 10, /nodata, xtitle = 'log(emission measure)', charsize = charsize
oplot, loc, hist, psym = 10, linestyle = 0

hist = histogram(sigma_map[0].data, binsize = 0.1, min = 0, max = 1.0, locations = loc)
plot, loc, hist, psym = 10, /nodata, xtitle = 'sigma map', charsize = charsize
oplot, loc, hist, psym = 10, linestyle = 0

hist = histogram(chisq_map[0].data, binsize = 0.1, min = 0, max = 10 , locations = loc)
plot, loc, hist, psym = 10, /nodata, xtitle = 'chi sq', charsize = charsize
oplot, loc, hist, psym = 10, linestyle = 0
oplot, [1,1], [0,1d10], linestyle = 1

f = file_search(totfile)
!P.multi = 0
IF keyword_set(totfile) AND f[0] NE '' THEN BEGIN
	restore, totfile
	hsi_linecolors
	
	yrange = [min(emlog), max(emlog)]
	plot, telog, emlog, xtitle = 'log(Temperature [MK])', ytitle = 'log(Emmission Measure [cm!U-5!N])', /nodata, yrange = yrange, /ystyle,title = cur_time, charsize = 2.0
	
	oplot, telog, emlog, thick = 2.0, linestyle = 1, psym = symcat(16), color = 6
	oploterr, telog, emlog, emlog_err, /nohat, color = 6, errcolor = 6
ENDIF

IF keyword_set(OUTPLOT) THEN pclose

END