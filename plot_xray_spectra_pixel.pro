PRO plot_xray_spectra_pixel, files = files, OUTPS = outps

IF keyword_set(OUTPS) THEN BEGIN
    loadct, 0
    hsi_linecolors
    popen, 'plot_xray_spectra_pixel', xsize = 7, ysize = 6
ENDIF

default, dir, '/Users/schriste/Dropbox/idl/aia_deem/'

restore, dir + 'teem_data_00020110215_000009_innercore_npix10.sav', /verbose

IF NOT keyword_set(files) THEN BEGIN
	s = size(te_map)
	dimx = s[1]
	dimy = s[2]
ENDIF ELSE BEGIN
	dimx = 1
	dimy = n_elements(files)
ENDELSE 

fexxi_fexvii_data = te_map
fexxiii_data = te_map
ni_vi_data = te_map

FOR i = 0, dimx-1 DO BEGIN
	FOR j = 0, dimy-1 DO BEGIN
	
		IF NOT keyword_set(files) THEN fname = dir + 'xray_spectrum_i' + num2str(i) + '_j' + num2str(j) + '.sav' $
			ELSE fname = files[j]
			
		print, 'Search for ', fname
		f = file_search(fname)
		IF f[0] NE '' THEN BEGIN
			
			restore, f[0], /verbose
			
			!P.MULTI = [0,1,4]
			
			plot, telog, emlog, xtitle = 'log(Temperature [K]', ytitle = 'log(Emission measure)', /nodata, charsize = 2.0, ytickf = 'exp1', title = 'i = ' + num2str(i) + ', j = ' + num2str(j)
			oplot, telog, emlog, thick = 1.0
			
			plot, r[0,*], r[1,*], /nodata, xtitle = 'Energy [keV]', ytitle = 'Flux', /ylog, charsize = 2.0
			oplot, r[0,*], r[1,*]
			
			xrange = [0.95, 1.15]
			con1 = r[0,*] GE xrange[0]
			con2 = r[0,*] LE xrange[1]
			index = where(con1 AND con2)
			
			plot, r[0,index], r[1,index], /nodata, xtitle = 'Energy [keV]', ytitle = 'Flux', /ylog, charsize = 2.0
			oplot, r[0,index], r[1,index], psym = 10
			
			xrange = [6.0, 7.0]
			con1 = r[0,*] GE xrange[0]
			con2 = r[0,*] LE xrange[1]
			index = where(con1 AND con2)
			
			plot, r[0,index], r[1,index], /nodata, xtitle = 'Energy [keV]', ytitle = 'Flux', /ylog, charsize = 2.0
			oplot, r[0,index], r[1,index], psym = 10
			
			fexxi_fexvii = [1.020, 1.024]
			fexxiii = [1.008, 1.014]
			ni_vi = [0.41, 0.44]
			
			xrange = fexxi_fexvii
			con1 = r[0,*] GE xrange[0]
			con2 = r[0,*] LE xrange[1]
			index = where(con1 AND con2)
			fexxi_fexvii_data[i,j] = total(r[1,index])
			
			xrange = fexxiii
			con1 = r[0,*] GE xrange[0]
			con2 = r[0,*] LE xrange[1]
			index = where(con1 AND con2)
			fexxiii_data[i,j] = total(r[1,index])
			
			xrange = ni_vi
			con1 = r[0,*] GE xrange[0]
			con2 = r[0,*] LE xrange[1]
			index = where(con1 AND con2)
			ni_vi_data[i,j] = total(r[1,index])
			
				
		ENDIF
	ENDFOR
ENDFOR
loadct, 0


!P.multi = 0

plot_map, emission_map, dmin = 22
plot_map, temperature_map, dmin = 6

map_fexxi_fexvii = aia_map_cube[0]
map_fexxi_fexvii.data = fexxi_fexvii_data
map_fexxi_fexvii.id = 'Fe XXI, FeXVII'

plot_map, map_fexxi_fexvii

map_fexxiii = aia_map_cube[0]
map_fexxiii.data = fexxiii_data
map_fexxiii.id = 'Fe XXIII'

plot_map, map_fexxiii

ni_vi = aia_map_cube[0]
ni_vi.data = ni_vi_data
ni_vi.id = 'Ni VI'

plot_map, ni_vi

stop

IF keyword_set(OUTPS) THEN pclose

END