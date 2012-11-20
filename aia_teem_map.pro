pro aia_teem_map,fileset = fileset,wave_,npix,teem_table,teem_map, FILElist = filelist, VERBOSE = verbose, filename_extra = FILENAME_extra, SAVE_DIR = save_dir, DEBUG = debug, XRANGE = xrange, YRANGE = yrange

;+
; Project     : AIA/SDO
;
; Name        : AIA_TEMPMAP 
;
; Category    : Data analysis   
;
; Explanation : calculates EM and Te temperature maps
;		based on single-Gaussian fits in each macropixel
;
; Syntax      : IDL> aia_teem_map,fileset,fov,wave_,npix,teem_table,teem_map
;
; Inputs      : fileset  = common filename part of 6 wavelength FITS images
;		wave_	 = strarr(6) with wavelengths in Angstroem
;		/// fov(4)	 = [i1,j1,i2,j2] pixel ranges of field-of-view - deprecated
;       xrange = fltarr(2) the xrange in arcsec for the subimage
;       yrange = fltarr(2) the yrange in arcsec for the subimage
;		npix	 = size of macropixel (spatial resolution)
;		teem_table = filename of DEM lookup table
;                            (calculated previously with AIA_TEEM_TABLE.PRO)
;		teem_map = savefile containing EM and Te maps 
;
; Outputs     : postscript file <plotname>_col.ps (if io=2)
;
; History     :  3-Mar-2011, Version 1 written by Markus J. Aschwanden
;				 16-May-2011, made fov optional. If not given, uses the
;								whole images. 
;								made fileset optional. Added filelist keyword to just give a list of files.
;		: 5-Oct-2011 - Added search for 1 sigma contour in chi^2 space to calculate errors in telog and tsig - A. Inglis
;		: 1-Feb-2011 - Changed definition of files (see lines 51-55) so that files = filelist directly. Solves issues with AIA2011xxx file formats.
;
; Contact     : aschwanden@lmsal.com
;-

default, save_dir, ''

;_________________________________________________________________________
t1 = systime(0,/seconds)
nwave = n_elements(wave_)
default, filename_extra, ''

files = filelist
texp_ = fltarr(nwave)

FOR iw = 0, nwave-1 DO BEGIN
	read_sdo,files[iw],index,data
    index2map,index,float(data),map
	
	IF keyword_set(xrange) AND keyword_set(yrange) THEN BEGIN
	    sub_map, map, smap, xrange = xrange, yrange = yrange
	    map = smap
	    data = smap.data
    ENDIF
    
	image = map.data
	texp = map.dur
	texp_[iw]=texp 
	dateobs = map.time
	
	IF (iw EQ 0) THEN BEGIN
		dim	= size(image)
		nx	= dim[1]
		ny	= dim[2]
		nxx	= (nx/npix)
		nyy	= (ny/npix)
		i3	= nxx*npix-1
		j3	= nyy*npix-1
		images = fltarr(nxx,nyy,nwave)
		x = get_map_xp(map, /oned)
		y = get_map_yp(map, /oned)
	ENDIF
	
	IF (npix EQ 1) THEN images[*,*,iw] = float(image)/texp

	IF (npix GT 1) THEN BEGIN
		; units of image are in dn/s
		images[*,*,iw] = rebin(float(image[0:i3,0:j3]),nxx,nyy)/texp
		sub_map, map, smap, xrange = [0,i3], yrange = [0,j3], /pixel
		rmap = rebin_map(smap, nxx, nyy)
		x = get_map_xp(rmap, /oned)
		y = get_map_yp(rmap, /oned)
		map = rmap
	ENDIF
	
    ;map.data=(map.data/texp)
    IF iw EQ 0 THEN aia_map_cube = replicate(map, nwave) ELSE aia_map_cube[iw] = map
        
	IF keyword_set(PLOT) THEN BEGIN
		aia_lct, rr, gg, bb, wavelnth=wave_[iw], /load
		plot_map, map, /limb, /log
	ENDIF
ENDFOR

;________________________TEMPERATURE MAP_________________________________ 
restore, save_dir + teem_table, /verbose

dim	= size(flux)
nte	= dim[1]
nsig = dim[2]
nwave = dim[3]
ntot = nte*nsig

; result arrays
telog_map = fltarr(nxx,nyy)
em_map = fltarr(nxx,nyy)
sig_map	= fltarr(nxx,nyy)
flux_dem_map = fltarr(nxx, nyy, nwave)
chi_map = fltarr(nxx,nyy)
em_errmap = fltarr(nxx,nyy)
telog_errmap = fltarr(nxx,nyy,2)
sig_errmap = fltarr(nxx,nyy,2)
telog_errmap_symmetric = fltarr(nxx,nyy)
sig_errmap_symmetric = fltarr(nxx,nyy)
chi_map6 = fltarr(nxx,nyy,6)

chi_map_hsi = fltarr(nxx,nyy)
telog_map_hsi = fltarr(nxx,nyy)
sig_map_hsi = fltarr(nxx,nyy)

telog_best	= 0.
em_best	= 0.
sig_best = 0.

telog_best_hsi = 0.
em_best_hsi = 0.
sig_best_hsi = 0.

; the number of fit parameters
nfree = 3

temp = get_sun(dateobs, sd = solar_radius)

FOR j = 0, nyy-1 DO BEGIN

	; the following code makes sure that regions with radial distance beyond r0
	; are not considered.
	ry = sqrt(x^2 + y[j]^2)
	index = where(ry LE solar_radius, count)
	
	if count NE 0 THEN BEGIN
		i1 = min(index)
		i2 = max(index) 
	ENDIF ELSE continue
	
	FOR i = i1, i2-1 DO BEGIN
	
		flux_obs = reform(images[i,j,*])
		counts = flux_obs*texp_
		;stop
		noise = sqrt(counts)/texp_

		chimin = 9999.
		chi6min = 9999.
		chimin_hsi=9999.
		flux_dem_best = 9999.
		chi2d = fltarr(nte,nsig)
		chi_hsi = fltarr(nte,nsig)
		
		telog_err = [0,0]
		sig_err = [0,0]
		sig_errsymmetric = 0.
		telog_errsymmetric = 0.
		
		FOR k=0, nte-1 DO BEGIN
   			FOR l = 0, nsig-1 DO BEGIN
				flux_dem1 = reform(flux[k,l,*])
				em1	= total(flux_obs)/total(flux_dem1)
				em1_err = sqrt(total(noise^2)) / total(flux_dem1)
				flux_dem = flux_dem1 * em1
				chi	= sqrt(total((flux_obs-flux_dem)^2/noise^2)/(nwave-nfree))
				chi6 = abs(flux_obs-flux_dem)/noise
				chi2d[k,l] = chi
				;add in RHESSI
				em_hsi=em1 /1e25
				em_hsi=em_hsi/1e24    ;
				em_hsi=em_hsi*39.*1e15 
				IF (k eq 0) THEN BEGIN
				em_hsi=em_hsi * (10^telog(1) - 10^telog(0))
				ENDIF ELSE BEGIN
				em_hsi=em_hsi * (10^telog(k) - 10^telog(k-1))
				ENDELSE	
			
				get_hsi_table_entry,[em_hsi,telog(k),tsig(l)],model_count_flux,real_count_flux,axis,summary, obj=obj
				chi_hsi[k,l] = summary.spex_summ_chisq			
				IF (chi le chimin) THEN BEGIN
					chimin = chi
					chi6min	= chi6
					em_best	= alog10(em1)
					em_best_err = em1_err/em_best
					telog_best	= telog[k]
					sig_best = tsig[l]
					flux_dem_best = flux_dem
				ENDIF
				IF (chi_hsi[k,l] le chimin_hsi) THEN BEGIN
					chimin_hsi=chi_hsi[k,l]
					em_best_hsi=alog10(em1)
					telog_best_hsi = telog[k]
					sig_best_hsi = tsig[l]
				ENDIF
			ENDFOR 
		ENDFOR

		; find errors in tsig and telog by finding the 1 sigma contour in chi^2 space for each pixel - A. Inglis 5-Oct-2011
		IF (chimin LT 9999.) THEN BEGIN

			;find the 1 sigma contour in chi^2-space to get errors in tsig and telog
			contour, chi2d, telog, tsig, levels = [chimin +2.3], path_xy=path_xy, path_info=path_info, /path_data_coords, closed=0
			IF exist(path_xy) THEN BEGIN
				telog_err = [telog_best - min(path_xy[0,*]), max(path_xy[0,*]) - telog_best]
				sig_err = [sig_best - min(path_xy[1,*]), max(path_xy[1,*]) - sig_best]

				;error bars are asymmetric - symmetrise!
				telog_errsymmetric = max(abs(telog_err))	
				sig_errsymmetric = max(abs(sig_err))
			ENDIF 
		ENDIF

		IF ((keyword_set(VERBOSE)) AND (chimin LT 9999.)) THEN BEGIN
			print, 'i, j =          ', float(i), float(j)
			print, 'chimin =        ', chimin
			print, 'EM =            ', em_best
			print, 'Log(Temp) = ' + num2str(telog_best) + ' +/- ' + num2str(telog_errsymmetric)
			print, 'Sigma = ' + num2str(sig_best) + ' +/- ' + num2str(sig_errsymmetric)

			print, 'Wavelength:     ', float(wave_)
			print, '----------      ', replicate('-------------', nwave)
			print, 'Observed flux:  ', flux_obs
			print, 'Fitted flux:    ', flux_dem_best
			print, '----------      ', replicate('-------------', nwave)
			print, 'Fitted/Obs flux:', flux_dem_best/flux_obs
		
			window, 0
			loadct, 0
			hsi_linecolors
			nlevels = 20
			levels = chimin * (findgen(nlevels)*0.1 + 1.0)
			anot = strarr(20)
			FOR k = 0, nlevels-1 DO anot[k] = num2str(levels[k])
			contour, chi2d, telog, tsig, levels = levels, c_annotation = anot, xtitle = 'log[Temp]', ytitle = 'sig'
			oplot, [telog_best], [sig_best], psym = symcat(16), color = 6
			leg = num2str(chimin) + '[' + num2str(telog_best) + ',' + num2str(sig_best) + ']'
			legend, leg, psym = 4
			contour,chi2d,telog, tsig, levels=[chimin + 2.3],/over,thick=2, color = 6
			
			window, 1
			emlog =em_best*exp(-(telog-telog_best)^2/(2.*sig_best^2))
			plot, telog, emlog, yrange=minmax(emlog),xtitle='Temperature  log(T)',$
   			ytitle='Emission measure  log(EM [cm!U-5!N K!U-1!N])'
   			;r = chianti_spec_from_dem(telog, emlog, /plot)
   			;save, emlog, telog, r, filename = 'xray_spectrum_i' + num2str(i) + '_j' + num2str(j) + '.sav'
		ENDIF
		
		flux_dem_map[i,j,*] = flux_dem_best
		em_map[i,j] = em_best
		telog_map[i,j] = telog_best
		sig_map[i,j] = sig_best
		chi_map[i,j] = chimin
		chi_map6[i,j,*] = chi6min
		telog_errmap[i,j,0:1] = telog_err[0:1]
		sig_errmap[i,j,0:1] = sig_err[0:1]
		telog_errmap_symmetric[i,j] = telog_errsymmetric
		sig_errmap_symmetric[i,j] = sig_errsymmetric
	
		chi_map_hsi[i,j] = chimin_hsi
		telog_map_hsi[i,j] = telog_best_hsi
		sig_map_hsi[i,j] = sig_best_hsi
	ENDFOR
	
	t2 = systime(0,/seconds)
	elapsed_time_min = (t2-t1)/60.0
	percent_done = float(j)/nyy * 100
	percent_rate = percent_done/elapsed_time_min
	remaining_time = (100 - percent_done)/percent_rate
 	IF (j mod 10) eq 0 THEN print, num2str(percent_done, length = 2) + '% done. ' + num2str(elapsed_time_min) + ' min elapsed, ' + num2str(remaining_time) + ' min left'
ENDFOR

aia_simul_map_cube = aia_map_cube

FOR i = 0, nwave-1 DO BEGIN
	aia_simul_map_cube[i].data = flux_dem_map[*,*,i]
	aia_simul_map_cube[i].id = 'Simulated SDO/AIA ' + num2str(wave_[i]) + ' A'
ENDFOR

temperature_map = map
temperature_map.data = telog_map
temperature_map.id = 'log(Temperature [MK])'
temperature_map.time = dateobs

emission_map = map
emission_map.data = em_map
emission_map.id = 'log(Emission Measure [cm!U-3!N]'
emission_map.time = dateobs

sigma_map = map
sigma_map.data = sig_map
sigma_map.id = 'Sigma'
sigma_map.time = dateobs

chisq_map = map
chisq_map.data = chi_map
chisq_map.id = textoidl('\chi^2')
chisq_map.time = dateobs

sig_error_map_minus = map
sig_error_map_minus.data = sig_errmap[*,*,0]
sig_error_map_minus.id = 'Temp sigma err -'
sig_error_map_minus.time = dateobs

telog_error_map_minus = map
telog_error_map_minus.data = telog_errmap[*,*,0]
telog_error_map_minus.id = 'log(Temp)  err -'
telog_error_map_minus.time = dateobs

sig_error_map_plus = map
sig_error_map_plus.data = sig_errmap[*,*,1]
sig_error_map_plus.id = 'Temp sigma err +'
sig_error_map_plus.time = dateobs

telog_error_map_plus = map
telog_error_map_plus.data = telog_errmap[*,*,1]
telog_error_map_plus.id = 'log(Temp)  err +'
telog_error_map_plus.time = dateobs

;________________________STATISTICS OF CHI-2 FITS_________________
print,'Statistics of chi2='
statistic,chi_map
print,'Statistics of chi2 (94 A)'
statistic,chi_map6(*,*,5)
print,'log(EM)-range = ',minmax(em_map)

;________________________SAVE MAPS________________________________
;save,filename=teem_map,telog_map,em_map,sig_map,chi_map,dateobs 
save,filename=save_dir + teem_map,telog_map,em_map,sig_map,chi_map,sig_errmap_symmetric,telog_errmap_symmetric,temperature_map,emission_map,sigma_map,chisq_map,aia_map_cube, $
dateobs , aia_simul_map_cube, chi_map_hsi, telog_best_hsi, sig_best_hsi
save,temperature_map,emission_map,sigma_map,chisq_map,sig_error_map_minus,sig_error_map_plus,telog_error_map_minus,telog_error_map_plus,aia_simul_map_cube, filename=save_dir + 'tempmap' + filename_extra + '.sav'

print,'TE+EM maps saved in file : ',teem_map
t2	=systime(0,/seconds)

print,'Computation time = ',(t2-t1)/60.0,' min'

IF keyword_set(DEBUG) THEN stop

END
