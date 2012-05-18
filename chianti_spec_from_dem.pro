FUNCTION chianti_spec_from_dem, telog, emlog, ENERGY_RANGE = energy_range, BIN_SIZE_KEV = bin_size_keV, preset_dem = preset_dem, PLOT = plot, EMISSION_AREA = emission_area

;Returns the spectrum observed at Earth. The DEM must be given in log with emission measure in units of cm^-5 and temperature in units of million Kelvin and logged. The DEM can be either an array of only two elements for an isothermal.
;
;KEYWORDS
;		preset_dem - load a preset dem. Choices are
;			flare
;			active_region
;			coronal_hole
;			quiet_sun
;			prominence
;		flare_area - only used for isothermal dem
;		lines - get a struct for the lines in the spectrum
;
;TODO 
;	Does not match with fvth, not sure why. close though.
;
;WRITTEN: Steven Christe (22-Sept-2011)

default, au, 1.49598d13	; in units of cm
default, emission_area, (712d5)^2	; 1 arcsec in cm^2
default, energy_range, [1,30]
default, dem, [alog10(1d49/emission_area), alog10(10d6)]
default, density, 1d9 ; in units of cm^-3
default, bin_size_kev, 0.3

hsi_linecolors

; if I use 10 times less flare_area than I think I would need
; then the fluxes match those gotten by fv_th[1, 23.0/11.6]
; not sure why that is.

factor = emission_area/au^2

;generate transitions
min_angstrom = 12.4/max(energy_range)
max_angstrom = 12.4/min(energy_range)

;for a isothermal emission measure
;logem is in units of cm^-5
;em = 1d49
;logem = alog10(em/flare_area)
;logt = alog10(10d6)

if keyword_set(preset_dem) THEN BEGIN
	filename = !xuvtop+'/dem/' + preset_dem + '.dem'
	read_dem, filename, logt, logem, ref
	dim = n_elements(logt)
	
	ch_synthetic, min_angstrom, max_angstrom, output=transitions, density = density, $
  		ioneq_name=concat_dir(concat_dir(!xuvtop,'ioneq'),'mazzotta_etal.ioneq'), $
  		/photons, /noprot, dem_name=filename, sngl_ion = 'mg_9'
ENDIF ELSE BEGIN
	dim = n_elements(telog)
	IF dim GT 1 THEN BEGIN
	
		;create a temporary file with data in it for chianti to read
		openw, lun, 'datafile.dem', /get_lun
		for i = 0, n_elements(telog)-1 DO printf, lun, telog[i], emlog[i]
		printf, lun, -1
		printf, lun, '%filename: datafile.dem'
		printf, lun, '%dem: blah'
		printf, lun, -1
		close,lun
		free_lun, lun
		logt = telog
		logem = emlog
		ch_synthetic, min_angstrom,max_angstrom, output=transitions, density = density, $
			ioneq_name=concat_dir(concat_dir(!xuvtop,'ioneq'),'mazzotta_etal.ioneq'), $
			/photons, /noprot, dem_name = 'datafile.dem', sngl_ion = 'ni_20'
			
		spawn, 'rm datafile.dem'
	ENDIF ELSE BEGIN
		ch_synthetic, min_angstrom,max_angstrom, output=transitions, density = density, $
			ioneq_name=concat_dir(concat_dir(!xuvtop,'ioneq'),'mazzotta_etal.ioneq'), $
			/photons, /noprot, /all, logem_isothermal = emlog, logt_isothermal = telog
	ENDELSE
ENDELSE

IF keyword_set(PLOT) AND dim NE 1 THEN BEGIN
		!P.multi = [0,1,2]
		if keyword_set(preset_dem) THEN title = preset_dem ELSE title = ''
		plot, logt, logem, /nodata, yrange = [min(logem), max(logem)], /ystyle, ytitle = 'log(Emission Measure [cm!U-5!N])', xtitle = 'log(Temperature [K])', charsize = 2.0, title = title
		oplot, logt, logem, psym = -5
ENDIF

;bin_size_keV = 0.001
bin_size_ang = 0.001
;instr_fwhm_kev = 0.0025
make_chianti_spec, transitions, lambda, struct, bin_size = bin_size_ang, $
	wrange = [min_angstrom, max_angstrom], $
	abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund'), $
 	/photons, /continuum

flux = struct.spectrum
lambda = struct.lambda
 
;convert units of flux from photons cm-2 sr-1 s-1 Angstroms-1
;to photons cm-2 sr-1 s-1 keV-1
energy_arr = 12.4/lambda

s = sort(energy_arr)
;this energy arr has unequally spaced bins
energy_arr = energy_arr[s]
flux = flux * lambda^2/12.4

;now also sort flux
flux = flux[s]

mine = energy_range[0]
maxe = energy_range[1]
nenergy_arr = (findgen((maxe - mine)/bin_size_kev)*bin_size_kev + mine)

factor = factor

;interpolate to equally spaced bins with bins 1 eV apart
nflux = interpol(flux, energy_arr, nenergy_arr)

photon_flux = nflux * factor

IF keyword_set(PLOT) THEN BEGIN
	ytitle = 'photons cm!U-2!N s!U-1!N keV!U-1!N'
	xtitle = 'Energy [keV]'
	yrange = minmax(photon_flux)
	xrange = minmax(nenergy_arr)
	IF keyword_set(preset_dem) THEN title = preset_dem ELSE title = ''
	
	plot, nenergy_arr, photon_flux, /ylog, /xlog, yrange = yrange, xrange = xrange, /xstyle, /ystyle, xtitle = xtitle, ytitle = ytitle, /nodata, charsize = 2.0, ytickf = 'exp1', title = title

	oplot, nenergy_arr, photon_flux, thick = 1.0, psym = 10
ENDIF

RETURN, [transpose(nenergy_arr), transpose(photon_flux)]

END