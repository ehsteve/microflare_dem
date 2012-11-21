PRO plot_dem_result, epstein = epstein

IF keyword_set(epstein) THEN BEGIN
	restore,'aia_fit_results_epstein.sav'
	a=aia_fit_results
	fname='aia_flux_ratios_for_aia_best_epstein.ps'
ENDIF ELSE BEGIN
	restore,'aia_fit_results.sav'
	a=aia_fit_results
	fname='aia_flux_ratios_for_aia_best.ps'
ENDELSE

IF a.model EQ 'Epstein' THEN epstein = 1 ELSE epstein = 0

temp_range = a.telog_best + a.sig_best*[-1,1]*1.2
plot_trange = [5.5, 7.5]
temp_dx = 0.01

t = temp_dx*findgen((temp_range[1] - temp_range[0])/temp_dx) + temp_range[0]

;params[0] - emission measure value
;params[1] - temperature value
;params[2] - width value

params = [a.em_best, a.telog_best, a.sig_best]
em = get_dem_from_params( t, params,n = 10, epstein=epstein)

plot, t, em, xrange = plot_trange, /nodata, xtitle = 'log(Temperature [K])', ytitle = 'log(Emission Measure [cm!U-5!N K!U-1!N])', charsize = 1.2
oplot, t, em, thick = 2

END