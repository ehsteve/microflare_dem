PRO aia_hsi_dem_results, filename

;FUNCTION to plot results from aia_hsi_dem_analysis

restore, filename
res = aia_fit_results

print, 'Analysis Summary'
print, '================'
print, ''
print, 'Model = ' + res.model
print, 'Fit Results'
print, '-----------'
print, 'log(Te) = ' + num2str(res.telog_best, length = 5) + ' (' + num2str(10^res.telog_best/1d6, length = 5) + ' MK)'
print, 'sigma = ' + num2str(res.sig_best, length = 5)
print, 'Tmin = ' + num2str(res.telog_best - res.sig_best, length = 5)  + ' (' + num2str(10^(res.telog_best - res.sig_best)/1d6, length = 5) + ' MK)'
print, 'Tmax = ' + num2str(res.telog_best + res.sig_best, length = 5)  + ' (' + num2str(10^(res.telog_best + res.sig_best)/1d6, length = 5) + ' MK)'
print, 'Emission measure = ' + num2str(10d^res.em_best) + ' cm^5'
print, 'chi^2 = ' + num2str(res.chimin, length = 5)

window, 0
plot_aia_flux_ratios, /epstein

window, 1
plot_dem_result, /epstein

plot_xray_spectrum_from_dem, res.em_best, res.telog_best, res.sig_best, n=10, flare_area=res.flare_area, /epstein, eph = eph, photons = ph

oplot, e.mean, ph_flux

ssw_legend, leg

foxsi_area = get_foxsi_effarea(energy_arr = e.mean, /per_module, /plot)

cts = ph_flux*foxsi_area.eff_area_cm2

plot, e.mean, ph_flux, xrange = [5,20], psym = 10, ytitle = 'FOXSI Count rate', yrange = [0,50]

END

