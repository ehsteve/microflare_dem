PRO aia_teem_hsi_aia_compare, OUTPS = outps

loadct, 0
hsi_linecolors

ospex_fit_results = '/Users/schriste/Desktop/flare_0_rhessi/ospex_results_20110603_0715_27_apr_2012.fits'
fits2map, hsi_image, hsimap

fits2map, '/Users/schriste/Desktop/flare_0_rhessi/hsi_image_20110603_071546.fits', map
index = where(map.data GE 0.5*(max(map.data)), count)
emission_area = (3.6*716d5)^2         ;3.6 arcsec in km

s = spex_read_fit_results(ospex_fit_results)

hsi_timerange = anytim(s.spex_summ_time_interval,/yoh)

deemtot_dir = '/Users/schriste/Desktop/aia_deem_files/'
f = file_search(deemtot_dir + 'teem_tot*_mask_.sav')
times = unbreak_time(strmid(file_basename(f), 12, 15))

;hsi_timerange = ['2011/06/03 06:55:11','2011/06/03 07:01:11']

con1 = anytim(times) GE anytim(hsi_timerange[0])
con2 = anytim(times) LE anytim(hsi_timerange[1])
index = where(con1 AND con2, count)

FOR i = 0, count-1 DO BEGIN

	restore, f[index[i]]
	IF i EQ 0 THEN BEGIN
		telog_total = telog
		emlog_total = emlog
		emlog_err_total = emlog_err
	ENDIF ELSE BEGIN
		emlog_total = alog10(10^emlog_total + 10^emlog)
		;emlog_err_total = sqrt(emlog_err_total^2 + emlog_err^2)
	ENDELSE
	
ENDFOR

ytitle = 'log(Emission Measure [cm!U-5!N])'
xtitle = 'log(Temperature [K])'
plot, telog, emlog_total, xtitle = xtitle, ytitle = ytitle, yrange = minmax(emlog_total)
FOR i = 0, count-1 DO BEGIN
	restore, f[index[i]], /verbose
	oplot, telog, emlog
ENDFOR

aiaflux = chianti_spec_from_dem(telog, emlog_total, emission_area = emission_area)

IF keyword_set(OUTPS) THEN popen, 'fig_proposal_aiarhessi', xsize = 7, ysize = 4

!P.MULTI = [0,2,1]

ytitle = 'log(DEM (cm!E-5!N log (K)!E-1!N)'
xtitle = 'log(Temperature [K])'
plot, telog, emlog_total, yrange = minmax(emlog_total), /nodata, xtitle = xtitle, ytitle = ytitle, title = 'AIA DEM'
oplot, telog, emlog_total, thick = 3

ytitle = 'photons s!U-1!N cm!U-1!N keV!U-1!N'
xtitle = 'Energy [keV]'

plot, aiaflux[0,*], aiaflux[1,*], /ylog, /xlog, /xstyle, /ystyle, xtitle = xtitle, ytitle = ytitle, /nodata, charsize = 1.0, ytickf = 'exp1', title = title, xrange = [3,30]

oplot, aiaflux[0,*], aiaflux[1,*], thick = 2.0, psym = 10, linestyle = 2

conv_fact = s.spex_summ_conv
area = s.spex_summ_area
ewidth = get_edges(s.spex_summ_energy,/width)
ct_flux = s.spex_summ_ct_rate  / area / rebin(ewidth, size(s.spex_summ_ct_rate, /dim))
ph_flux = ct_flux / conv_fact

oplot, s.spex_summ_energy, ph_flux, thick = 2.0, psym = 10

legend, ['rhessi', 'aia'], thick = [2,2], linestyle = [0,2], /right, charsize = 0.8

IF keyword_set(OUTPS) THEN pclose

stop

END