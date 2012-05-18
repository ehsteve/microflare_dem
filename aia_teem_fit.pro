FUNCTION aia_teem_fit, flux_obs, flux_obs_error, FLUX_fit = flux_fit, fit_params = fit_params, PLOT = plot

restore, '/Users/schriste/Dropbox/idl/aia_deem/teem_table.sav', /verbose

dim	= size(flux)
nte	= dim[1]
nsig = dim[2]
nwave = dim[3]
ntot = nte*nsig
nfree = 3

chi2d = fltarr(nte,nsig)
chimin = 9999.
chi6min = 9999.
flux_dem_best = 9999.

FOR i = 0, nte-1 DO BEGIN
	FOR j = 0, nsig-1 DO BEGIN
	
		flux_dem1 = reform(flux[i,j,*])
		
		em1	= total(flux_obs)/total(flux_dem1)
		em1_err = sqrt(total(flux_obs_error^2)) / total(flux_dem1)
		flux_dem = flux_dem1 * em1
		
		chi	= sqrt(total((flux_obs-flux_dem)^2/flux_obs_error^2)/(nwave-nfree))
		chi6 = abs(flux_obs-flux_dem)/flux_obs_error
		chi2d[i,j] = chi
			
		IF (chi LE chimin) THEN BEGIN
			chimin = chi
			chi6min	= chi6
			em_best	= alog10(em1)
			em_best_err = em1_err/em_best
			telog_best	= telog[i]
			sig_best = tsig[j]
			flux_dem_best = flux_dem
		ENDIF
		
		print, chi
	ENDFOR 
ENDFOR

flux_fit = flux_dem_best
residuals = flux_obs - flux_fit

;calculate the error on the parameters
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

IF keyword_set(PLOT) THEN BEGIN
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
	
	window, 2
	!p.multi= 0

	plot,indgen(nwave),flux_obs,/ylog,psym=6,$
	xrange=[-1,nwave],xtickf='(a1)',xticks=nwave+1,ystyle=16,thick=3,$
	ytit='DN s!U-1!N',/nodata,$
	yrange=[0.9*min(flux_obs),1.1*max(flux_obs)], xtickn=[' ', wave_, ' '], xtitle = 'Filter'
	oplot,indgen(nwave),flux_obs,psym=6,color=5,thick=1
	oplot,indgen(nwave),flux_fit,psym=7,color=2,thick=2
	for i=0, nwave-1 do oplot, [i,i], flux_obs[i] + [-flux_obs_error[i], flux_obs_error[i]], thick=5,color=5
	legend, ['flux_obs', 'flux_fit'], psym = [6,7], color = [5,2]
	legend, 'chisq='+num2str(chimin), /right
	print, wave_
	;maxr=1.1*max(abs(residuals))
	;plot,indgen(nwave),residuals,xrange=[-1,nwave],xtickn=[' ',wave_,' '],$
	;xticks=nwave+1,ystyle=17,thick=1,yrange=maxr*[-1,1],psym=6,$
	;ytit='Residuals',xtit='Filter'
	;oplot,[-2,nwave],[0,0],lines=1
	;legend, 'chisq='+num2str(chimin)
ENDIF

fit_params = [em_best, telog_best, sig_best]
dem =em_best*exp(-(telog-telog_best)^2/(2.*sig_best^2))

RETURN, [transpose(telog), transpose(dem)]

END
