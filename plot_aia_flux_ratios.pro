PRO plot_aia_flux_ratios, epstein=epstein, OUTPS = outps

;plots modelled vs observed AIA flux for a given model DEM.


;restore,'aia_fit_results_manual.sav',/verbose
;a_man=aia_fit_results

;set_plot,'ps'
;device,encaps=1,filename='aia_flux_ratios_for_hsi_best.ps'

;plot,findgen(6)+1,a_man.flux_dem_3d[15,18,*]/a_man.flux_obs,psym=4,symsize=2,charsize=1.2,xrange=[0,7],xticks=7, $
;xtickname=[' ','131A','171A','193A','211A','335A','94A',' '],ytitle='flux_model/flux_obs',yran=[-0.2,2.5]

;oploterr,findgen(6)+1,a_man.flux_dem_3d[15,18,*]/a_man.flux_obs,findgen(6)*0 + 0.2,psym=4,symsize=2

;oplot,findgen(10),findgen(10)*0 + 1,linestyle=2,thick=2
;leg = num2str(a_man.chimin) + '[' + num2str(22.04) +',' + num2str(6.75) + ',' + num2str(0.1) + ']'

;leg = num2str(a_man.chimin) + '[' + num2str(6.75) + ',' + num2str(0.1) + ']'
;legend, leg, psym = 4,charsize=1.2

;device,/close
;set_plot,'x'



IF keyword_set(epstein) THEN BEGIN
restore,'aia_fit_results_epstein.sav'
a=aia_fit_results
fname='aia_flux_ratios_for_aia_best_epstein.ps'
ENDIF ELSE BEGIN
restore,'aia_fit_results.sav'
a=aia_fit_results
fname='aia_flux_ratios_for_aia_best.ps'
ENDELSE

tt=value_locate(a.telog,a.telog_best)
ss=value_locate(a.tsig,a.sig_best)

IF keyword_set(OUTPS) THEN BEGIN
	set_plot,'ps'
	device,encaps=1,filename=fname
ENDIF

!P.MULTI = [0,1,2]
hsi_linecolors
plot,findgen(6)+1,a.flux_dem_best,charsize=1.2,xrange = [0,7],xticks=7, $
xtickname=[' ','131A','171A','193A','211A','335A','94A',' '], ytitle='AIA DN Numbers', /nodata

oploterr,findgen(6)+1,a.flux_obs,a.flux_obs*0.2, psym = symcat(14),symsize=2, /nohat
oplot, findgen(6)+1, a.flux_dem_best, psym = 7, symsize = 1

ssw_legend, ['AIA observed', 'FIT'], psym = [symcat(14), symcat(14)], charsize = 1.2

plot,findgen(6)+1,a.flux_dem_best/a.flux_obs,psym=4,symsize=2,charsize=1.2,xrange=[0,7],xticks=7, $
xtickname=[' ','131A','171A','193A','211A','335A','94A',' '],ytitle='flux_model/flux_obs',yran=[-0.2,2.5]

oploterr,findgen(6)+1,a.flux_dem_3d[tt,ss,*]/a.flux_obs,findgen(6)*0 + 0.2,psym=symcat(14),symsize=2
oplot,findgen(10),findgen(10)*0 + 1,linestyle=2,thick=2

leg = textoidl('\chi^2 = ') + num2str(a.chimin, length = 5) + ' [log(EM) = ' + num2str(a.em_best, length = 5) + ' cm!U-5!N, log(Tmax) = ' + num2str(a.telog_best, length = 5) + ' (' + num2str(10d^a.telog_best/1d6, length=5) + ' MK), sigma = ' + num2str(a.sig_best) + ']'
ssw_legend, leg, psym = symcat(14),charsize=1.2
!P.MULTI = 0
IF keyword_set(OUTPS) THEN BEGIN
	device,/close
	set_plot,'x'
ENDIF

print,' '
print,'------------------------------'
print,'Wrote file: ',fname
print,'------------------------------'
print,' '


f = file_search('aia_hsi_fit_results_epstein')

IF f[0] NE '' THEN BEGIN
	IF keyword_set(epstein) THEN BEGIN
	restore,'aia_hsi_fit_results_epstein.sav'
	ah=aia_hsi_fit_results
	fname='aia_flux_ratios_for_best_combo_epstein.ps'
	ENDIF ELSE BEGIN
	restore,'aia_hsi_fit_results.sav'
	ah=aia_hsi_fit_results
	fname='aia_flux_ratios_for_best_combo.ps'
	ENDELSE
	
	chi_combo=ah.chi_2d+ah.chi2d_hsi
	m= min(chi_combo,pos)
	pos=array_indices(chi_combo,pos)
	
	IF keyword_set(OUTPS) THEN BEGIN
		set_plot,'ps'
		device,encaps=1,filename=fname
	ENDIF
	
	plot,findgen(6)+1,a.flux_dem_3d[pos[0],pos[1],*]/a.flux_obs,psym=4,symsize=2,charsize=1.2,xrange=[0,7],xticks=7, $
	xtickname=[' ','131A','171A','193A','211A','335A','94A',' '],ytitle='flux_model/flux_obs',yran=[-0.2,4.0]
	
	oploterr,findgen(6)+1,a.flux_dem_3d[pos[0],pos[1],*]/a.flux_obs,findgen(6)*0 + 0.2,psym=4,symsize=2
	oplot,findgen(10),findgen(10)*0 + 1,linestyle=2,thick=2
	
	leg = num2str(min(chi_combo)) + '[' + num2str(alog10(ah.em_2d[pos[0],pos[1]])) +',' + num2str(ah.telog[pos[0]]) + ',' + num2str(ah.tsig[pos[1]]) + ']'
	legend, leg, psym = 4,charsize=1.2
	
	IF keyword_set(OUTPS) THEN BEGIN
		device,/close
		set_plot,'x'
	ENDIF

print,' '
print,'------------------------------'
print,'Wrote file: ',fname
print,'------------------------------'
print,' '

ENDIF


END