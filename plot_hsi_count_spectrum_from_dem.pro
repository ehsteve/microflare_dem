PRO plot_hsi_count_spectrum_from_dem,epstein=epstein,ps=ps,type=type

default,type,'combo'

IF keyword_set(epstein) THEN BEGIN
   restore,'aia_hsi_fit_results_epstein.sav',/verbose
ENDIF ELSE BEGIN
   restore,'aia_hsi_fit_results.sav',/verbose
ENDELSE

a=aia_hsi_fit_results

tloc=value_locate(a.telog,a.telog_best)
sigloc=value_locate(a.tsig,a.sig_best)

aia_spec=a.model_count_flux_hsi[tloc,sigloc,*]


chi_combo=a.chi_2d + a.chi2d_hsi
m=min(chi_combo,pos)
pos=array_indices(chi_combo,pos)

combo_spec=a.model_count_flux_hsi[pos[0],pos[1],*]

;xticks=[' ',' ','5',' ',' ',' ',' ','10']

IF (type eq 'hsi') THEN BEGIN
   model=a.model_count_flux_hsi_best
   yrange=[0.01,max(model) + 1]
ENDIF ELSE IF (type eq 'aia') THEN BEGIN
   model=aia_spec
   yrange=[0.01,max(model) + 1]
ENDIF ELSE IF (type eq 'combo') THEN BEGIN
   model=combo_spec
   yrange=[0.01,max(model) + 1]
ENDIF ELSE BEGIN
   print,' '
   print,'-----------------'
   print,'Unknown value of TYPE keyword. Aborting.'
   return
ENDELSE

IF keyword_set(ps) THEN BEGIN
   set_plot,'ps'
   IF keyword_set(epstein) THEN BEGIN
      device,encaps=1,filename='rhessi_count_spectrum_vs_model_'+type+'_epstein.ps'
   ENDIF ELSE BEGIN
      device,encaps=1,filename='rhessi_count_spectrum_vs_model_'+type+'.ps'
   ENDELSE
ENDIF

plot,a.axis[0,*],a.real_count_flux_hsi[1,1,*],/xlog,/ylog,thick=3,xthick=3,ythick=3,xrange=[4,15],yrange=yrange,linestyle=3,charsize=1.2, $
xtitle='energy (keV)',ytitle='counts s!U-1!N kev!U-1!N',xstyle=1,ystyle=1,charthick=3;,xticks=5;xtickname=xticks     

oplot,a.axis[0,*],model,thick=3

;show the fitting range on the plot
oplot,(findgen(10000)*0. + 5),(findgen(10000)*0.01),linestyle=1,thick=2
oplot,(findgen(10000)*0. + 12),(findgen(10000)*0.01),linestyle=1,thick=2

ssw_legend,['RHESSI count flux','Model count flux'],linestyle=[3,0],thick=[3,3],/top,/right,charsize=1.2,charthick=3

IF keyword_set(ps) THEN BEGIN
   device,/close
   set_plot,'x'
ENDIF

END
