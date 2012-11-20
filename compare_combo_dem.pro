PRO compare_combo_dem,aia_dem,hsi_dem,combo_dem,telog

restore,'aia_hsi_fit_results_fineres.sav',/verbose

telog=findgen(30)*0.05 + 6.
tsig=findgen(78)*0.005 + 0.01

a=aia_hsi_fit_results



chi_combo=a.chi_2d+a.chi2d_hsi

combo_min=MIN(chi_combo,pos)
pos=array_indices(chi_combo,pos)

telog_best_combo=telog(pos[0])
sig_best_combo=tsig(pos[1])
em_best_combo=a.em_2d[pos]


;plot the dems
aia_dem=a.em_best*exp(-(telog-a.telog_best)^2/((2.*a.sig_best)^2))

hsi_dem=a.em_best_hsi*exp(-(telog-a.telog_best_hsi)^2/((2.*a.sig_best_hsi)^2))

combo_dem=a.em_best*exp(-(telog-telog_best_combo)^2/((2.*sig_best_combo)^2))


;window,1,xsize=756
set_plot,'ps'
device,encaps=1,filename='combo_dem.ps'
plot,telog,aia_dem,thick=2,linestyle=0,yrange=[10,30],ytitle='log EM (cm^-5 K-1)',xtitle = 'log T (K)'

oplot,telog,hsi_dem,thick=2,linestyle=2

oplot,telog,combo_dem,thick=2,linestyle=3

ssw_legend,['AIA','HSI','AIA and HSI'],linestyle=[0,2,3],thick=[2,2,2],/right,charsize=1.2
device,/close
;set_plot,'x'

device,color=1,bits_per_pixel=16,filename='combo_chi_surface.ps'
loadct,3
shade_surf,chi_combo,telog,tsig,charsize=2,/zlog,ax=40,az=70,yrange=[0,0.5],xtitle='log T (K)',ytitle='sigma',ztitle='chi_combo'

device,/close
;set_plot,'x'
device,color=1,bits_per_pixel=16,filename='combo_chi_contour.ps'
;window,1
loadct,0

nlevels = 200
levels = a.chimin_hsi * (findgen(nlevels) + 1.0)
contour,chi_combo,telog,tsig,charsize=2,levels=levels,thick=1,yrange=[0,0.5]

device,/close
set_plot,'x'
END
