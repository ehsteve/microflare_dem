PRO plot_surfaces,ps=ps

restore,'aia_hsi_fit_results_fineres.sav',/verbose

restore,'hsi_teem_table.sav',/verbose

a=aia_hsi_fit_results

h=hsi_teem_table

loadct,3


IF keyword_set(ps) THEN BEGIN
set_plot,'ps'
device,encaps=1,color=1,bits_per_pixel=16,xsize=30,ysize=25,filename='surfaces1.ps'
ENDIF ELSE BEGIN
window,1,xsize=1024,ysize=1024
ENDELSE
!p.multi=[0,2,2]

shade_surf,a.chi2d_hsi,findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,charsize=3,/zlog,az=50,ax=30,xtitle='log T (K)',ytitle='sigma',ztitle='chi^2';,title='RHESSI chi^2 map'
shade_surf,a.chi2d_hsi,findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,charsize=3,/zlog,az=-130,ax=30,xtitle='log T (K)',ytitle='sigma',ztitle='chi^2';,title='RHESSI chi^2 map'

shade_surf,a.chi_2d,findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,charsize=3,az=50,ax=30,/zlog,xtitle='log T (K)',ytitle='sigma',ztitle='chi^2';,title='AIA chi^2 map'
shade_surf,a.chi_2d,findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,charsize=3,az=-130,ax=30,/zlog,xtitle='log T (K)',ytitle='sigma',ztitle='chi^2';,title='AIA chi^2 map'

xyouts,0.43,0.95,'RHESSI chi^2 map',/norm,charsize=2
xyouts,0.43,0.45,'AIA chi^2 map',/norm,charsize=2

IF keyword_set(ps) THEN BEGIN
device,/close
device,encaps=1,color=1,bits_per_pixel=16,xsize=30,ysize=25,filename='surfaces2.ps'
ENDIF ELSE BEGIN
window,2,xsize=1024,ysize=1024
ENDELSE

shade_surf,h.chi_best_hsi,h.telog,h.tsig,charsize=3,az=50,ax=30,zrange=[20,80],xtitle='log T (K)',ytitle='sigma',ztitle='chi^2';,title='RHESSI chi^2 map'
shade_surf,h.chi_best_hsi,h.telog,h.tsig,charsize=3,az=-130,ax=30,zrange=[20,80],xtitle='log T (K)',ytitle='sigma',ztitle='chi^2';,title='RHESSI chi^2 map'

shade_surf,h.em_best_hsi,h.telog,h.tsig,charsize=3,az=50,ax=30,/zlog,xtitle='log T (K)',ytitle='sigma',ztitle='EM (cm^-3)';,title='RHESSI EM map'
shade_surf,h.em_best_hsi,h.telog,h.tsig,charsize=3,az=-130,ax=30,/zlog,xtitle='log T (K)',ytitle='sigma',ztitle='EM (cm^-3)';,title='RHESSI EM map'

xyouts,0.43,0.95,'RHESSI chi^2 map',/norm,charsize=2
xyouts,0.43,0.45,'RHESSI EM map',/norm,charsize=2

IF keyword_set(ps) THEN BEGIN
device,/close
device,encaps=1,color=1,bits_per_pixel=16,xsize=25,ysize=15,filename='surfaces3.ps'
ENDIF ELSE BEGIN
window,3,xsize=1024,ysize=512
ENDELSE


!p.multi=[0,2,1]

loadct,39
nlevels = 50
levels = a.chimin_hsi * (findgen(nlevels)*0.2 + 1.0)
contour, a.chi2d_hsi,findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,levels=levels,xtitle='log T (K)',ytitle='sigma',title='RHESSI chi^2 map'
oplot, [a.telog_best_hsi], [a.sig_best_hsi], psym = 2, color = 240,symsize=2
contour,a.chi2d_hsi,findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01, levels=[a.chimin_hsi + 2.3],/over,thick=2, color = 240

levels = a.chimin * (findgen(nlevels)*0.2 + 1.0)
contour, a.chi_2d,findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,levels=levels,xtitle='log T (K)',ytitle='sigma',title='AIA chi^2 map'
oplot, [a.telog_best], [a.sig_best], psym = 2, color = 240,symsize=2
contour,a.chi_2d,findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01, levels=[a.chimin + 2.3],/over,thick=2, color = 240


IF keyword_set(ps) THEN BEGIN
device,/close
device,encaps=1,color=1,bits_per_pixel=16,xsize=20,ysize=15,filename='surfaces4.ps'
ENDIF ELSE BEGIN
window,4
ENDELSE


nlevels = 50
levels = MIN(h.chi_best_hsi,pos) * (findgen(nlevels)*0.2 + 1.0)
contour, h.chi_best_hsi,h.telog,h.tsig,levels=levels,xtitle='log T (K)',ytitle='sigma',title='RHESSI chi^2 map'
pos=array_indices(h.chi_best_hsi,pos)
oplot, [h.telog[pos[0]]], [h.tsig[pos[1]]], psym = 2, color = 240,symsize=2
contour, h.chi_best_hsi,h.telog,h.tsig,levels=[MIN(h.chi_best_hsi) + 2.3],/over,thick=2, color = 240

nlevels = 50
levels = MIN(h.em_best_hsi,pos) * (findgen(nlevels)*0.2 + 1.0)
contour, h.em_best_hsi,h.telog,h.tsig,levels=levels,xtitle='log T (K)',ytitle='sigma',title='RHESSI EM map'
pos=array_indices(h.em_best_hsi,pos)
;oplot, [h.telog[pos[0]]], [h.tsig[pos[1]]], psym = 2, color = 240,symsize=2
;contour, h.em_best_hsi,h.telog,h.tsig,levels=[a.chimin_hsi + 2.3],/over,thick=2, color = 240

IF keyword_set(ps) THEN BEGIN
device,/close
device,encaps=1,color=1,bits_per_pixel=16,xsize=20,ysize=15,filename='surfaces5.ps'
ENDIF ELSE BEGIN
window,5
ENDELSE

loadct,3
!p.multi=[0,1,1]
shade_surf,a.em_2d_hsi/h.em_best_hsi,findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,charsize=3,az=50,ax=30,/zlog,zrange=[1e-2,1e4]

xyouts,0.43,0.95,'AIA / RHESSI EM ratio',/norm,charsize=2


IF keyword_set(ps) THEN BEGIN
device,/close
device,encaps=1,color=1,bits_per_pixel=25,xsize=30,ysize=15,filename='surfaces6.ps'
ENDIF ELSE BEGIN
window,5
ENDELSE
restore,'aia_hsi_fit_results_ratio_test.sav',/verbose
r=aia_hsi_fit_results


!p.multi=[0,3,2]
shade_surf,r.ratio_test[*,*,0],findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,charsize=3,az=50,ax=40,/zlog,xtitle='log T (K)',ytitle='sigma',ztitle='FLUX_OBS/FLUX_MOD'
shade_surf,r.ratio_test[*,*,1],findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,charsize=3,az=50,ax=40,/zlog,xtitle='log T (K)',ytitle='sigma',ztitle='FLUX_OBS/FLUX_MOD'
shade_surf,r.ratio_test[*,*,2],findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,charsize=3,az=50,ax=40,/zlog,xtitle='log T (K)',ytitle='sigma',ztitle='FLUX_OBS/FLUX_MOD'
shade_surf,r.ratio_test[*,*,3],findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,charsize=3,az=50,ax=40,/zlog,xtitle='log T (K)',ytitle='sigma',ztitle='FLUX_OBS/FLUX_MOD'
shade_surf,r.ratio_test[*,*,4],findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,charsize=3,az=50,ax=40,/zlog,xtitle='log T (K)',ytitle='sigma',ztitle='FLUX_OBS/FLUX_MOD'
shade_surf,r.ratio_test[*,*,5],findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,charsize=3,az=50,ax=40,/zlog,xtitle='log T (K)',ytitle='sigma',ztitle='FLUX_OBS/FLUX_MOD'

xyouts,0.2,0.95,'131A',/norm,charsize=1.5
xyouts,0.5,0.95,'171A',/norm,charsize=1.5
xyouts,0.8,0.95,'193A',/norm,charsize=1.5
xyouts,0.2,0.45,'211A',/norm,charsize=1.5
xyouts,0.5,0.45,'335A',/norm,charsize=1.5
xyouts,0.8,0.45,'94A',/norm,charsize=1.5

IF keyword_set(ps) THEN BEGIN
device,/close
device,encaps=1,color=1,bits_per_pixel=25,xsize=30,ysize=15,filename='surfaces7.ps'
ENDIF ELSE BEGIN
window,5
ENDELSE
;restore,'aia_hsi_fit_results_ratio_test.sav',/verbose
;r=aia_hsi_fit_results


contour,alog10(r.ratio_test[*,*,0]),findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,xtitle='log T (K)',ytitle='sigma',title='FLUX_OBS/FLUX_MOD',nlevels=100
ssw_legend,'131A',/right
contour,alog10(r.ratio_test[*,*,1]),findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,xtitle='log T (K)',ytitle='sigma',title='FLUX_OBS/FLUX_MOD',nlevels=100
ssw_legend,'171A',/right
contour,alog10(r.ratio_test[*,*,2]),findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,xtitle='log T (K)',ytitle='sigma',title='FLUX_OBS/FLUX_MOD',nlevels=100
ssw_legend,'193A',/right
contour,alog10(r.ratio_test[*,*,3]),findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,xtitle='log T (K)',ytitle='sigma',title='FLUX_OBS/FLUX_MOD',nlevels=100
ssw_legend,'211A',/right
contour,alog10(r.ratio_test[*,*,4]),findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,xtitle='log T (K)',ytitle='sigma',title='FLUX_OBS/FLUX_MOD',nlevels=100
ssw_legend,'335A',/right
contour,alog10(r.ratio_test[*,*,5]),findgen(30)*0.05 + 6.,findgen(78)*0.005 + 0.01,xtitle='log T (K)',ytitle='sigma',title='FLUX_OBS/FLUX_MOD',nlevels=100
ssw_legend,'94A',/right
IF keyword_set(ps) THEN BEGIN
device,/close
set_plot,'x'
ENDIF

END