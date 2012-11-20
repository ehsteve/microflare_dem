; NAME:
;		PLOT_AIA_RESPONSE
;
; PURPOSE:
; Plots the AIA response curves for each EUV wavelength as a function of temperature
;
; CATEGORY:
;	AIA, EUV
;
; CALLING SEQUENCE:
;	plot_aia_response
;
;INPUTS:
;	NONE
;
;OUTPUTS:
;	aia_response_curves.ps
;
;WRITTEN: Andrew Inglis, 2012/07/30



PRO plot_aia_response

loadct,0
resp=aia_get_response(/temp)

set_plot,'ps'
device,encaps=1,color=1,bits_per_pixel=8,filename='aia_response_curves.ps'

plot,resp.logte,resp.all(*,0),/ylog,yrange=[1e-28,1e-23],xrange=[5,8],thick=3,charsize=1.2,xtitle='log T (K)',ytitle='phot cm!U5!N s!U-1!N pix!U-1!N',/nodata

aia_lct,wave='94',/load
oplot,resp.logte,resp.all(*,0),thick=3,color=150
al_legend,'94A',linestyle=[0],color=[150],thick=3,position=[7.3,4.5e-24],linsize=0.5,box=0
aia_lct,wave='131',/load
oplot,resp.logte,resp.all(*,1),thick=3,color=150
al_legend,'131A',linestyle=[0],color=[150],thick=3,position=[7.3,2.9e-24],linsize=0.5,box=0
aia_lct,wave='171',/load
oplot,resp.logte,resp.all(*,2),thick=3,color=150
al_legend,'171A',linestyle=[0],color=[150],thick=3,position=[7.3,1.9e-24],linsize=0.5,box=0
aia_lct,wave='193',/load
oplot,resp.logte,resp.all(*,3),thick=3,color=150
al_legend,'193A',linestyle=[0],color=[150],thick=3,position=[7.3,1.2e-24],linsize=0.5,box=0
aia_lct,wave='211',/load
oplot,resp.logte,resp.all(*,4),thick=3,color=150
al_legend,'211A',linestyle=[0],color=[150],thick=3,position=[7.3,8e-25],linsize=0.5,box=0
aia_lct,wave='304',/load
oplot,resp.logte,resp.all(*,5),thick=3,color=150
al_legend,'304A',linestyle=[0],color=[150],thick=3,position=[7.3,5e-25],linsize=0.5,box=0
aia_lct,wave='335',/load
oplot,resp.logte,resp.all(*,6),thick=3,color=150
al_legend,'335A',linestyle=[0],color=[150],thick=3,position=[7.3,3e-25],linsize=0.5,box=0



device,/close
set_plot,'x'

END