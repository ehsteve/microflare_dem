;
; NAME:
; 		PLOT_XRAY_SPECTRUM_FROM_DEM
;
; PURPOSE:
; This procedure plots the X-ray spectrum that would be expected from a given DEM. The DEM is provided as a set of fitted parameters
; which come from either a Gaussian or an Epstein fitting model. See aia_hsi_dem_analysis.pro for details of the model fits.
;
; CATEGORY:
;       X-RAYS, DEM ANALYSIS
;
; CALLING SEQUENCE:
;       plot_xray_spectrum_from_dem, em, tlog, sigma, n=n, flare_area=flare_area, epstein=epstein
;;
; CALLS:
;	f_multi_therm_epstein.pro
;	f_multi_therm_gauss.pro
;	f_multi_therm.pro
;
; PREREQUISITES:
;	The SPEX and XRAY packages within SSW must be installed
;	
;
; INPUTS:
;	EM	:	The log of the differential emission measure at the centre of the fitted DEM distribution
;	TLOG	:	the temperature (in log T) marking the centre point of the fitted DEM distribution
;	sigma	:	The width parameter (in log T units) of the fitted DEM distribution
;
; KEYWORD INPUTS:
;	EPSTEIN	:	if set, the input parameters are interpreted as an Epstein DEM profile. If this is not set, a Gaussian profile is assumed.
;	N	:	This needs to be set if the EPSTEIN keyword is set. N is the steepness parameter of the Epstein profile, such that
;			n=1 is the classical profile and n -> infinity corresponds to a boxcar function.
;	PS	:	if set, plot to postscript instead of the X window
;
; KEYWORD OUTPUTS:
;	EPH	:	Optionally, return the energy array
;	PHOTONS :	Optionally, return the photons array
;
; LIMITATION: The input parameter units are assumed to be in the 'AIA UNITS' form. See CONVERT_DEM_UNITS.PRO for details
;
; WRITTEN: Andrew Inglis, 2012/10/31
;

PRO plot_xray_spectrum_from_dem,em, tlog, sigma, n=n, flare_area=flare_area,epstein=epstein,ps=ps, eph=eph, photons=photons 


IF NOT keyword_set(flare_area) THEN BEGIN
print,' '
print,'-----------------------------------'
print,'Keyword FLARE_AREA not given. Aborting'
print,'-----------------------------------'
print,' '
return
ENDIF


;create a 2xN dummy energy array similar to the RHESSI native energy bins. Roughly between 1 - 100 keV with 0.35 keV binning.
eph = fltarr(2,300)

eph[0,*]=(findgen(300)*0.352982)  +1.07638
eph[1,*]=(findgen(300)*0.352982)  +(1.07638 + 0.352982)

;convert EM from cm^-5 K^-5 to 10^49 cm^-3 keV^-1
em=10^em
em=em/1e24
em=em/1e25
em=em*flare_area
em=em*11.6e6

;convert t from log K to kev
tinkev=(10^tlog)/11.6e6

;set integration bounds. If a steep distribution, set the upper bound for integration
lower_bound=0.087
IF (n gt 1) THEN BEGIN
upper_bound=(tinkev + sigma + 0.3)
ENDIF ELSE BEGIN
upper_bound=8.5
ENDELSE

;now plot the photon spectrum
IF keyword_set(epstein) THEN BEGIN
photons=f_multi_therm_epstein(eph,[em,lower_bound, upper_bound,sigma,tinkev,n,1],_extra=_extra)
ENDIF ELSE BEGIN
photons=f_multi_therm_gauss(eph,[em,lower_bound, upper_bound,sigma,tinkev,1],_extra=_extra)
ENDELSE

IF keyword_set(ps) THEN BEGIN
set_plot,'ps'
device,encaps=1,filename='xray_spectrum_from_dem.ps'
ENDIF

plot,eph[0,*],photons,thick=2,charsize=1.2,/xlog,/ylog,xrange=[1,30],/xstyle,xtitle='energy [keV]',ytitle='photons cm!U-2!N keV!U-1!N'

IF keyword_set(epstein) THEN BEGIN
ssw_legend,['EM = ' + num2str(em) + ' cm!U-3!N kev!U-1!N','T = ' + num2str(tinkev) + ' keV','sigma = ' + num2str(sigma), 'n = ' + num2str(n),$
 'flare area = ' + num2str(flare_area) + ' cm^2', 'Model: Epstein'],/right,charsize=1.2
ENDIF ELSE BEGIN
ssw_legend,['EM = ' + num2str(em) + ' cm^-3 kev^-1','T = ' + num2str(tinkev) + ' keV','sigma = ' + num2str(sigma),$
 'flare area = ' + num2str(flare_area) + ' cm^2', 'Model: Gaussian'],/right,charsize=1.2
ENDELSE

IF keyword_set(ps) THEN BEGIN
device,/close
set_plot,'x'
ENDIF


END