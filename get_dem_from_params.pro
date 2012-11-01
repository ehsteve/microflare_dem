;
; NAME: GET_DEM_FROM_PARAMS
;
; PURPOSE: Utility function which returns a model DEM function given a
;set of three input parameters. These parameters are expected in 'AIA
;UNITS' (see convert_dem_units.pro for more details).
;
; CALLING SEQUENCE:
; result = get_dem_from_params(t, params,n=n,epstein=epstein)
;
; INPUTS:
;      t:      an input temperature array
;      params: the DEM fit parameters. Assumed to consist of 3
;              components.
;              params[0] - emission measure value at the DEM centre
;              params[1] - temperature where the DEM centre is located
;              params[2] - width parameter associated with the DEM
;                          distribution
;
; KEYWORDS:
;      n:        The steepness parameter of the Epstein distribution. Needed
;                if EPSTEIN keyword is set.
;      EPSTEIN:  If set, the input DEM parameters are interpreted as
;                an Epstein profile. Otherwise a Gaussian profile is assumed.
;
; WRITTEN:
;      Andrew Inglis, 2012/11/01
;



FUNCTION get_dem_from_params, t, params,n=n,epstein=epstein

;assume params has three components
;params[0] - emission measure value
;params[1] - temperature value
;params[2] - width value

em=params[0]
Tc = params[1]
sigma = params[2]

IF keyword_set(epstein) THEN BEGIN
   dem= em * (1/cosh([  (tc-t)/(sigma)]^n))^2
ENDIF ELSE BEGIN
   dem =em * exp(-(tc - t)^2 / (2. * sigma^2))
ENDELSE

return,dem

END
