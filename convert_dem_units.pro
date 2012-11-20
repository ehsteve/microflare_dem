;Name:
;   CONVERT_DEM_UNITS
;PURPOSE:
;   This utility function converts between two sets of DEM units, 'AIA units', and 'RHESSI units'. 
;   These units are considered to be:
;	'AIA UNITS': EM - log of emission measure in cm^-5 K^-1 units, e.g. EM = 21.5
;		     T  - temperature in log T (K), e.g. T = 6.3
;		   sigma - width parameter in log T (K) units, e.g. sigma = 0.2
;
;	'RHESSI UNITS':	EM - emission measure in cm^-3 keV^-1 units, e.g. EM = 0.02
;			T - temperature in keV, e.g. T = 0.086 for a 1MK plasma temperature
;		     sigma - width parameter in log T (K) units, e.g. sigma = 0.2
;			
;CATEGORY:
;   DEM_ANALYSIS, XRAYS
;
;CALLING SEQUENCE:
;	result=convert_dem_units(params,flare_area=flare_area,aia_to_hsi=aia_to_hsi,hsi_to_aia=hsi_to_aia)
;
;INPUTS:
;	PARAMS - consisting of three elements:
;		params[0] - an emission measure value
;		params[1] - a temperature value
;		params[2] - a width value
;
; KEYWORD INPUTS:
;	FLARE_AREA - required keyword. Gives the area measurement to be used in the unit conversion, in cm^2
;	AIA_TO_HSI - if set, assumes that PARAMS is in 'AIA units' and converts to 'RHESSI units'
;	HSI_TO_AIA - if set, assumes that PARAMS is in 'RHESSI units' and converts to 'AIA units'
;
; WRITTEN:
;	Andrew Inglis, 2012/10/31


FUNCTION convert_dem_units, params,flare_area=flare_area,aia_to_hsi=aia_to_hsi,hsi_to_aia=hsi_to_aia

;assume params has three values:
;a[0] = EM
;a[1] = T
;a[2] = sigma

;if no keywords set don't know which way to convert
IF NOT keyword_set(aia_to_hsi) and NOT keyword_set(hsi_to_aia) THEN BEGIN
print,''
print,'-------------------------------------------------------------------------------'
print,'Neither /AIA_TO_HSI or HSI_TO_AIA is set. Conversion is ambiguous. Returning -1'
print,'-------------------------------------------------------------------------------'
print,' '
return,-1
ENDIF

;can't convert both ways at once
IF keyword_set(aia_to_hsi) AND keyword_set(hsi_to_aia) THEN BEGIN
print,''
print,'-----------------------------------------------------------'
print,'/AIA_TO_HSI and HSI_TO_AIA cannot both be set! Returning -1'
print,'-----------------------------------------------------------'
print,' '
return,-1
ENDIF

IF NOT keyword_set(flare_area) THEN BEGIN
print,''
print,'-----------------------------------------------------------------'
print,'Keyword FLARE_AREA needs to be given for conversion. Returning -1'
print,'-----------------------------------------------------------------'
print,' '
return,-1
ENDIF

em=params[0]
t=params[1]
sigma=params[2]

IF keyword_set(aia_to_hsi) THEN BEGIN
;convert EM from cm^-5 K^-1 to 10^49 cm^-3 keV^-1
em=10^em
em=em/1e24
em=em/1e25
em=em*flare_area
em=em*11.6e6

;convert t from log K to kev
tinkev=(10^t)/11.6e6

params[0]=em
params[1]=tinkev
params[2]=sigma

return,params
ENDIF

IF keyword_set(hsi_to_aia) THEN BEGIN

;convert EM from 10^49 cm^-3 keV^-1 to cm^-5 K^-1
em=em/11.6e6
em=em/flare_area
em=em*1e24
em=em*1e25
em=alog10(em)

;convert from keV to log T
logt = alog10(t*11.6e6)

params[0]=em
params[1]=logt
params[2]=sigma

return,params
ENDIF

print,'Should never get here. Something probably went wrong!'
return,params

END