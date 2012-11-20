function f_mexp, t, alpha,_extra=_extra
return, exp((2.0-t)/alpha)
end

function f_mpow, t, alpha,_extra=_extra
return, (2.0/t)^(alpha)
end

function f_gauss,t,alpha,t_gauss=t_gauss,_extra=_extra
t_logk=alog10(t*1e6/0.08617)
tmax_logk=alog10(t_gauss*1e6/0.08617)
return, exp(-(t_logk-tmax_logk)^2/(2.*alpha)^2)
end

function f_epstein,t,alpha,t_epstein=t_epstein,n_epstein=n_epstein,_extra=_extra
t_logk=alog10(t*1e6/0.08617)
tmax_logk=alog10(t_epstein*1e6/0.08617)
n=n_epstein
return,(1/cosh([  (t_logk-tmax_logk)/(alpha)]^n))^2
end

function f_therm_dem, eph, t, emission, alpha, abun_params, $
	func=func, rel_abun=rel_abun, sep_abun=sep_abun, _extra=_extra

;return, f_vth( eph, [emission, 0., rel_abun[1,0]], $
;   multi_temp=reform(t[0,*]), _extra=_extra) * call_function( func,t,alpha)

if keyword_set(sep_abun) then return, f_vth_abun( eph, [emission, 0., abun_params], rel_abun=rel_abun, $
   multi_temp=reform(t[0,*]), _extra=_extra) * call_function( func,t,alpha,_extra=_extra) else $
     return, f_vth( eph, [emission, 0., abun_params], rel_abun=rel_abun, $
     multi_temp=reform(t[0,*]), _extra=_extra) * call_function( func,t,alpha,_extra=_extra)
     
   end

;+
;
; NAME:
; 		F_MULTI_THERM
;
; PURPOSE:
; This function returns the photon spectrum seen at the Earth
; for a multithermal model (optically thin thermal bremsstrahlung function,
; normalized to the Sun-Earth distance)
; The differential emission measure can have an exponential, power law
; or Gaussian dependence on temperature.
;
; CATEGORY:
;       SPECTRA, XRAYS
;
; CALLING SEQUENCE:
;       Flux = f_multi_therm( eph, a, pow=pow, rel_abun=rel_abun, _extra=_extra)
;
; CALLED BY: f_mth_exp_bpow
;            f_multi_therm_gauss
;            f_multi_therm_pow
;            f_multi_therm_exp
;            f_multi_therm_abun_pow
;            f_multi_therm_abun_exp
;
; CALLS:
;       f_therm_dem
;       brm_gauleg
;
; INPUTS:
; eph - energy vector in keV, 2N edges or N mean energies
; a -   model parameters defined below
;   a(0) = differential emission measure at T = 2 keV in units of 10^49 cm^(-3) keV^(-1)
;   a(1) = minimum plasma temperature in keV
;   a(2) = maximum plasma temperature in keV
;   a(3) = temperature scale length in keV for calculating the differential emission measure
;          at temperature T:  DEM(T) = a(0) * exp( (2. - T) / a(3) )
;   a(4)  Relative abundance for Iron and Nickel
;            Relative to coronal abundance for chianti
;            Relative to solar abundance for mewe
;           (unless user selects a different abundance table manually)
;   If rel_abun keyword is used, apar(4) is set to keyword value.
;   If a[4] is not passed in, the abundances are set to 1.
;
; KEYWORD INPUTS:
;	SEP_ABUN - if set, means we have separate controls for abundances, so call f_vth_abun, not f_vth
;   REL_ABUN - 2x2 array giving Fe, Ni abundance [ 26,x],[28,x] ],  If rel_abun keyword not used,
;     the value of x is taken from apar(4) to make 2x2 array.  If that's not there either, x is 1.
;   In _extra, can pass /full, /continuum, /lines, /mewe, /chianti keywords (see f_vth)
;  EXPONENTIAL  - if set use the exponential dem
;  POW - if set use the powerlaw dem, default is to use POW
;  GAUSS - if set use the Gaussian DEM as a function of log T (log E). Default is to use POW
;
; OUTPUTS:
; Array of photon fluxes at photon energies determined by eph
;
; PROCEDURE:
; Thermal bremsstrahlung fluxes are computed by integrating the isothermal bremsstrahlung
; from plasma in the temperature range a(1) (t_min) to a(2) (t_max) with differential
; emission measure DEM(T).  The DEM and bremsstrahlung emissivity are provided in F_THERM_DEM.
; The integration is performed using Gaussian Quadrature.  The Gaussian Quadrature abscissas
; and weights are provided by BRM_GAULEG.
;
; WRITTEN: Linhui Sui, 2003/08/28
;
; REVISED: Gordon Holman, 2003/09/04, Enhanced documentation.
; Linhui Sui, 2004/03/04: check a[0] before start computation. If a[0]= 0., then
;						  output zero fluxes at all energy bands
;
; Kim Tolbert, 2004/03/04 - added _extra so won't crash if keyword is used in call
; Kim Tolbert, 2006/03/21 - pass _extra through to f_therm_dem_exp
; Kim, 19-Apr-2006.  Now has 5 params - added abundance (for Fe,Ni) as 5th param.
; 18-jul-2006, richard.schwartz@gsfc.nasa.gov; changed () to [] as appropriate
;	l, l1, l2 to L, L1, L2 to prevent confusion with 1
;	forces all bins with hi edge below 10 keV to be computed each time,
;	this prevents bad interaction with chianti_kev_lines which attempts
;	to spread narrow lines out over adjacent bins.
;	Combined f_multi_therm_exp and f_multi_therm_pow into a single
;	procedure differentiated only by the dem function passed through
;	to the integrating routine
; 12-May-2008, Kim. Pass any extra elements in a as abundance parameters to f_therm_dem, and pass
;	rel_abun keyword through, instead of checking and setting abun values here (f_vth will do it)
;	Also, added sep_abun keyword.
; 29-Aug-2012, Andrew Inglis. Added f_gauss function to fit a DEM Gaussian as a function of log T.
;       All functions (f_gauss,f_mexp,f_mpow) now accept the _extra keyword. 
;       The top level routine f_multi_therm_gauss_beta.pro was created, which calls f_multi_therm 
;       with the /gauss keyword and an additional keyword t_gauss.
; 29-Aug-2012, Andrew Inglis. Updated header documentation.
;
;

;-

function f_multi_therm, eph, a, rel_abun=rel_abun, sep_abun=sep_abun, $
	exponential=exponential, $
	pow=pow, gauss=gauss, epstein=epstein, $
	_extra=_extra

emission = a[0]
t_min = a[1]
t_max = a[2]
alpha = a[3]
abun_params = n_elements(a) eq 4 ? 1. : a[4:*]
;abun = n_elements(a) eq 4 ? 1. : a[4]
;rel_abun = keyword_set(rel_abun) ? rel_abun : reform( [26,abun,28,abun], 2,2)

default, func, 'f_mpow'
func = keyword_set(exponential) ? 'f_mexp' : func
func = keyword_set(gauss) ? 'f_gauss' : func
func = keyword_set(epstein) ? 'f_epstein' : func

estep = (size(eph))[2]

;	Create arrays for integral sum

intsum = FltArr(estep)


;check whether a[0] eq 0, if so, return zero flux
if a[0] eq 0. then goto, conv


tlowarr = fltarr(estep) + t_min

;   Specify the maximum order of the Gaussian quadrature integration
;   before ending without convergence.  MAXFCN should be less than
;   2^(NLIM).

maxfcn = 2048

;   Specify the maximum relative error allowed for convergence

rerr = 0.001

L = Where(tlowarr)

;	Create arrays for integral sum and error flags.

intsum = FltArr(estep)

;	Maximum possible order of the Gaussian quadrature integration is 2^nlim.

nlim = 12
ires_start = 1
FOR ires = ires_start, nlim DO BEGIN

	npoint = 2L^ires
	IF (npoint GT maxfcn) THEN GOTO, noconv

	eph1 = eph[*, L]

	;	Create arrays for abscissas and weights.

	x = DblArr(N_Elements(L),npoint)
	w = DblArr(N_Elements(L),npoint)

	;	Use Brm_GauLeg to perform the outer integration.

	Brm_GauLeg, tlowarr[L], t_max, npoint, x, w

    ;   Save the last approximation to the integral.

	lastsum = intsum

    ;   Obtain the numerical approximation to the integral.

	intsum[L] = (w*f_therm_dem(eph1, x, emission, alpha, abun_params, sep_abun=sep_abun, $
	             rel_abun=rel_abun, _extra=_extra, func=func)) $
		         # (Make_Array(npoint, Value=1))

	;	Use last two calculations to check for convergence.

	L1 = Abs(intsum-lastsum)
	L2 = rerr*Abs(intsum)


    ;   Obtain the indices of those integrals that have not yet converged.

	WL = Where(L1 GT L2, nwl)

    ;  Return the computed fluxes when all of the integrals have converged.

	IF (nwl EQ 0) THEN GOTO, conv
	L  = Where((L1 gt L2) or (eph1[1,*] lt 10.))


ENDFOR


noconv:

	print, '!!!!!!!!!!!!!!! no convergence'

conv:

	mmflux = intsum

RETURN, mmflux

END



