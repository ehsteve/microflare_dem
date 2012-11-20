;+
;
; NAME:
; 		F_MULTI_THERM_GAUSS
;
; PURPOSE:
; This function returns the photon spectrum seen at the Earth
; for a multithermal model (optically thin thermal bremsstrahlung function,
; normalized to the Sun-Earth distance)
; The differential emission measure has a Gaussian dependence on log temperature.
;
; CATEGORY:
;       SPECTRA, XRAYS
;
; CALLING SEQUENCE:
;       Flux = f_multi_therm_gauss( eph, a )
;
; CALLED BY: NONE
;
; CALLS:
;       f_multi_therm
;
; INPUTS:
; eph - energy vector in keV, 2N edges or N mean energies
; a -   model parameters as defined below:
;   a(0) = differential emission measure at the maximum point in the
;          Gaussian (i.e. the point defined in a[4]), in units of 10^49 cm^(-3) keV^(-1)
;   a(1) = minimum plasma temperature to use in the integration, in keV
;   a(2) = maximum plasma temperature to use in the integration, in keV
;   a(3) = Width of the differential emission measure Gaussian, in
;          units of *log K* (e.g. a[3] = 0.2.)
;   a(4) = Temperature at which the maximum of the Gaussian DEM is located, in keV
;   a(5)  Relative abundance for Iron and Nickel
;            Relative to coronal abundance for chianti
;            Relative to solar abundance for mewe
;           (unless user selects a different abundance table manually)
;
; KEYWORD INPUTS:
;   REL_ABUN - 2x2 array giving Fe, Ni abundance [ 26,x],[28,x] ],  If rel_abun keyword not used,
;     the value of x is taken from apar[5] to make 2x2 array.
;     If that's not there either, x is 1.
;
;   In _extra, can pass /full, /continuum, /lines, /mewe, /chianti keywords (see f_vth)
;
; OUTPUTS:
; Array of photon fluxes at photon energies determined by eph
;
; PROCEDURE:
; Thermal bremsstrahlung fluxes are computed by integrating the isothermal bremsstrahlung
; from plasma in the temperature range a(1) (t_min) to a(2) (t_max) with differential
; emission measure DEM(T).  The DEM and bremsstrahlung emissivity are provided in F_MULTI_THERM.
; The integration is performed using Gaussian Quadrature.  The Gaussian Quadrature abscissas
; and weights are provided by BRM_GAULEG.
;
; WRITTEN: Andrew Inglis, 2012/08/30
;
; NOTES: Intended to supercede original f_multi_therm_gauss.pro, which
;was not integrated with f_multi_therm_pow.pro and
;f_multi_therm_exp.pro. 

function f_multi_therm_gauss, eph, a, rel_abun=rel_abun, _extra=_extra

;f_multi_therm expects a[4], if it exists, to refer to abundances. To avoid this, we
;pass in only a[0:3] as parameters. The temperature associated with
;the Gaussian maximum is passed into f_multi_therm through the keyword
;t_gauss.

;T_gauss is stored by f_multi_therm in _extra until it is called by
;the f_gauss function (see f_multi_therm.pro)


b=[a[0:3],a[5]]
t_gauss=a[4]

;call f_multi_therm with the /gauss keyword to perform the main procedure
return, f_multi_therm(/gauss, eph, b, rel_abun=rel_abun,t_gauss=t_gauss,_extra=_extra)
end
