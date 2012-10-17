function f_multi_therm_epstein_defaults
;+
; NAME:
;	F_MULTI_THERM_EPSTEIN_DEFAULTS
;
; PURPOSE: Function to return default values for
;   parameters, minimum and maximum range, and free parameter mask when
;   fitting to f_multi_therm_epstein.
;
; CALLING SEQUENCE: defaults = f_multi_therm_gauss_defaults()
;
; INPUTS:
;	None
; OUTPUTS:
;	Structure containing default values
;
; MODIFICATION HISTORY:
; Andrew Inglis, 18-Sep-2012, first created.

;
;-
;------------------------------------------------------------------------------

  defaults = { $
  fit_comp_params:           [1e-2,  0.087, 8.5, 0.30, 1.,    1,    1.], $
  fit_comp_minima:           [1e-12, 0.087, 8.5, 0.01, 0.087, 1,    0.01], $
  fit_comp_maxima:           [1e8,   0.087, 8.5, 1.00, 8.5,   100,  10.], $
  fit_comp_free_mask:        [1b,    0b,    0b,  1b,   1b,    0b,   0b] $
             }

  return, defaults
end
