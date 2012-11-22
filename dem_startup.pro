;Name:
;   	DEM_STARTUP
;PURPOSE:
;   	Utility to ensure that all necessary local modules are compiled and OSPEX_MODELS_DIR is set
;	correctly prior to initiation of joint DEM analysis
;			
;CATEGORY:
;   DEM_ANALYSIS, XRAYS
;
;CALLING SEQUENCE:
;	dem_startup,ospex_models_dir=ospex_models_dir
;
;INPUTS:
;	None
;
; KEYWORD INPUTS:
;	OSPEX_MODELS_DIR  - provides the directory to search for the fit_model_components.txt file that
;	is used in the joint RHESSI and AIA DEM analysis. If not set, the local directory is searched instead.
;	if the file sucessfully found in the directory, OSPEX_MODELS_DIR is set to this directory.
;
; WRITTEN:
;	Andrew Inglis, 2012/11/21


PRO dem_startup,ospex_models_dir=ospex_models_dir

;compile needed modules locally to override files present in OSPEX
;resolve_routine,'f_mexp',/is_function,/compile_full_file
;resolve_routine,'f_mpow',/is_function
;resolve_routine,'f_gauss',/is_function
;resolve_routine,'f_epstein',/is_function
resolve_routine,'f_multi_therm',/is_function,/compile_full_file
resolve_routine,'f_multi_therm_gauss',/is_function
resolve_routine,'f_multi_therm_gauss_defaults',/is_function
resolve_routine,'f_multi_therm_epstein',/is_function
resolve_routine,'f_multi_therm_epstein_defaults',/is_function


;need to set the environment variable OSPEX_MODELS_DIR to point to the version of
;fit_model_components.txt that comes with these codes. Otherwise Epstein analysis
;with RHESSI will fail.
IF NOT keyword_set(ospex_models_dir) THEN BEGIN
	spawn,'pwd',ospex_models_dir
	print,' '
	print,'---------------------------------------------------------'
	print,'ospex_models_dir was not passed in. Trying local directory.'
	print,'---------------------------------------------------------'
	print,' '
ENDIF
	
;search the local directory for the fit_model_components.txt file
spawn,'\ls ' + ospex_models_dir + '/fit_model_components.txt',model_file
IF keyword_set(model_file) THEN BEGIN
	setenv,'OSPEX_MODELS_DIR=' + ospex_models_dir
	print,' '
	print,'----------------------------------------------------------------------'
	print,'Notice: fit_model_components.txt is present in ' + ospex_models_dir
	print,'Notice: set environment variable OSPEX_MODELS_DIR to ' + ospex_models_dir
	print,'----------------------------------------------------------------------'
	print,' '
ENDIF ELSE BEGIN
	print,' '
	print,'--------------------------------------------------------------------'
	print,'Notice: fit_model_components.txt not found in ' + ospex_models_dir
	print,'Notice: environment variable OSPEX_MODELS_DIR has not been changed.'
	print,'Warning: Epstein analysis with RHESSI may not run correctly.'
	print,'--------------------------------------------------------------------'
	print,' '
ENDELSE		

END