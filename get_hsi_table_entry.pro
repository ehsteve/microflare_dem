;; NAME:
; 		GET_HSI_TABLE_ENTRY
;
; PURPOSE:
; This procedure returns the photon spectrum at the Earth, d(Flux)(eph)/dt,
; from a differential emission measure distribution which has a Gaussian dependence on temperature,
; for a given set of Gaussian parameters. 
; 
; CATEGORY:
;       SPECTRA, XRAYS
;
; CALLING SEQUENCE:
;       get_hsi_table_entry, params,model_count_flux,real_count_flux,axis,summary, obj=obj
;
; CALLS:
;       f_therm_dem_gauss.pro
;
; INPUTS:
;       params -   parameter array which should consist of:
;                  a[0] = emission measure (in 10^49 cm-3) at the maximum
;                  point in the Gaussian DEM distribution, i.e. at tmax
;                  a[1] = the temperature where the peak of the
;                  Gaussian DEM distribution is located, in units of
;                  log T (K). e.g. a[1] = 7.0 for 10 MK.
;                  a[2] = The width of the DEM Gaussian distribution
;                  in log T (K) units. e.g. a[2] = 0.3.
;
; OUTPUTS:
;	model_count_flux 	- the count RHESSI count flux resulting from the model DEM function
;	real_count_flux		- the actual RHESSI count flux observed (taken from the spectrum file used)
;	axis			- the array of count energy bins used.
;	summary			- a structure containing all of the summary data from the model fitting procedure.
;
; LIMITATIONS:
;	Requires a hardcoded spectrum file and srm file.
;
; WRITTEN: Andrew Inglis, 2012/07/30
;


                                    
;  params is of the form [a0,a1,a2,a3,a4,a5] where:
;  ao - emission measure at tmax, in 10^49 cm-3 keV-1
; 
;                                                                                       
pro get_hsi_table_entry, params,model_count_flux,real_count_flux,axis,summary, obj=obj,spec_file=spec_file, drm_file=drm_file, $
fit_time=fit_time,bkg_time=bkg_time, epstein=epstein                                        
if not is_class(obj,'SPEX',/quiet) then obj = ospex();/no_gui)                                     
;obj-> set, $                                                                              
 ;spex_specfile= 'hsi_spectrum_20110621_175120_d4.fits'    
 obj-> set, $                                                                                                                       
 spex_specfile= spec_file;'/home/ainglis/physics/steven_code/hsi_flare_20110716_170350/hsi_spectrum_20110716_161036_d1.fits'
obj-> set, spex_drmfile= drm_file;'hsi_srm_20110716_161036_d1.fits'
;obj-> set, spex_source_angle= 23.8350                                                     
;obj-> set, spex_source_xy= [284.317, 263.935]                                             
obj-> set, spex_erange= [5.0000000D, 12.000000D]                                          
obj-> set, spex_fit_time_interval= fit_time ;['16-Jul-2011 17:02:00.000', $                         
; '16-Jul-2011 17:03:00.000']   
obj-> set, spex_bk_time_interval=bkg_time ;['16-Jul-2011 17:36:00.000', '16-Jul-2011 17:39:00.000']                                                           
obj-> set, spex_uncert= 0.0200000

IF keyword_set(epstein) THEN BEGIN
	obj -> set, fit_function= 'multi_therm_epstein'
	obj-> set, fit_comp_params= [params[0], params[1], params[2],params[3],params[4],params[5],params[6]]                              
	obj-> set, fit_comp_minima= [1.00000e-12, 0.087, 8.5, 0.005, 0.09, 0.01, 1]                                   
	obj-> set, fit_comp_maxima= [1.00000e+2, 0.087, 8.5, 0.4, 6.0, 10., 100]                             
	obj-> set, fit_comp_free_mask= [0B, 0B, 0B, 0B, 0B, 0B, 0B]                                               

ENDIF ELSE BEGIN                                                          
	obj-> set, fit_function= 'multi_therm_gauss'
        obj-> set, fit_comp_params= [params[0], params[1], params[2],params[3],params[4],params[5]]                              
	obj-> set, fit_comp_minima= [1.00000e-12, 0.087, 8.5, 0.005, 1., 0.09]                                   
	obj-> set, fit_comp_maxima= [1.00000e+2, 0.087, 8.5, 0.4, 1., 6.0]                             
	obj-> set, fit_comp_free_mask= [0B, 0B, 0B, 0B, 0B, 0B]          
ENDELSE
                                     
obj-> set, fit_comp_spectrum= ''                                                          
obj-> set, fit_comp_model= ''                                                             
obj-> set, spex_autoplot_units= 'Flux'                                                    
obj-> set, spex_eband= [[3.00000, 6.00000], [6.00000, 12.0000], [12.0000, 25.0000], $     
 [25.0000, 50.0000], [50.0000, 100.000], [100.000, 300.000]]                              
obj-> set, spex_tband= fit_time;['16-Jul-2011 17:02:00.000', '16-Jul-2011 17:03:00.000']  
obj->set, spex_fit_manual=0, spex_autoplot_enable=0, spex_fitcomp_plot_resid=0, spex_fit_progbar=0        
set_logenv, 'OSPEX_NOINTERACTIVE', '1'                     
     obj -> dofit, /all
 
;add the commands to get the real and model count flux, plus the X^2
model_count_flux = obj -> calc_func_components(this_interval=0, /use_fitted, photons=0, spex_units='flux')
real_count_flux= obj -> calc_summ(item='data_count_flux');obj ->getdata(class='spex_data',spex_units='flux')
axis=obj ->getaxis(/ct_energy,/edges_2)

summary = obj->get(/spex_summ) 


;return this information
         
end                                                                                       
