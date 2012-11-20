; OSPEX script created Thu Jul 12 17:11:57 2012 by OSPEX writescript method.              
                                    
;  params is of the form [a0,a1,a2] where:
;  ao - emission measure at tmax, in 10^49 cm-3
;  a1 - tmax, in log T (K) units, e.g. a1 = 6.6
;  a2 - sigma, in log T (K) units, e..g a2 = 0.3
                                                                                       
pro get_hsi_table_entry_iso, params,model_count_flux,real_count_flux,axis,summary, obj=obj                                                     
if not is_class(obj,'SPEX',/quiet) then obj = ospex();/no_gui)                                     
;obj-> set, $                                                                              
 ;spex_specfile= 'hsi_spectrum_20110621_175120_d4.fits'    
 obj-> set, $                                                                                                                       
 spex_specfile= '/home/ainglis/physics/steven_code/hsi_flare_20110621_182202/hsi_spectrum_20110621_175120_d4_bksub_bksub_test.fits'
obj-> set, spex_drmfile= 'hsi_srm_20110621_175120_d4.fits'
obj-> set, spex_source_angle= 23.8350                                                     
obj-> set, spex_source_xy= [284.317, 263.935]                                             
obj-> set, spex_erange= [3.0000000D, 12.000000D]                                          
obj-> set, spex_fit_time_interval= ['21-Jun-2011 18:22:00.000', $                         
 '21-Jun-2011 18:23:00.000']                                                              
obj-> set, spex_uncert= 0.100000                                                          
obj-> set, fit_function= 'vth'                                                
obj-> set, fit_comp_params= [params[0], params[1], params[2]]                              
obj-> set, fit_comp_minima= [1.00000e-8, 5.50000, 1.000000]                                   
obj-> set, fit_comp_maxima= [1.00000e+6, 7.50000, 1.000000]                             
obj-> set, fit_comp_free_mask= [0B, 0B, 0B]                                               
obj-> set, fit_comp_spectrum= ''                                                          
obj-> set, fit_comp_model= ''                                                             
obj-> set, spex_autoplot_units= 'Flux'                                                    
obj-> set, spex_eband= [[3.00000, 6.00000], [6.00000, 12.0000], [12.0000, 25.0000], $     
 [25.0000, 50.0000], [50.0000, 100.000], [100.000, 300.000]]                              
obj-> set, spex_tband= ['21-Jun-2011 18:22:00.000', '21-Jun-2011 18:23:00.000']  
obj->set, spex_fit_manual=0, spex_autoplot_enable=0, spex_fitcomp_plot_resid=0, spex_fit_progbar=0        
set_logenv, 'OSPEX_NOINTERACTIVE', '1'                     
     obj -> dofit, /all
 
;add the commands to get the real and model count flux, plus the X^2
model_count_flux = obj -> calc_func_components(this_interval=0, /use_fitted, photons=0, spex_units='flux')
real_count_flux= obj -> calc_summ(item='data_count_flux');obj ->getdata(class='spex_data',spex_units='flux')
axis=obj ->getaxis(/ct_energy)

summary = obj->get(/spex_summ) 


;return this information
         
end                                                                                       
