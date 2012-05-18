; OSPEX script created Fri Apr 27 16:42:22 2012 by OSPEX writescript method.                
;                                                                                           
;  Call this script with the keyword argument, obj=obj to return the                        
;  OSPEX object reference for use at the command line as well as in the GUI.                
;  For example:                                                                             
;     ospex_script_27_apr_2012, obj=obj                                                     
;                                                                                           
;  Note that this script simply sets parameters in the OSPEX object as they                 
;  were when you wrote the script, and optionally restores fit results.                     
;  To make OSPEX do anything in this script, you need to add some action commands.          
;  For instance, the command                                                                
;     obj -> dofit, /all                                                                    
;  would tell OSPEX to do fits in all your fit time intervals.                              
;  See the OSPEX methods section in the OSPEX documentation at                              
;  http://hesperia.gsfc.nasa.gov/ssw/packages/spex/doc/ospex_explanation.htm                
;  for a complete list of methods and their arguments.                                      
;                                                                                           
pro ospex_script_27_apr_2012, obj=obj                                                       
if not is_class(obj,'SPEX',/quiet) then obj = ospex()                                       
obj-> set, $                                                                                
 spex_specfile= '/Users/schriste/Desktop/flare_0_rhessi/hsi_spectrum_20110603_064600.fits'  
obj-> set, $                                                                                
 spex_drmfile= '/Users/schriste/Desktop/flare_0_rhessi/hsi_srm_20110603_064600.fits'        
obj-> set, spex_source_angle= 22.2012                                                       
obj-> set, spex_source_xy= [19.0353, -362.247]                                              
obj-> set, spex_fit_time_interval= [' 3-Jun-2011 07:15:48.000', $                           
 ' 3-Jun-2011 07:16:56.000']                                                                
obj-> set, spex_bk_time_interval=[' 3-Jun-2011 07:10:24.000', ' 3-Jun-2011 07:13:32.000']   
obj-> set, fit_function= 'vth+thick2'                                                       
obj-> set, fit_comp_params= [0.00100097, 0.965864, 1.00000, 0.488585, 6.03198, 100000., $   
 11.1289, 8.74125, 32000.0]                                                                 
obj-> set, fit_comp_minima= [1.00000e-20, 0.500000, 0.0100000, 1.00000e-10, 1.10000, $      
 1.00000, 1.10000, 1.00000, 100.000]                                                        
obj-> set, fit_comp_maxima= [1.00000e+20, 8.00000, 10.0000, 1.00000e+10, 20.0000, 100000., $
 20.0000, 1000.00, 1.00000e+07]                                                             
obj-> set, fit_comp_free_mask= [1B, 1B, 0B, 1B, 1B, 0B, 1B, 1B, 0B]                         
obj-> set, fit_comp_spectrum= ['full', '']                                                  
obj-> set, fit_comp_model= ['chianti', '']                                                  
obj-> set, spex_autoplot_units= 'Flux'                                                      
obj-> set, spex_eband= [[3.00000, 6.00000], [6.00000, 12.0000], [12.0000, 25.0000], $       
 [25.0000, 50.0000], [50.0000, 100.000], [100.000, 300.000]]                                
obj-> set, spex_tband= [[' 3-Jun-2011 06:46:00.000', ' 3-Jun-2011 07:00:47.000'], $         
 [' 3-Jun-2011 07:00:47.000', ' 3-Jun-2011 07:15:34.000'], [' 3-Jun-2011 07:15:34.000', $   
 ' 3-Jun-2011 07:30:21.000'], [' 3-Jun-2011 07:30:21.000', ' 3-Jun-2011 07:45:08.000']]     
end                                                                                         
