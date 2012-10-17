;; NAME:
; 		HSI_TEEM_TABLE_EPSTEIN
;
; PURPOSE:
; This procedure returns the count spectrum by RHESSI (at the Earth)
; from a differential emission measure distribution which has a Gaussian dependence on temperature,
; for a range of values of peak temperature T and Gaussian width sigma. 
; 
; CATEGORY:
;       SPECTRA, XRAYS
;
; CALLING SEQUENCE:
;       hsi_teem_table
;
; CALLS:
;       get_hsi_table_entry.pro
;	f_therm_dem_gauss.pro
;
; INPUTS:
;       NONE
;
; OUTPUTS:
;	HSI_TEEM_TABLE - a structure which is saved in the file 'hsi_teem_table.sav'
;
; LIMITATIONS:
;	T and sig range are hardcoded.
;
; WRITTEN: Andrew Inglis, 2012/07/30
;


PRO hsi_teem_table_epstein

;params has 7 elements
;a0 - differential emission measure at the maximum of the Gaussian distribution (cm-3 kev-1)
;a1 - lower bound for the temperature integral, in keV
;a2 - upper bound for the temperature integral, in keV
;a3 - width of the DEM Gaussian, in log T (e.g. 0.3)
;a4 - location of the maximum of the Gaussian distribution, in keV
;a5 - steepness parameter
;a6 - value for the relative abundance of iron, nickel etc. (usually 1)
  

em_def=0.02
t_low = 0.09
t_high = 6.0
abun = 1.
n=10.


t_min = 6.0
t_max = 7.5
t_d = 0.025
telog = t_d * findgen((t_max - t_min)/t_d) + t_min
t = (10^telog)/11.6e6



tsig_min = 0.01
tsig_max = 0.40
tsig_d = 0.01
tsig = tsig_d * findgen((tsig_max - tsig_min)/tsig_d) + tsig_min

nsig=n_elements(tsig)
nte=n_elements(telog)

em_best_hsi=fltarr(nte,nsig)
chi_best_hsi=fltarr(nte,nsig)

		FOR k=0, nte-1 DO BEGIN
   			FOR l = 0, nsig-1 DO BEGIN
				;t_low2= (10^(telog(k) - tsig(l)))/11.6e6
				;t_high2=(10^(telog(k) + tsig(l)))/11.6e6
				;IF (t_low2 lt 0.087) THEN t_low2 = 0.087
				;IF (t_high2 gt 8.5) THEN t_high = 8.5
				get_hsi_table_entry_epstein,[em_def,t_low,t_high,tsig(l),t(k),n,abun],model_count_flux,real_count_flux,axis,summary, obj=obj
				em_best_hsi[k,l] = summary.spex_summ_params[0]
				chi_best_hsi[k,l] = summary.spex_summ_chisq
			endfor
		endfor
				
hsi_teem_table_epstein=create_struct('telog',telog,'tsig',tsig,'em_best_hsi',em_best_hsi,'chi_best_hsi',chi_best_hsi)

SAVE,hsi_teem_table_epstein,filename='hsi_teem_table_epstein.sav',/verbose


END