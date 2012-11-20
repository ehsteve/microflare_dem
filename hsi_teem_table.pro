;; NAME:
; 		HSI_TEEM_TABLE
;
; PURPOSE:
; This procedure returns the best-fit DEM to a RHESSI count spectrum (at the Earth)
; using a differential emission measure distribution which has a either a Gaussian or Epstein dependence on temperature,
; for a range of values of peak temperature T and widths sigma. 
; 
; CATEGORY:
;       SPECTRA, XRAYS, DEM ANALYSIS
;
; CALLING SEQUENCE:
;       hsi_teem_table,epstein=epstein,n=n,spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time
;
; CALLS:
;       get_hsi_table_entry.pro
;	f_multi_therm.pro
;	f_multi_therm_gauss.pro
;	f_multi_therm_epstein.pro
;
; PREREQUISITES:
;	The SPEX and XRAY packages within SSW must be installed
;
; INPUTS:
;       NONE
;
; KEYWORD INPUTS:
;	EPSTEIN		: if set, the Epstein DEM profile is used. By default, the Gaussian profile is used.
;	N		: This needs to be set if the EPSTEIN keyword is set. N controls the steepness of the Epstein profile, such that
;			n=1 is the classical profile and n -> infinity corresponds to a boxcar function.
;	SPEC_FILE	: the RHESSI spectrum file to use for RHESSI DEM analysis, e.g. hsi_spectrum_20110621_182202_d1.fits
;	DRM_FILE	: the corresponding RHESSI DRM file to use for the RHESSI DEM analysis, e.g. hsi_srm_20110621_182202_d1.fits
;	FIT_TIME	: a 2-element string indicating the time interval between which the RHESSI DEM analysis should be performed, e.g.
;			fit_time = ['16-Jul-2011 17:02:00.000', '16-Jul-2011 17:03:00.000']
;	BKG_TIME	: a 2-element string indicating the time interval to use for the RHESSI spectral background subtraction,e.g.
;			bkg_time = ['16-Jul-2011 17:36:00.000', '16-Jul-2011 17:39:00.000']
;
; OUTPUTS:
;	HSI_TEEM_TABLE - a structure which is saved in the file 'hsi_teem_table.sav' or 'hsi_teem_table_epstein.sav'
;
; LIMITATIONS:
;	T and sig range are hardcoded.
;
; WRITTEN: Andrew Inglis, 2012/07/30
;	   Andrew Inglis, 2012/10/23 - modified to accept the input keywords EPSTEIN, N, SPEC_FILE, DRM_FILE, $
;				       FIT_TIME, and BKG_TIME. Now calls get_hsi_table_entry with the /EMFREE keyword.
;


PRO hsi_teem_table,epstein=epstein,n=n,spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time

;params used in get_hsi_table_entry has 6 (or 7 if /epstein) elements
;a0 - differential emission measure at the maximum of the Gaussian distribution (cm-3 kev-1)
;a1 - lower bound for the temperature integral, in keV
;a2 - upper bound for the temperature integral, in keV
;a3 - width of the DEM Gaussian, in log T (e.g. 0.3)
;a4 - value for the relative abundance of iron, nickel etc. (usually 1)
;a5 - location of the maximum of the Gaussian distribution, in keV

em_def=0.02
t_low = 0.09
t_high = 6.0
abun = 1.



t_min = 6.0
t_max = 7.5
t_d = 0.05
telog = t_d * findgen((t_max - t_min)/t_d) + t_min
tinkev = (10^telog)/11.6e6



tsig_min = 0.01
tsig_max = 0.80
tsig_d = 0.01
tsig = tsig_d * findgen((tsig_max - tsig_min)/tsig_d) + tsig_min

nsig=n_elements(tsig)
nte=n_elements(telog)

em_best_hsi=fltarr(nte,nsig)
chi_best_hsi=fltarr(nte,nsig)

		FOR k=0, nte-1 DO BEGIN
   			FOR l = 0, nsig-1 DO BEGIN

			IF keyword_set(epstein) THEN BEGIN
				get_hsi_table_entry,[em_def,t_low,t_high,tsig(l),tinkev(k),n,1],model_count_flux,real_count_flux,axis,summary, obj=obj, $
				spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time,epstein=epstein,/emfree
			ENDIF ELSE BEGIN
				get_hsi_table_entry,[em_def,t_low,t_high,tsig(l),tinkev(k),1],model_count_flux,real_count_flux,axis,summary, obj=obj, $
				spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time,/emfree
			ENDELSE

			;get_hsi_table_entry,[em_def,t_low,t_high,tsig(l),t(k),abun],model_count_flux,real_count_flux,axis,summary, obj=obj
			em_best_hsi[k,l] = summary.spex_summ_params[0]
			chi_best_hsi[k,l] = summary.spex_summ_chisq
			endfor
		endfor


chimin=MIN(chi_best_hsi,pos)
pos=array_indices(chi_best_hsi,pos)
telog_best=telog[pos[0]]
sig_best=tsig[pos[1]]

				


IF keyword_set(epstein) THEN BEGIN
	hsi_teem_table=create_struct('telog',telog,'tsig',tsig,'em_best_hsi',em_best_hsi,'chi_best_hsi',chi_best_hsi,'telog_best',telog_best,$
	'sig_best',sig_best,'chimin',chimin,'epstein',epstein,'n',n)
	SAVE,hsi_teem_table,filename='hsi_teem_table_epstein.sav',/verbose
ENDIF ELSE BEGIN
	hsi_teem_table=create_struct('telog',telog,'tsig',tsig,'em_best_hsi',em_best_hsi,'chi_best_hsi',chi_best_hsi,'telog_best',telog_best,'sig_best',sig_best,'chimin',chimin)
	SAVE,hsi_teem_table,filename='hsi_teem_table.sav',/verbose
ENDELSE

END