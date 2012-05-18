PRO aia_teem_integrate, dir, HSI_FITS = hsi_fits, time_range = time_range

; PURPOSE: Integrate (avearge) a set of DEM measurements from AIA (to match the integration
; 			time of RHESSI for example)
;
; KEYWORDS: Give it a time range of a hsi_fits file which will be used to figure out the
;			time range.
;
; WRITTEN: Steven Christe (20-Dec-2011)

default, dir, '~/Desktop/' + 'flare0_deem/'
default, hsi_fits, '/Users/schriste/Desktop/flare_0_rhessi/ospex_results_20110603_071556_071656_mtherm.fits'

IF keyword_set(HSI_FITS) THEN BEGIN
	result = spex_read_fit_results(hsi_fits)
	p = result.spex_summ_params
		
	;now find the fit that is closest in time to time t
	time_range = result.spex_summ_time_interval
ENDIF	

f = file_search(dir + 'teem_tot_*.sav')

fname = filename(f)

times = unbreak_time(strmid(fname, 12,18))