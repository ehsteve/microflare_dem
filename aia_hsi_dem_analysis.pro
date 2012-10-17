;
; NAME:
; 		AIA_HSI_DEM_ANALYSIS
;
; PURPOSE:
; This procedure performs differential emission measure (DEM) analysis of both AIA
; and RHESSI data, using either a simple Gaussian distribution as the DEM model, or
; an Epstein function as the DEM model.
;
; CATEGORY:
;       DEM ANALYSIS
;
; CALLING SEQUENCE:
;       aia_hsi_dem_analysis, DIR = dir, HSI_IMAGE = hsi_image, FILESET = fileset, FORCE_TABLE = force_table, aia_only=aia_only,epstein=epstein,n=n, $
;       spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time
;
; CALLED BY: do_dem_analysis
;
; CALLS:
;	get_aia_file_list.pro
;	aia_teem_pixel_area.pro
;       aia_teem_table.pro
;	aia_teem_table_epstein.pro
;	get_hsi_table_entry.pro
;
; PREREQUISITES:
;	The SPEX and XRAY packages within SSW must be installed
;	
;
; INPUTS:
;	None - all keywords
;
; KEYWORD INPUTS:
;	DIR 		: directory to look for AIA fits files. Default is the current directory
;	HSI_IMAGE 	: filename of a RHESSI image to be used in the DEM analysis. The 50% contour from the RHESSI image is used
;			to define the region of interest over which the DEM is calculated in the AIA images.
;	FILESET		: the prefix describing the AIA files to search for. The options are 'AIA' (finds all files in AIA20110101* format)
;			 or 'ssw_cutout' (finds all files in ssw_cutout* format). Currently only AIA is recommended.
;	FORCE_TABLE	: forces aia_teem_table.pro to be re-run, and the AIA DEM lookup table to be recalculated.
;	AIA_ONLY	: if set, DEM analysis is only performed on AIA - RHESSI data is ignored.
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
;
; OUTPUTS:
; Depending on keywords, will produce some of the following:
; 	AIA_FIT_RESULTS.SAV
;	AIA_HSI_FIT_RESULTS.SAV
;	AIA_FIT_RESULTS_EPSTEIN.SAV
;	AIA_HSI_FIT_RESULTS_EPSTEIN.SAV
;
; PROCEDURE:
; 1) For an event, download to a directory a set of AIA data files of the 6 optically thin wavelengths (94A, 131A, 171A, 193A, 211A, 335A)
; 2) Create a RHESSI image of the event.
; 3) For HSI analysis, create a set of RHESSI spectrum files. Use seperate detectors, native energy bins.
; 4) Execute aia_hsi_dem_analysis, or alternatively call it with the DO_DEM_ANALYSIS wrapper.
;
; WRITTEN: Andrew Inglis, 2012/10/17
;




PRO aia_hsi_dem_analysis, DIR = dir, HSI_IMAGE = hsi_image, FILESET = fileset, FORCE_TABLE = force_table, aia_only=aia_only,epstein=epstein,n=n, $
spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time


IF (n_elements(fit_time) ne 2) THEN BEGIN
	print,'fit_time must be a 2-element string. Aborting'
	return
ENDIF ELSE IF (n_elements(bkg_time) ne 2) THEN BEGIN
	print,'bkg_time must be a 2-element string. Aborting'
	return
ENDIF ELSE IF NOT keyword_set(spec_file) THEN BEGIN
	print,'No RHESSI spectrum file selected. Aborting'
	return
ENDIF ELSE IF NOT keyword_set(drm_file) THEN BEGIN
	print,'No RHESSI DRM file selected. Aborting'
	return
ENDIF ELSE BEGIN
	print,' '
	print,'----------------------------------------------------------'
	print,'About to run DEM analysis with the following parameters:'
	print,' '
	print,'DIR: ',dir
	IF keyword_set(force_table) THEN print,'FORCE TABLE'
	IF keyword_set(aia_only) THEN print,'AIA ONLY'
	IF keyword_set(epstein) THEN print,'EPSTEIN'
	IF keyword_set(n) THEN print,'N = ',n
	print,'RHESSI image: ',hsi_image
	print,'RHESSI fit time: ',fit_time
	print,'RHESSI background time: ',bkg_time
	print,'RHESSI spectrum file: ',spec_file
	print,'RHESSI DRM file: ',drm_file
	print,'----------------------------------------------------------'
	print,' '
	wait,5
ENDELSE




;first section should be to get the AIA fluxes based on the RHESSI image. Have some code that already does this.

wave_ =['131','171','193','211','335','94'] 
nwave =n_elements(wave_) 
nfiles = fltarr(nwave)

IF keyword_set(epstein) THEN BEGIN
teem_table='teem_table_epstein.sav'
ENDIF ELSE BEGIN
teem_table='teem_table.sav'
ENDELSE


file_list = get_aia_file_list(dir, fileset = fileset)
FOR i = 0, nwave-1 DO nfiles[i] = n_elements(file_list[*,i])

;hardcoded this for 21 June 2011 flare
;marker=120;120
;hardcode for the 15 July 2011 flare
;marker=122

aia_time=anytim(fit_time[0],/utime)

filetimes=anytim(aiaprep_to_time(file_list[*,0]),/utime)
marker=value_locate(filetimes,aia_time)

print,'AIA selected time is: ',anytim(filetimes(marker),/vms)
;print,marker
;stop
;area of 1 AIA pixel in cm^2
aia_pixel_area = 1.85589e+15


t_min = 6.0
t_max = 7.5
t_d = 0.05
telog = t_d * findgen((t_max - t_min)/t_d) + t_min

tsig_min = 0.05
tsig_max = 0.80
tsig_d = 0.01
tsig = tsig_d * findgen((tsig_max - tsig_min)/tsig_d) + tsig_min

f = file_search(teem_table)

area = aia_teem_pixel_area(file_list[0,0])

IF keyword_set(epstein) THEN BEGIN
IF f[0] EQ '' OR keyword_set(FORCE_TABLE) THEN aia_teem_table_epstein, wave_, tsig, telog = telog, q94 = q94, teem_table, save_dir = save_dir, area = area, n=n
ENDIF ELSE BEGIN
IF f[0] EQ '' OR keyword_set(FORCE_TABLE) THEN aia_teem_table, wave_, tsig, telog = telog, q94 = q94, teem_table, save_dir = save_dir, area = area
ENDELSE



; if a hsi_image was given then create a mask out of it
IF hsi_image[0] NE '' THEN BEGIN
	fits2map, file_list[marker,0], aiamap
	fits2map, hsi_image, hsimap
	; interpolate the rhessi map to the aia map
	
	mask_map = inter_map(hsimap,aiamap)
	mask_map = drot_map(mask_map, time = aiamap.time)
	m = max(mask_map.data)
	; set the mask at everything above 50% contour
	index = where(mask_map.data LE m*0.5, complement = complement)
	mask_map.data[index] = 0
	mask_map.data[complement] = 1
	
	; now define the inverse mask
	invmask_map = mask_map
	invmask_map.data[index] = 1
	invmask_map.data[complement] = 0
ENDIF

num_aia_pixels=total(mask_map.data)
print,num_aia_pixels

flare_area=aia_pixel_area*num_aia_pixels
;stop

;;;;;;



default, save_dir, ''

nwave = n_elements(wave_)
flux_ = fltarr(nwave)
texp_ = fltarr(nwave)

file_iw = reform(file_list[marker,*])

;file_bk_iw = reform(file_list(marker-100,*))

FOR iw = 0, nwave-1 DO BEGIN
	 
	read_sdo,file_iw[iw],index,data
;	read_sdo,file_bk_iw[iw],index_bk,data_bk

        index2map,index,float(data),map
;	index2map,index_bk,float(data_bk),map_bk
			
	IF keyword_set(xrange) AND keyword_set(yrange) THEN BEGIN
	    sub_map, map, smap, xrange = xrange, yrange = yrange
;	    sub_map, map_bk,smap_bk,xrange=xrange,yrange=yrange
	    map = smap
;	    map_bk = smap_bk
	    data = smap.data
 ;           data_bk = smap_bk.data
    ENDIF
    
	s = size(data)
	i1=0 & j1=0 & i2=s[1]-1 & j2=s[2]-1

	image = data
	texp = map.dur
	texp_[iw]=texp 
	dateobs = map.time	
	
	IF keyword_set(mask_map) then begin
		mask = mask_map.data
		; zero out everything that is not in the mask
		FOR k = 0, i2 DO BEGIN
			FOR l = 0, j2 DO BEGIN
				data[k,l] = data[k,l]*mask[k,l]
;				data_bk[k,l] = data_bk[k,l] * mask[k,l]
			ENDFOR
		ENDFOR	
	ENDIF
	
	flux_[iw] = total(data[i1:i2,j1:j2])/texp; - (total(data_bk[i1:i2,j1:j2])/map_bk.dur)
	print,'Total flux in ', wave_(iw),flux_(iw)
ENDFOR

;stop

;then do the mapping using the AIA fluxes. Since RHESSI is basically 1 pixel for spectrum we only need one summed AIA pixel too.




restore, save_dir + teem_table, /verbose

dim	= size(flux)
nte	= dim[1]
nsig = dim[2]
nwave = dim[3]
ntot = nte*nsig

nfree=3


;flux_obs = reform(images[i,j,*])
		flux_obs = flux_ / num_aia_pixels
		counts = flux_obs*texp_
		;stop
		;when we sum pixels up the noise level gets very low - but not realistic, due to uncertainty in T response functions. Set a baseline uncertainty level.
		noise = 0.2*flux_obs;sqrt(counts)/texp_

		chimin = 9999.
		chi6min = 9999.
		chimin_hsi=9999.
		flux_dem_best = 9999.
		chi2d = fltarr(nte,nsig)
		chi2d_hsi = fltarr(nte,nsig)
		em_2d = fltarr(nte,nsig)
		em_2d_hsi = fltarr(nte,nsig)
		flux_dem_3d = fltarr(nte,nsig,6)
		model_count_flux_hsi = fltarr(nte,nsig,300)
		real_count_flux_hsi = fltarr(nte,nsig,300)
		
		telog_err = [0,0]
		sig_err = [0,0]
		sig_errsymmetric = 0.
		telog_errsymmetric = 0.
		
		FOR k=0, nte-1 DO BEGIN
   			FOR l = 0, nsig-1 DO BEGIN
				flux_dem1 = reform(flux[k,l,*])
				em1	= total(flux_obs)/total(flux_dem1) 
				em1_err = sqrt(total(noise^2)) / total(flux_dem1)
				flux_dem = flux_dem1 * em1
				flux_dem_3d[k,l,*]=flux_dem
				chi	= sqrt(total((flux_obs-flux_dem)^2/noise^2)/(nwave-nfree))
				chi6 = abs(flux_obs-flux_dem)/noise
				chi2d[k,l] = chi
				em_2d[k,l] = em1
				;add in RHESSI
				em_hsi=em1 /1e25
				em_hsi=em_hsi/1e24    ;
				em_hsi=em_hsi*flare_area;39.*2e15 
				em_hsi=em_hsi*11.6e6 ; convert to cm^-3 kev^-1
				;IF (k eq 0) THEN BEGIN
				;em_hsi=em_hsi * (10^telog(1) - 10^telog(0))
				;ENDIF ELSE BEGIN
				;em_hsi=em_hsi * (10^telog(k) - 10^telog(k-1))
				;ENDELSE
				tinkev=(10^telog(k)) / 11.6e6	
				em_2d_hsi[k,l]=em_hsi
				IF NOT keyword_set(aia_only) THEN BEGIN
					IF keyword_set(epstein) THEN BEGIN
						;get_hsi_table_entry_epstein,[em_hsi,0.09,8.0,tsig(l),tinkev,10,1],model_count_flux,real_count_flux,axis,summary, obj=obj, $
						;spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time
						get_hsi_table_entry,[em_hsi,0.09,8.0,tsig(l),tinkev,n,1],model_count_flux,real_count_flux,axis,summary, obj=obj, $
						spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time,epstein=epstein
					ENDIF ELSE BEGIN
						get_hsi_table_entry,[em_hsi,0.09,6.0,tsig(l),tinkev,1],model_count_flux,real_count_flux,axis,summary, obj=obj, $
						spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time
					ENDELSE
					;if the model function is all zeros then deal with this
					typ=datatype(model_count_flux)
					;find the data type of model_count_flux. If there is no model flux at all then an integer of -1 will be
					;returned in place of the usual structure. 
					IF (typ eq 'INT') then begin
						model_count_flux_hsi[k,l,*] = 0.
						real_count_flux_hsi[k,l,*] = real_count_flux
						chi2d_hsi[k,l] = 9999.
					ENDIF ELSE BEGIN
						model_count_flux_hsi[k,l,*] = model_count_flux.yvals
						real_count_flux_hsi[k,l,*] = real_count_flux
						chi2d_hsi[k,l] = summary.spex_summ_chisq		
					ENDELSE
				ENDIF
				;stop
				IF (chi le chimin) THEN BEGIN
					chimin = chi
					chi6min	= chi6
					em_best	= alog10(em1)
					em_best_err = em1_err/em_best
					telog_best	= telog[k]
					sig_best = tsig[l]
					flux_dem_best = flux_dem
					
				ENDIF
				IF NOT keyword_set(aia_only) THEN BEGIN
					IF (chi2d_hsi[k,l] le chimin_hsi) THEN BEGIN
						chimin_hsi=chi2d_hsi[k,l]
						em_best_hsi=em_2d_hsi[k,l]
						telog_best_hsi = telog[k]
						sig_best_hsi = tsig[l]
						model_count_flux_hsi_best=model_count_flux_hsi[k,l,*]
					ENDIF
				ENDIF
			ENDFOR 
		ENDFOR

; find errors in tsig and telog by finding the 1 sigma contour in chi^2 space for each pixel - A. Inglis 5-Oct-2011
		IF (chimin LT 9999.) THEN BEGIN

			;find the 1 sigma contour in chi^2-space to get errors in tsig and telog
			contour, chi2d, telog, tsig, levels = [chimin +2.3], path_xy=path_xy, path_info=path_info, /path_data_coords, closed=0
			IF exist(path_xy) THEN BEGIN
				telog_err = [telog_best - min(path_xy[0,*]), max(path_xy[0,*]) - telog_best]
				sig_err = [sig_best - min(path_xy[1,*]), max(path_xy[1,*]) - sig_best]

				;error bars are asymmetric - symmetrise!
				telog_errsymmetric = max(abs(telog_err))	
				sig_errsymmetric = max(abs(sig_err))
			ENDIF 
		ENDIF
		
		
		
			;window, 0
			;loadct, 0
			;hsi_linecolors
			;nlevels = 20
			;levels = chimin * (findgen(nlevels)*0.1 + 1.0)
			;anot = strarr(20)
			;FOR k = 0, nlevels-1 DO anot[k] = num2str(levels[k])
			;contour, chi2d, telog, tsig, levels = levels, c_annotation = anot, xtitle = 'log[Temp]', ytitle = 'sig'
			;oplot, [telog_best], [sig_best], psym = symcat(16), color = 6
			;leg = num2str(chimin) + '[' + num2str(telog_best) + ',' + num2str(sig_best) + ']'
			;legend, leg, psym = 4
			;contour,chi2d,telog, tsig, levels=[chimin + 2.3],/over,thick=2, color = 6
			
			;window, 1
			;emlog =em_best*exp(-(telog-telog_best)^2/(2.*sig_best^2))
			;plot, telog, emlog, yrange=minmax(emlog),xtitle='Temperature  log(T)',$
   			;ytitle='Emission measure  log(EM [cm!U-5!N K!U-1!N])'
   			;r = chianti_spec_from_dem(telog, emlog, /plot)
   			;save, emlog, telog, r, filename = 'xray_spectrum_i' + num2str(i) + '_j' + num2str(j) + '.sav'

			 ; plot the RHESSI chi^2 map 
			;hsi_levels = chimin_hsi * (findgen(nlevels)*0.1 + 1.0)
			;FOR k = 0, nlevels-1 DO anot[k] = num2str(levels[k])
			;window,5
			;contour,chi_hsi,telog,tsig, levels = hsi_levels, c_annotation = anot, xtitle = 'log[Temp]', ytitle = 'sig',title='RHESSI X^2 map'
			;oplot, [telog_best_hsi], [sig_best_hsi], psym = symcat(16), color = 6
			;
			;legend, leg, psym = 4
			;contour,chi_hsi,telog, tsig, levels=[chimin_hsi + 2.3],/over,thick=2, color = 6
		;ENDIF

;SAVE,chi2d,chi_hsi,axis,em_2d,em_2d_hsi,real_count_flux_hsi,model_count_flux_hsi,flux_dem_best,flux_dem,em_best,em_best_hsi,telog_best,$
;telog_best_hsi,sig_best,sig_best_hsi,chimin,chimin_hsi,flare_area,filename='aia_hsi_fit_results_bk.sav',/verbose
;		;stop

IF keyword_set(aia_only) THEN BEGIN
	aia_fit_results=create_struct('flux_obs',flux_obs,'flux_dem_3d',flux_dem_3d,'chi_2d',chi2d,'em_2d',em_2d,'em_2d_hsi',em_2d_hsi,$
	'flux_dem_best',flux_dem_best,'flux_dem',flux_dem,'em_best',em_best,'telog_best',telog_best,$
	'sig_best',sig_best,'chimin',chimin,'flare_area',flare_area,'telog',telog,'tsig',tsig)

	IF keyword_set(epstein) THEN BEGIN
		SAVE,aia_fit_results,filename='aia_fit_results_epstein.sav',/verbose

		print,' '
		print,'-----------------------------------'
		print,'Saved results in file: aia_fit_results_epstein.sav'
		print,'-----------------------------------'
		print,' '
	ENDIF ELSE BEGIN
		SAVE,aia_fit_results,filename='aia_fit_results.sav',/verbose

		print,' '
		print,'-----------------------------------'
		print,'Saved results in file: aia_fit_results.sav'
		print,'-----------------------------------'
		print,' '
	ENDELSE
ENDIF ELSE BEGIN

	aia_hsi_fit_results=create_struct('flux_obs',flux_obs,'flux_dem_3d',flux_dem_3d,'chi_2d',chi2d,'chi2d_hsi',chi2d_hsi,'axis',axis,'em_2d',em_2d,'em_2d_hsi',em_2d_hsi,$
	'real_count_flux_hsi',real_count_flux_hsi,'model_count_flux_hsi',model_count_flux_hsi,'model_count_flux_hsi_best',model_count_flux_hsi_best,$
	'flux_dem_best',flux_dem_best,'flux_dem',flux_dem,'em_best',em_best,'em_best_hsi',em_best_hsi,'telog_best',telog_best,'telog_best_hsi',telog_best_hsi,$
	'sig_best',sig_best,'sig_best_hsi',sig_best_hsi,'chimin',chimin,'chimin_hsi',chimin_hsi,'flare_area',flare_area,'telog',telog,'tsig',tsig)
	
	IF keyword_set(epstein) THEN BEGIN
		SAVE,aia_hsi_fit_results,filename='aia_hsi_fit_results_epstein.sav',/verbose
		
		print,' '
		print,'-----------------------------------'
		print,'Saved results in file: aia_hsi_fit_results_epstein.sav'
		print,'-----------------------------------'
		print,' '
	ENDIF ELSE BEGIN
		SAVE,aia_hsi_fit_results,filename='aia_hsi_fit_results.sav',/verbose

		print,' '
		print,'-----------------------------------'
		print,'Saved results in file: aia_hsi_fit_results.sav'
		print,'-----------------------------------'
		print,' '
	ENDELSE
ENDELSE







END