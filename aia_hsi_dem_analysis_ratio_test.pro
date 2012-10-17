PRO aia_hsi_dem_analysis_ratio_test, DIR = dir, HSI_IMAGE = hsi_image, FILESET = fileset, FORCE_TABLE = force_table

;first section should be to get the AIA fluxes based on the RHESSI image. Have some code that already does this.

;restore,'hsi_teem_table.sav',/verbose
;h=hsi_teem_table
wave_ =['131','171','193','211','335','94'] 
nwave =n_elements(wave_) 
nfiles = fltarr(nwave)

teem_table='teem_table.sav'


file_list = get_aia_file_list(dir, fileset = fileset)
FOR i = 0, nwave-1 DO nfiles[i] = n_elements(file_list[*,i])

;hardcoded this for 21 June 2011 flare
marker=120;120

;area of 1 AIA pixel in cm^2
aia_pixel_area = 1.85589e+15


t_min = 5.5
t_max = 8.0

t_min = 6.0
t_max = 7.5
t_d = 0.05
telog = t_d * findgen((t_max - t_min)/t_d) + t_min

tsig_min = 0.01
tsig_max = 0.40
tsig_d = 0.005
tsig = tsig_d * findgen((tsig_max - tsig_min)/tsig_d) + tsig_min

f = file_search(teem_table)

area = aia_teem_pixel_area(file_list[0,0])

IF f[0] EQ '' OR keyword_set(FORCE_TABLE) THEN aia_teem_table, wave_, tsig, telog = telog, q94 = q94, teem_table, save_dir = save_dir, area = area


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
		flux_obs = flux_
		counts = flux_obs*texp_
		;stop
		;when we sum pixels up the noise level gets very low - but not realistic, due to uncertainty in T response functions. Set a baseline uncertainty level.
		noise = 0.1*flux_obs;sqrt(counts)/texp_

		chimin = 9999.
		chi6min = 9999.
		chimin_hsi=9999.
		flux_dem_best = 9999.
		chi2d = fltarr(nte,nsig)
		chi2d_hsi = fltarr(nte,nsig)
		em_2d = fltarr(nte,nsig)
		em_2d_hsi = fltarr(nte,nsig)
		model_count_flux_hsi = fltarr(nte,nsig,301)
		real_count_flux_hsi = fltarr(nte,nsig,301)
		ratio_test=fltarr(nte,nsig,6)		

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
				chi	= sqrt(total((flux_obs-flux_dem)^2/noise^2)/(nwave-nfree))
				chi6 = abs(flux_obs-flux_dem)/noise
				chi2d[k,l] = chi
				em_2d[k,l] = em1
				
				ratio_test[k,l,*]=flux_obs/flux_dem
				
				;add in RHESSI
				;em_hsi=em1 /1e25
				;em_hsi=em_hsi/1e24    ;
				;em_hsi=em_hsi*flare_area;39.*2e15 
				;IF (k eq 0) THEN BEGIN
				;em_hsi=em_hsi * (10^telog(1) - 10^telog(0))
				;ENDIF ELSE BEGIN
				;em_hsi=em_hsi * (10^telog(k) - 10^telog(k-1))
				;ENDELSE	
				;em_2d_hsi[k,l]=h.em_best_hsi[k,l]
				;get_hsi_table_entry,[em_hsi,telog(k),tsig(l)],model_count_flux,real_count_flux,axis,summary, obj=obj
				
				;if the model function is all zeros then deal with this
				;typ=datatype(model_count_flux)
				;find the data type of model_count_flux. If there is no model flux at all then an integer of -1 will be
				;returned in place of the usual structure. 
				;IF (typ eq 'INT') then begin
				;	model_count_flux_hsi[k,l,*] = 0.
				;	real_count_flux_hsi[k,l,*] = real_count_flux
				;	chi2d_hsi[k,l] = 9999.
				;ENDIF ELSE BEGIN
				;	model_count_flux_hsi[k,l,*] = model_count_flux.yvals
				;	real_count_flux_hsi[k,l,*] = real_count_flux
				;	chi2d_hsi[k,l] = summary.spex_summ_chisq		
				;ENDELSE
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
				;IF (chi2d_hsi[k,l] le chimin_hsi) THEN BEGIN
				;	chimin_hsi=chi2d_hsi[k,l]
				;	em_best_hsi=alog10(em1)
				;	telog_best_hsi = telog[k]
				;	sig_best_hsi = tsig[l]
				;	model_count_flux_hsi_best=model_count_flux_hsi[k,l,*]
				;ENDIF
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
			;leg = num2str(chimin_hsi) + '[' + num2str(telog_best_hsi) + ',' + num2str(sig_best_hsi) + ']'
			;legend, leg, psym = 4
			;contour,chi_hsi,telog, tsig, levels=[chimin_hsi + 2.3],/over,thick=2, color = 6
		;ENDIF

;SAVE,chi2d,chi_hsi,axis,em_2d,em_2d_hsi,real_count_flux_hsi,model_count_flux_hsi,flux_dem_best,flux_dem,em_best,em_best_hsi,telog_best,$
;telog_best_hsi,sig_best,sig_best_hsi,chimin,chimin_hsi,flare_area,filename='aia_hsi_fit_results_bk.sav',/verbose
;		;stop

aia_hsi_fit_results=create_struct('chi_2d',chi2d,'chi2d_hsi',chi2d_hsi,'em_2d',em_2d,'em_2d_hsi',em_2d_hsi,$
$
'flux_dem_best',flux_dem_best,'flux_dem',flux_dem,'em_best',em_best,'telog_best',telog_best,$
'sig_best',sig_best,'chimin',chimin,'chimin_hsi',chimin_hsi,'flare_area',flare_area,'ratio_test',ratio_test)

SAVE,aia_hsi_fit_results,filename='aia_hsi_fit_results_ratio_test.sav',/verbose

SAVE,wave_,flux_obs,num_aia_pixels,filename='aia_fluxes_from_20110621_182200.sav'



END