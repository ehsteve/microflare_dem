PRO flux_dem_test


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




END