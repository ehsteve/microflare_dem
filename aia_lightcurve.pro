FUNCTION aia_lightcurve, DIR = dir, HSI_FITS = hsi_fits, DEBUG = debug, PLOT = plot, ISTART = istart, IEND = iend, MASK_LEVEL = mask_level, WAVE_LIST = wave_list

; PURPOSE: Given a directory full of AIA images, generate a light curve
; of the integrated flux over the image.
;
; KEYWORDS: HSI_FITS - defines a mask. Only includes pixels above masklevel
;           MASK_LEVEL - defines the level (in percent of max) above which
;           pixels will be included if using a mask
;           ISTART - begin the light curve from the istart file (to speed up)
;           IEND - end the light curve at the iend file (to speed up)
;           WAVE_LIST - set the wavelength 
;
; WRITTEN: Steven Christe (8-Oct-2011)
; MODIFIED: Steven Christe (26-Oct-2011)

default, dir, '/Users/schriste/Desktop/flare_0_aia/'
;fileset = file_search(dir + 'ssw_cutout_20110716_*')
default, hsi_fits, '/Users/schriste/Desktop/flare_0_rhessi/hsi_image_20110716_170350.fits'
default, istart, 0
default, mask_level, 50

IF NOT KEYWORD_SET(WAVE_LIST) THEN wave_list =['131','171','193','211','335','94'] 
nwave = n_elements(wave_list)

file_list = get_aia_file_list(dir, wave_list = wave_list)
nfiles = n_elements(file_list[*,0])

IF keyword_set(HSI_FITS) THEN BEGIN
	fits2map, file_list[0,0], aiamap
	fits2map, hsi_fits, map
	; interpolate the rhessi map to the aia map
	imap = inter_map(map,aiamap)
	m = max(imap.data)
	; set the mask at everything above 50% contour
	index = where(imap.data LE m*mask_level/100.0, complement = complement)
	imap.data[index] = 0
	imap.data[complement] = 1
	mask = imap.data
	inverse_mask = NOT mask
ENDIF

result = dblarr(nwave,2,nfiles)

IF NOT keyword_set(iend) THEN iend = nfiles

FOR j = 0, nwave-1 DO BEGIN
    FOR i = istart, iend-1 DO BEGIN
        fits2map, file_list[i,j], aiamap    
	    result[j,0,i] = anytim(aiamap.time)
    	IF keyword_set(hsi_fits) THEN BEGIN
            imap_rot = drot_map(imap, ref_map = aiamap)
            aiamap.data = aiamap.data * imap_rot.data
        ENDIF	    
        result[j,1,i] = total(aiamap.data)
    ENDFOR
ENDFOR

IF keyword_set(PLOT) THEN BEGIN
	loadct,0
	hsi_linecolors
	P.MULTI = [0,1,nwave]
	charsize = 1.5
	FOR j = 0, nwave-1 DO BEGIN    
        yrange = minmax(result[j,1,*])
        utplot, anytim(result[j,0,*],/yoh), result[j,1,*], /nodata, yrange = yrange, charsize = charsize
        outplot, anytim(result[j,0,*],/yoh), result[j,1,*], color = 6, thick = 2.0
        legend, wave_list[j] + ' A', /right, charsize = charsize, color = 6
    ENDFOR
ENDIF

IF keyword_set(DEBUG) THEN stop

END