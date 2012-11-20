PRO get_mask_script,mask_map,invmask_map

	fits2map, 'AIA20110621_182200_0171.fits', aiamap
	fits2map, 'hsi_image_20110621_182202.fits', hsimap
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

END