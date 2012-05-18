FUNCTION aia_teem_pixel_area, file

dsun = 1.49598d13

; grab the pixel size from any of the files
read_sdo,file[0], index, data
cdelt1 = index.cdelt1
dsun   = index.dsun_obs
arcsec = 2.*!pi*dsun/(1.e3*360.*60.*60.)         ;1 arcsec in km
pix    = cdelt1*arcsec*1.e5                     ;pixel in cm
area   = pix^2                                  ;pixel area in cm^2

RETURN, area

END