PRO do_dem_analysis,image_file=image_file,spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time

default,image_file,'hsi_image_20110716_170350.fits'
default,spec_file, 'hsi_spectrum_20110716_161036_d1.fits'
default,drm_file, 'hsi_srm_20110716_161036_d1.fits'
default,fit_time,['16-Jul-2011 17:02:00.000', '16-Jul-2011 17:03:00.000'] 
default,bkg_time,['16-Jul-2011 17:36:00.000', '16-Jul-2011 17:39:00.000']    


aia_hsi_dem_analysis,dir='',fileset='AIA',/force_table,hsi_image=image_file,spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time,/aia_only

aia_hsi_dem_analysis,dir='',fileset='AIA',/force_table,hsi_image=image_file,spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time,/epstein,/aia_only,n=10

aia_hsi_dem_analysis,dir='',fileset='AIA',/force_table,hsi_image=image_file,spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time

aia_hsi_dem_analysis,dir='',fileset='AIA',/force_table,hsi_image=image_file,spec_file=spec_file,drm_file=drm_file,fit_time=fit_time,bkg_time=bkg_time,/epstein,n=10


;add plotting routines here
plot_aia_flux_ratios
plot_aia_flux_ratios,/epstein




END

