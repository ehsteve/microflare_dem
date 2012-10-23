DEM analysis codes
=====

These IDL routines are designed to perform Differential Emission Measure (DEM) analysis using data from the AIA and RHESSI instruments. IDL and SSW must be available in order to use these codes.

Usage
-----

To use the RHESSI analysis tools in this package, the XRAY and SPEX branches of SSW must be installed.

There are a few different main executables in the bundle. In general, the following data are required to use them:

 * A set of level 1.5 AIA data files in the 6 optically thin wavelengths (131A, 171A, 193A, 211A, 335A, 94A)
 * A RHESSI FITS image of the event
 * RHESSI spectral and SRM files for the time of interest

An example of using one of the main executables is below:

    aia_hsi_dem_analysis,dir='',fileset='AIA',/force_table,hsi_image='hsi_image_20110716_170350.fits',spec_file='hsi_spectrum_20110716_161036_d1.fits', $
    drm_file='hsi_srm_20110716_161036_d1.fits',fit_time=['16-Jul-2011 17:02:00.000', '16-Jul-2011 17:03:00.000'], $
    bkg_time=['16-Jul-2011 17:36:00.000', '16-Jul-2011 17:39:00.000']

This will attempt to perform a DEM analysis of AIA and RHESSI data at 17:02 UT on 2011-Jul-16, provided the files listed are present. For detailed installation instructions, see the individual file headers.

