
;standard AR

xrange = [-28.5000, 571.500]
yrange = [-508.5, 121.500]

aia_teem_movie, xrange = xrange, yrange = yrange, fileset = 'AIA', npix = 1, macro_dem = 10

restore, '/Users/schriste/Dropbox/idl/aia_deem/teem_data_00020110215_000009_npix1.sav', /verbose

wave_ =['131','171','193','211','335','94'] 

aia_lct, rr, gg, bb, wavelnth=wave_[0], /load
plot_map, AIA_MAP_CUBE[0], /limb, /log

aia_lct, rr, gg, bb, wavelnth=wave_[1], /load
plot_map, AIA_MAP_CUBE[1], /limb, /log

aia_lct, rr, gg, bb, wavelnth=wave_[2], /load
plot_map, AIA_MAP_CUBE[2], /limb, /log

aia_lct, rr, gg, bb, wavelnth=wave_[3], /load
plot_map, AIA_MAP_CUBE[3], /limb, /log

aia_lct, rr, gg, bb, wavelnth=wave_[4], /load
plot_map, AIA_MAP_CUBE[4], /limb, /log

aia_lct, rr, gg, bb, wavelnth=wave_[5], /load
plot_map, AIA_MAP_CUBE[5], /limb, /log

loadct, 5
plot_map, temperature_map, dmin = 6, /cbar
loadct, 1
plot_map, emission_map, /cbar, dmin = 20



aia_teem_movie, xrange = xrange, yrange = yrange, fileset = 'AIA', npix = 10, /verbose

restore, '/Users/schriste/Dropbox/idl/aia_deem/teem_data_00020110215_000009_npix10.sav', /verbose

restore, '/Users/schriste/Dropbox/idl/aia_deem/teem_tot_00020110215_000009_q94_67.sav'
r = chianti_spec_from_dem(telog, emlog, /plot)

xrange = [274.20, 280.20]
yrange = [-283.80, -277.80]

aia_teem_movie, xrange = xrange, yrange = yrange, fileset = 'AIA', npix = 1

xrange = [130, 190]
yrange = [-283.80, -187.80]

aia_teem_movie, xrange = xrange, yrange = yrange, fileset = 'AIA', npix = 10, /verbose, /force_table

oplot, [xrange[0], xrange[0]], [yrange[0], yrange[1]]
oplot, [xrange[1], xrange[1]], [yrange[0], yrange[1]]
oplot, [xrange[0], xrange[1]], [yrange[0], yrange[0]]
oplot, [xrange[0], xrange[1]], [yrange[1], yrange[1]]