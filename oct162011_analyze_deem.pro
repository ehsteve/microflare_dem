PRO oct162011_analyze_deem

workdir ='~/idlpro/schriste/aiatest/' 

wave_ =['131','171','193','211','335','94'] 
nwave =n_elements(wave_) 

files = file_search(workdir + 'aia_test*')

fileset = 'ssw_cutout'

te_range=[0.5,20]*1.e6 ;   ([K], valid temperature range for DEM solutions) 
tsig=0.1*(1+1*findgen(10)) ;   (values of Gaussian logarithmic temperature widths) 
q94=6.7 ;   (correction factor for low-temperature 94 A response) 

fov=[0,0,4095,4095] ;   (pixel range [i1,j1,i2,j2] of full image) 

npix=4 ;   (macropixel size=4x4 pixels, yields 512x512 map) 

;not sure what the following does
vers='a' ;   (version number of label in filenames used for full images) 

; Savefile that contains DEM loopup table
teem_table='teem_table.sav' 
; Savefile that contains EM and Te maps
teem_map =fileset+vers+'_teem_map.sav' 
;  jpg-file that shows EM and Te maps
teem_jpeg=fileset+vers+'_teem_map.jpg' 

aia_teem_table,wave_,tsig,te_range,q94,fileset=fileset, teem_table
for u=0,150 do begin
aia_teem_map,wave_,npix,teem_table,teem_map,fileset='ssw_cutout',u
print,u

;aia_teem_disp,teem_map,te_range,t1,teem_jpeg 

restore, teem_map, /verbose

restore, 'ssw_cutouta_teem_map.sav', /verbose

hist = histogram(te_map, binsize = 0.01, min = 5.5, max = 7, locations = loc)
plot, loc, hist, psym = 10, xtitle = 'log(Temperature [K])', title = 'Histogram'


hist = histogram(em_map, binsize = 0.01, min = 18, max = 22, locations = loc)
plot, loc, hist, psym = 10, xtitle = 'log(Emission Measure [cm!U-3!N])', title = 'Histogram'

!P.MULTI = [0,2,1]

plot_map, temperature_map, dmin = 5.5, dmax = 7, /cbar, /limb
plot_map, emission_map, dmin = 18, dmax = 22, /cbar, /limb

xrange = [350, 700]
yrange = [-500, -300]

sub_map, temperature_map, stemmap, xrange = xrange, yrange = yrange
sub_map, emission_map, semmap, xrange = xrange, yrange = yrange
sub_map, aia_map_cube, saia_map_cube, xrange = xrange, yrange = yrange

plot_map, stemmap, /cbar, /limb, grid = 10
plot_map, semmap, /cbar, /limb, grid = 10
plot_map, saia_map_cube[0], /limb, grid = 10, /log

xrange = [475, 750]
yrange = [-450, -300]

sub_map, temperature_map, stemmap, xrange = xrange, yrange = yrange
sub_map, emission_map, semmap, xrange = xrange, yrange = yrange
sub_map, aia_map_cube, saia_map_cube, xrange = xrange, yrange = yrange

plot_map, stemmap, /cbar, /limb, grid = 10
plot_map, semmap, /cbar, /limb, grid = 10
plot_map, saia_map_cube[0], /limb, grid = 10, /log

roi = findgen(n_elements(stemmap))

; savefile that contains total DEM distribution
teem_tot=fileset+vers+'_teem_tot.sav' 

mask_map = oct162011_ribbon_mask()
help,mask_map
;imap=inter_map(mask_map,temperature_map)
;stop
aia_teem_total,'ssw_cutout',npix,wave_,q94,teem_table,teem_map,teem_tot,mask=mask_map,u
endfor
stop

END