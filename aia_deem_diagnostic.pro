PRO aia_deem_diagnostic,te_min=te_min,te_max=te_max,empeak_log=empeak_log,gauss_params1=gauss_params1,gauss_params2=gauss_params2,$
type=type,q94=q94,diagnostic_struct=diagnostic_struct,ps=ps,quiet=quiet

; Project     : AIA/SDO
;
; Name        : AIA_DEEM_DIAGNOSTIC 
;
; Category    : Data analysis/simulation   
;
; Explanation : compares a known input DEM function with the result of a best-fit Guassian DEM, by searching over peak temperature and width
;		c.f. Markus' aia_teem_map.pro and aia_teem_total.pro. The user inputs parameters in order to generate a desired DEM function. 
;		'Real' flux values are then generated for each AIA wavelength based on this DEM, using the AIA response function. 
;		Then, a best-fit single Gaussian model is found which best reproduces those flux values. The current options are
;		a single gaussian profile, or a double gaussian profile, of variable peak T, width, and amplitude. 
;
; Syntax      : IDL> aia_deem_diagnostic
;
; Inputs      : No mandatory inputs
;
; Outputs     : No mandatory outputs
;
; Keywords    : TE_MIN = minimum (log) temperature value to consider, e.g. 5.5
;		TE_MAX = maximum (log) temperature value to consider, e.g. 8.0
;		GAUSS_PARAMS1 = a vector containing the Amplitude, peak temperature (in log space) and width (log space) of the first Gaussian
;		in the user-generated DEM function. For example, GAUSS_PARAMS1=[1.0,6.2,0.3] for a Gaussian peaking at log T = 6.2 with width 0.3
;		GAUSS_PARAMS2 = a vector containing the Amplitude, peak temperature (in log space) and width (log space) of the second Gaussian
;		in the user-generated DEM function. For example, GAUSS_PARAMS2=[1.0,7.0,0.2] for a Gaussian peaking at log T = 7.0 with width 0.2
;		TYPE = the type of function to generate. Current options are:	'gaussian'
;										'double_gaussian'
;		EMPEAK_LOG = baseline emission measure value to use for the user-generated DEM, e.g. 22
; 		Q94 = the desired empirical boost factor to use on the 94A channel at log(T) < 6.3. Default is currently 1
; 		DIAGNOSTIC_STRUCT = a structure containing all of the input and output parameters for reference
;		QUIET = If set, plots are not sent to the window, and most print info is omitted
;		
; 
; Discontinued: TPEAK1 = peak position in (log) temperature of the first Gaussian in the user-generated DEM function
; Keywords	TPEAK2 = peak position in (log) temperature of the second Gaussian in the user-generated DEM function, if using a double-gaussian model
;		SIG1 = the width of the first Gaussian in the user-generated DEM function, in units of log T
;		SIG2 = the width of the second Gaussian in the user-generated DEM function, in units of log T, if using a double-gaussian model
;		AMP1 = amplitude modifier of the first Gaussian in the user-generated DEM function. Default is 1.
;		AMP2 = amplitude modifier of the second Gaussian in the user-generated DEM function, if using a double-gaussian model
;		FORCE_TABLE = generate aia_teem_table.sav or overwrite existing aia_teem_table.sav. If set, then DIR is also needed
; 		DIR = the directory used to search for AIA files to generate an ordered list. Only needed if FORCE_TABLE is set
;
;
; History     :  9-Apr-2012, Version 1 written by Andrew R. Inglis
;		 11-Apr-2012, Version 1.1, ARI:	moved the teem_table lookup in-house, and removed the call to aia_teem_table. Force_table keyword removed.
;					        altered input keywords. Switched to using GAUSS_PARAMS1=[Amp1,Tpeak1,sig1],
;						GAUSS_PARAMS2=[Amp2,tpeak2,sig2]	
;					       	discontinued FORCE_TABLE,DIR,AMP1,AMP2,TPEAK1,TPEAK2,SIG1,SIG2
;						added DIAGNOSTIC STRUCT,QUIET keywords	
;
; Example     : aia_teem_diagnostic,type='double_gaussian',gauss_params1=[1.0,6.2,0.3],gauss_params2=[0.5,7.0,0.2],empeak_log=22.,diagnostic_struct=diagnostic_struct
;
;				 
; Issues      : /ps keyword currently broken
;
; Contact     : andrew.inglis@nasa.gov
;-

;--------INITIALISATION----------

;Set up some default parameters so that the user doesn't have to specify everything/anything. 

default,t_min,5.5
default,t_max,8.0
default,gauss_params1,[1.0,6.2,0.3]
default,gauss_params2,[1.0,7.0,0.3]
;default,tpeak1,6.2
;default,tpeak2,7.1
default,empeak_log,22.0
;default,amp1,1.0
;default,amp2,1.0
;default,sig1,0.3
;default,sig2,0.3
default,type,'gaussian'
default,q94,1.0
default,dir,'/home/ainglis/physics/steven_code/aschwanden_test/'

texp_=1.


;grab the input parameters and rename them for the rest of the code
amp1=gauss_params1[0]
tpeak1=gauss_params1[1]
sig1=gauss_params1[2]
amp2=gauss_params2[0]
tpeak2=gauss_params2[1]
sig2=gauss_params2[2]


wave_ =['131','171','193','211','335','94']
nwave =n_elements(wave_)

;generate the temperature array based on the given t_min and t_max. Currently, the binning of temperature is hardcoded at t_d=0.05
t_d = 0.05
telog = t_d * findgen((t_max - t_min)/t_d) + t_min
nte = n_elements(telog)

tsig_min = 0.01
tsig_max = 1.0
tsig_d = 0.01
tsig = tsig_d * findgen((tsig_max - tsig_min)/tsig_d) + tsig_min
nsig=n_elements(tsig)
nfree = 3

;convert the emission measure baseline value from a log to a linear value
empeak=(10^empeak_log)



;--------GET RESPONSE FUNCTION-------

;here we obtain the AIA response functions for each wavelength and subsequently use them to generate 'real' flux values for each wavelength


;get the AIA response function with desired keywords
tresp = aia_get_response(/temp,/full,/chiantifix,/evenorm,/phot)
telog_  =tresp.logte
; put the response in the right order
FOR iw = 0, nwave-1 DO BEGIN
	IF iw EQ 0 THEN BEGIN
		nte_ = n_elements(tresp.tresp[*,0])
		nwave_ = n_elements(tresp.tresp[0,*])
		resp_ = fltarr(nte_, nwave)
		logte_tmp = fltarr(nte_,nwave)
	ENDIF
	filter ='A'+wave_(iw)
	IF (wave_[iw] eq '094') THEN filter='A94'
	ichan = where(tresp.channels EQ filter)
	resp_[*,iw] = tresp.tresp[*,ichan]
	logte_tmp[*,iw]=tresp.logte[*,ichan]
ENDFOR

	nte = n_elements(telog)
	resp    =fltarr(nte,nwave)
	FOR iw = 0, nwave-1 DO resp[*,iw] = interpol(resp_[*,iw], telog_, telog, /quadratic)

;EMPIRICAL CORRECTION 94 A
ind1 = where(telog le 6.3)
resp_corr = resp
resp_corr[ind1,5] = resp[ind1,5] * q94


;plot the response functions for each wavelength to make sure they are in the right order
IF NOT keyword_set(quiet) THEN BEGIN
!p.multi=[0,3,2]
window,2
for i=0,nwave-1 do begin
	plot,telog,alog10(resp[*,i]),charsize=2.5,xran=[4,8],yran=[-28,-24],/ynozero
	al_legend,[wave_(i)],/right
endfor
ENDIF
;help,logte_tmp
;help,resp_
;print,resp_(*,1)

!p.multi=[0,1,1]

;------CREATE MODEL DEM CURVE------------
;Generate the user's DEM curve. Currently there are two options, a single or a double-gaussian function. The single gaussian is mainly for testing purposes.

IF (type eq 'gaussian') THEN BEGIN
	dem_curve=amp1*empeak*exp(-(telog-tpeak1)^2/(2.*sig1^2))
ENDIF

IF (type eq 'double_gaussian') THEN BEGIN
	dem_curve=(amp1*empeak*exp(-(telog-tpeak1)^2/(2.*sig1^2))) + (amp2*empeak*exp(-(telog-tpeak2)^2/(2.*sig2^2)))
ENDIF


;-----CALCULATE FLUX IN EACH WAVELENGTH----------

;get the flux values for each wavelength
flux_obs=fltarr(nwave)

dte1	=10.^telog(1:nte-1)-10.^telog(0:nte-2)
dte	=[dte1(0),dte1]

for i=0, nwave-1 do begin
	flux_obs(i)=total(dem_curve*resp(*,i)*dte)
endfor


;if teem_table.sav does not exist or needs to be overwritten then FORCE_TABLE is used here. Note that DIR is also needed here otherwise file_list is undefined.
;IF keyword_set(FORCE_TABLE) THEN aia_teem_table,wave_,tsig,telog=telog,q94,filelist=file_list[0,*], 'teem_table.sav'


;---------------CALCULATES LOOKUP TABLES----------------------
dte1	=10.^telog(1:nte-1)-10.^telog(0:nte-2)
dte	=[dte1(0),dte1]
em1	=1.
nsig	=n_elements(tsig)
flux	=fltarr(nte,nsig,nwave)
for i=0,nte-1 do begin
 for j=0,nsig-1 do begin
  em_te =em1*exp(-(telog-telog(i))^2/(2.*tsig(j)^2))
  for iw=0,nwave-1 do flux(i,j,iw)=total(resp_corr(*,iw)*em_te*dte)
 endfor
endfor
;save,filename=save_dir + teem_table,wave_,q94,area,resp_corr,telog,dte,tsig,flux
;print,'Lookup table created : ',teem_table
;end


;--------NOW USE THE CODE TO FIND THE BEST FIT---------

;Now search X^2 space using aia_teem_table.sav to find the values of T and sig that best fit the 'real' flux values.
;restore,'teem_table.sav',/verbose

		flux_obs=flux_obs
		counts=flux_obs*texp_
		noise=sqrt(counts)/texp_
		chimin = 9999.
		chi6min = 9999.
		flux_dem_best = 9999.
		chi2d = fltarr(nte,nsig)
		FOR k=0, nte-1 DO BEGIN
   			FOR l = 0, nsig-1 DO BEGIN
				flux_dem1=reform(flux(k,l,*))
				em1	=total(flux_obs)/total(flux_dem1)
				em1_err = sqrt(total(noise^2))/total(flux_dem1)
				flux_dem=flux_dem1*em1
				chi	=sqrt(total((flux_obs-flux_dem)^2/noise^2)/(nwave-nfree))
				chi6=abs(flux_obs-flux_dem)/noise
				chi2d[k,l] = chi				
				IF (chi le chimin) THEN BEGIN
					chimin	=chi
					chi6min	=chi6
					em_best	=alog10(em1)
					em_best_err = em1_err/em_best
					te_best	=telog(k)
					sig_best	=tsig(l)
					flux_dem_best = flux_dem
				ENDIF
			ENDFOR 
		ENDFOR

;---------------CALCULATE THE ERROR IN T and SIG-------------------------

IF (chimin NE 9999.) THEN BEGIN
		;chimin defaults to 9999. if bad pixel (i.e. no data), so skip those
		;find the 1 sigma contour in chi^2-space to get errors in tsig and telog
		contour, chi2d, telog, tsig, levels = [chimin +2.3],path_xy=path_xy,path_info=path_info,/path_data_coords,closed=0
		IF exist(path_xy) THEN BEGIN
			testxmax=MAX(path_xy(0,*))
			testxmin=MIN(path_xy(0,*))
			testymax=MAX(path_xy(1,*))
			testymin=MIN(path_xy(1,*))
		
			telog_err=[testxmin,testxmax] - te_best
			;error bars are asymmetric - symmetrise!
			telog_errsymmetric=max(abs(telog_err));sqrt(telog_err(0)^2 + telog_err(1)^2)
		
			tsig_err=[sig_best - testymin,testymax - sig_best]
			;error bars are asymmetric - symmetrise!
			;tsig_errsymmetric=sqrt(tsig_err(0)^2 + tsig_err(1)^2)
			tsig_errsymmetric = max(abs(tsig_err))
			print,telog_errsymmetric
			print,tsig_errsymmetric
			print,telog_err
			print,tsig_err
		ENDIF ELSE BEGIN
			;contingency statements in case bad pixel
			telog_err=[0,0]
			tsig_err=[0,0]
			tsig_errsymmetric=0.
			telog_errsymmetric=0.
		ENDELSE
	ENDIF ELSE BEGIN
		;contingency statements in case bad pixel
		telog_err=[0,0]
		tsig_err=[0,0]
		tsig_errsymmetric=0.
		telog_errsymmetric=0.
	ENDELSE 
		;plot the X^2 contour
		IF NOT keyword_set(quiet) THEN BEGIN
			IF chimin LT 9999. THEN BEGIN
			window,7
			loadct, 0
			hsi_linecolors
			nlevels = 20
			levels = chimin * (findgen(nlevels)*0.1 + 1.0)
			anot = strarr(20)
			FOR k = 0, nlevels-1 DO anot[k] = num2str(levels[k])
			contour, chi2d, telog, tsig, levels = levels, c_annotation = anot, xtitle = 'log[Temp]', ytitle = 't_sig'
			oplot, [te_best], [sig_best], psym = symcat(16), color = 6
			leg = num2str(chimin) + '[' + num2str(te_best) + ',' + num2str(sig_best) + ']'
			legend, leg, psym = 4
			contour,chi2d,telog,tsig,levels=[chimin + 2.3],/over,thick=2, color = 6
			ENDIF
		ENDIF
;-------------MAKE COMPUTED DEM CURVE BASED ON BEST-FIT PARAMETERS--------------
qflux=flux_dem_best/flux_obs
uncert=(noise/flux_obs)*(flux_dem_best/flux_obs)

computed_curve=10^(em_best)*exp(-(telog-te_best)^2/(2.*sig_best^2))

		;now work out the uncertainty on each data point in the best-fit DEM curve. Need to check this.

		;error1=double(tsig_errsymmetric)
  		;error2=double(telog_errsymmetric)
		;error3=(2.*(error2/(telog-te_best))*((telog-te_best)^2.)) ; error in (telog-te)^2 - careful with minus sign.
		;error4=(2.*(error1/sig_best))*(sig_best^2.) ; error in sig^2
		;error5=2.* error4 ; error in 2sig^2
		;error6=sqrt( (error3/(telog-te_best)^2.) + (error5/(2.*sig_best^2.)) ) * ((telog-te_best)^2. / (2.*sig_best^2.)) ;error in exponent: {-(telog-te)^2/(2*sig^2)} 
		;error7=error6 * exp(-(telog-te_best)^2./(2.*sig_best^2.)) ;error once the exp is done.
		;error8=error7*(10^(em_best));/(nmacro*mask_frac)) ; final error in em at i,j 01/27/12 - added mask_frac - ARI
		;check FOR -NaN values
		;p=finite(error8)
		;index = WHERE (p gt 0,count,complement=index_c)
		;replace NaN values with 0.
		;error8(index_c)=0.
		;dem_error=error8

		expression1=telog-te_best
		expression2=expression1^2
		expression3=sig_best^2
		expression4=2*expression3
		expression5=expression2/expression4
		expression6=exp(-expression5)
		expression7=10^(em_best)*expression6

		error1=double(tsig_errsymmetric)
		error2=double(telog_errsymmetric)
		error3=(2*error2/expression1) * expression2
		error4=(2*error1/sig_best) * expression3
		error5=error4*2
		error6=sqrt( (error3/expression2)^2 + (error5/expression4)^2) * expression5
		error7=error6*expression7
		error8=error7
		dem_error=error8

		;print,'error4:',error4
		;print,'error5:',error5
		;print,'error6:',error6
		

		lower_bound_curve=computed_curve-dem_error
		upper_bound_curve=computed_curve+dem_error

		;when error is larger than DEM curve value lower bound can be lt 0. Unphysical. For plotting purposes set to 1.
		index = WHERE (lower_bound_curve lt 0,count,complement=index_c)
		lower_bound_curve(index)=1.
;stop
;------DISPLAY RESULTS------------

;plot the user-generated DEM function with an overlay of the best-fit single Gaussian result.
IF NOT keyword_set(quiet) THEN BEGIN

	IF keyword_set(ps) THEN BEGIN
		set_plot,'ps'
		device,encaps=1,decomposed=0,bits_per_pixel=8,filename='aia_deem_diagnostics.ps'
	ENDIF ELSE BEGIN 
		window,5
	ENDELSE

	loadct,39
	plot,(telog),alog10(dem_curve),thick=2,charsize=1.5,xtitle='log T (K)',ytitle='emission measure log(EM [cm!U-5!N K!U-1!N])',ystyle=1,yrange=[max(alog10(dem_curve)) - 4.,max(alog10(dem_curve))+1.]
	oplot,(telog),alog10(computed_curve),thick=2,linestyle=2,color=220
	;oploterr, telog, (computed_curve), (dem_error),color=220,errcolor=220
	oplot,(telog),alog10(lower_bound_curve),thick=2,linestyle=0,color=220
	oplot,(telog),alog10(upper_bound_curve),thick=2,linestyle=0,color=220

	;overlay the upper and lower bounds of the DEM curve, grabbed from the error calculations
	errplot,(telog),alog10(lower_bound_curve),alog10(upper_bound_curve),width=0,color=220,$
	clip=[telog(0),max(alog10(dem_curve)) - 4.,telog(nte-1),max(alog10(dem_curve))+1.0],noclip=0
	al_legend,['T = ' +num2str(te_best),'sig = ' +num2str(sig_best),'em = ' +num2str(em_best),'chimin = ' +num2str(chimin)],/right,/top,charsize=1.5
	;pattern=[[255,0,0,0,0,0,0,0],[255,0,0,0,0,0,0,0]]
	;polyfill,[telog,reverse(telog)],[alog10(upper_bound_curve),alog10(reverse(lower_bound_curve))],color=220, $
	;pattern=pattern,clip=[telog(0),max(alog10(dem_curve)) - 4.,telog(nte-1),max(alog10(dem_curve))+1.],noclip=0

	;some information statements
	print,flux_obs
	print,flux_dem_best
	print,flux_dem_best/flux_obs
	print,'T = ',te_best
	print,'sig = ',sig_best
	print,'chimin = ',chimin
	print,'dem error = ',dem_error
	;print,10^(dem_error)



;plot the ratio of the user generated flux and the best-fit flux in each wavelength

	IF NOT keyword_set(ps) THEN BEGIN
	window,4,xsize=756
	ENDIF ELSE BEGIN
	ERASE
	ENDELSE

	loadct,0
	;xtickname=['131A','171A','193A','211A','335A','94A']
	;qflux=flux_dem_best/flux_obs
	plot,findgen(6),qflux,psym=4,charsize=1.2,/ynozero,symsize=2,thick=2,ytitle='sim. flux/obs.flux', $
	xrange=[-0.5,5.5],/xstyle,xtickname=['131A','171A','193A','211A','335A','94A']
	oplot,findgen(8)-1,findgen(8)*0. +1.,linestyle=2

	;overplot the uncertainty in the flux ratio - comes from 'noise'

	oploterr,findgen(6),qflux,uncert,psym=4

	al_legend,['T = ' +num2str(te_best),'sig = ' +num2str(sig_best),'em = ' +num2str(em_best),'chimin = ' +num2str(chimin)],/left,/top,charsize=1.5
	al_legend,['131A : ' +num2str(qflux(0)),'171A : ' +num2str(qflux(1)),'193A : ' +num2str(qflux(2)),'211A : ' $ 
	+num2str(qflux(3)),'335A : ' +num2str(qflux(4)),'94A : ' +num2str(qflux(5))],/right,/top,charsize=1.2

	IF keyword_set(ps) THEN BEGIN
	device,/close
	set_plot,'x'
	ENDIF

ENDIF

;create a structure containing all of the inputs, outputs and dem curves generated in this routine
diagnostic_struct=create_struct('dem_curve',dem_curve,'computed_curve',computed_curve,'upper_bound_curve',upper_bound_curve,$
'lower_bound_curve',lower_bound_curve,'gauss_params1',gauss_params1,'gauss_params2',gauss_params2,'empeak_log',empeak_log,'amp_ratio',amp1/amp2,$
'T',te_best,'sig',sig_best,'chimin',chimin,'wave_',wave_,'qflux',qflux,'qflux_error',uncert,'type',type)


END