PRO aia_fluxes_from_rhessi, params

em=params[0]
t=params[1]
sig=params[2]


restore,'teem_table.sav',/verbose


aia_model_flux=flux[t,sig,*] * em



END