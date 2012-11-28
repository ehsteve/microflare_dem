FUNCTION get_energy_from_dem,temperatures,params,area,n=n,epstein=epstein

;params=[a0,a1,a2]

dem=get_dem_from_params(temperatures,params,n=n,epstein=epstein)

temp_lin=10^temperatures
area=double(area)

;get the radiative loss function
rad_loss,rad_temps,loss,density=1e11

;interpolate the radiative loss function before performing integral
rlossfunction=interpol(rad_loss,rad_temps,temp_lin,/quadratic)

total_energy_per_sec=int_tabulated(temp_lin,(rlossfunction*dem*area))

print,'total radiated energy per second: ' + num2str(total_energy_per_sec) + ' (erg/s)'

return,total_energy_per_sec

END
