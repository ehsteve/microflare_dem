.r ~/idlpro/schriste/f_multi_therm

e1 = findgen(300)/10.0 + 5

edge_products, e1, edges_2=e2, gmean= em

rel_abun = [26,28]

low_t = 0.02
high_t = 2
dt = 0.01

hsi_linecolors

emission_measure = 1

alpha = 10.0
params = [emission_measure, low_t, high_t, alpha]

flux1 = f_multi_therm(e2, params, /pow, /cont)
flux2 = f_multi_therm(e2, params, anydem = [transpose(telog), transpose(emlog)], /cont)

plot, em, flux1, psy=10, xrange=[3,100.], ytitle = ytitle, charsize=1.5, /xstyle, ytickf = "exp1", xtickf = "exp1", /ylog, /xlog, /nodata, xtitle = 'Energy [keV]'

oplot, em, flux1, color = 2, psym = 10
oplot, em, flux2, color = 3, psym = 10

