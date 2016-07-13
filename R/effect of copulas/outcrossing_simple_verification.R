rm(list=ls(all=TRUE))

library(CDVine)
library(pracma)
# library(RcppArmadillo)


cop_type= 4
Ps_comp = 0.9
u1      = Ps_comp
u2      = Ps_comp
tau_Fv  = 1/365

nn      = 20
frac    = logspace(-4,-1,nn)

par     = 10
BiCopCDF(u1, u2, cop_type, par)
BiCopCDF(u1, u2, cop_type+10, par)

u1 + u2 - 1 + BiCopCDF(1-u1, 1-u2, cop_type, par)


BiCopPDF(u1, u2, cop_type, par)
BiCopPDF(u1, u2, cop_type+10, par)

BiCopPDF(1-u1, 1-u2, cop_type, par)

# nu_p    = vector(nn, mode = 'double')
# nu_pr   = vector(nn, mode = 'double')
# for (ii in 1:nn){
#   tau_F       = tau_Fv
#   delta_t     = tau_F*frac[ii]
#   k_tau       = 2/pi*asin(exp(-(delta_t/tau_F)^2))
#   par         = 1/(1-k_tau)
# 
#   
#   # BiCopPDF(u1, u2, cop_type, par)
#   # BiCopPDF(u1, u2, cop_type+10, par)
#   
#   # BiCopPDF(1-u1, 1-u2, cop_type, par)
#   
#   Pf_sys     = Ps_comp - BiCopCDF(u1, u2, cop_type, par)
#   Pf_sysr    = Ps_comp - BiCopCDF(u1, u2, cop_type+10, par)
#   
#   nu_p[ii]   = Pf_sys/delta_t
#   nu_pr[ii]  = Pf_sysr/delta_t
# }
# 
# plot(frac, nu_p, type = 'b', log = 'x')
# # plot(frac, nu_pr, type = 'b')

