# Check the effect of finite difference step on approximation of
# derivative! TO VERIFY MATLAB
#
# Time-variant reliability analysis with continuous stochastic process
# investigation of the effect of copula assumption
#
#
# g = R - S(t)
#
# R    - constant
# S(t) - continuous stochastic process

rm(list=ls(all=TRUE))

library(copula)
library(ggplot2)
library(reshape2)
library(pracma)

#--------------------------------------------------------------------------
# OPTIONS
#--------------------------------------------------------------------------

# COPULA
# copulatype = c('Gaussian', 't', 'Clayton', 'Gumbel')
# copulatype = c('Gaussian', 't', 'Gumbel')
copulatype = 'Gaussian'
# copulatype = 't'
# copulatype = 'Clayton'
# copulatype = 'Gumbel' # numerical problems!!
# copulatype = 'Plackett'
# copulatype = 'AMH'
# copulatype = 'Frank'
# copulatype = 'Husler-Reiss'
copulatype = 'Joe'

# survival probability of a component
Ps_comp     = 1 - 1e-8

# correlation 'length' for S process
tau_Fv      = 10/365

# increments
# frac        = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
frac        = logspace(-2,-1,20)

#--------------------------------------------------------------------------
# ANALYSIS/CALCULATION
#--------------------------------------------------------------------------
mm  = length(copulatype)
nn  = length(frac)

Pf      = matrix(nrow = mm, ncol = nn)
NU_p    = matrix(nrow = mm, ncol = nn)
PF_SYS  = matrix(nrow = mm, ncol = nn)
# loop over copula types
for (jj in 1:mm){ 
  
  copula = copulatype[jj]
  
  # loop over correlation increments
  for (ii in 1:nn) {
    
    tau_F   = tau_Fv
    # time increment
    delta_t = tau_F*frac[ii]
    
    # ktau correlation between 'adjactent' F realizations
    k_tau   = 2/pi*asin(exp(-(delta_t/tau_F)^2)) # to get back the typical correlation function for pearson rho
    
    # Pearson correlation for Gauss and t copulas
    ro_F    = sin(k_tau*pi/2)
    
    #------------------------------------------------------------------
    # DIRECT INTEGRATION - INITIALIZATION
    #------------------------------------------------------------------
    # copula
    if (copula == 'Gaussian'){
      theta = ro_F
      cop = normalCopula(param = theta, dim = 2)
    }else if (copula == 't'){
      theta = ro_F
      cop = tCopula(param = theta, dim = 2, df = 2)
    }else if (copula == 'Clayton'){
      theta = 2*k_tau/(1-k_tau)
      cop = claytonCopula(param = theta, dim = 2)
    }else if (copula == 'Gumbel'){
      theta = 1/(1-k_tau)
      cop = gumbelCopula(param = theta, dim = 2)
    }else if (copula == 'Plackett'){
      cop = plackettCopula()
      theta = iTau(cop, k_tau)
      cop = plackettCopula(param = theta)
    }else if (copula == 'AMH'){
      cop = amhCopula()
      theta = iTau(cop, k_tau)
      cop = amhCopula(param = theta, dim = 2)
    }else if (copula == 'Frank'){
      cop = frankCopula()
      theta = iTau(cop, k_tau)
      cop = frankCopula(param = theta, dim = 2)
    }else if (copula == 'Husler-Reiss'){
      cop = huslerReissCopula()
      theta = iTau(cop, k_tau)
      cop = huslerReissCopula(param = theta)
    }
    else if (copula == 'Joe'){
      cop = joeCopula()
      theta = iTau(cop, k_tau)
      cop = joeCopula(param = theta)
    }
    
    #------------------------------------------------------------------
    # DIRECT INTEGRATION - CALCULATION - BUILT-IN COPULACDF
    #------------------------------------------------------------------
    if (1-pCopula(c(Ps_comp, Ps_comp), cop) == 0){
      warning('trouble with copula function!')
    }
    Pf_sys        = Ps_comp - pCopula(c(Ps_comp, Ps_comp), cop)
    
    PF_SYS[jj,ii] = Pf_sys
    
    nu_p          = Pf_sys/delta_t
    NU_p[jj,ii]   = nu_p
    
    Pf0           = 1 - Ps_comp
    Pf[jj,ii]     = Pf0 + nu_p*50
  }
}

# plot(frac, PF_SYS, type = 'b')
# R = apply(Pf, 2, function(x) x/Pf[,1])
R = Pf
if (length(R) == nn){
  plot(frac, R, type = 'b', log = 'x')
} else {
  for (ii in 1:mm){
    if(ii == 1){
      plot(frac, R[ii,], type = 'b', log = 'x', ylim = c(0,1.2))
    } else {
      lines(frac, R[ii,], type = 'b')
    }
  }
}



# plot(frac, Pf/Pf[1], type = 'b')

# semilogx(frac, PF_SYS, '-o')
# loglog(frac, PF_SYS, '-o')
# legend(copulatype, 'Location', 'best')
# ylabel('P_{f,sys}(frac)')
# xlabel('frac = \Delta t/corr length')
# grid on


# figure
# nPf = bsxfun(@rdivide, Pf, Pf(:,1))
# semilogx(frac, nPf, '-o')
# legend(copulatype, 'Location', 'best')
# ylabel('P_f(frac)/P_f(min(frac))')
# xlabel('frac = \Delta t/corr length')
# grid on