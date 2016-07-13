## Tang et al. (2013). Impact of copula selection on geotechnical reliability under incomplete probability information. Computers and Geotechnics
## 
## SWE   - annual max snow water equivalent [mm]
## WS10  - accompanying (same day as annual swe max) horizontal wind speed measured 3 times a day [m/s]

rm(list=ls(all=TRUE))

# library(ggplot2)
# library(VineCopula)
library(copula)
library(FAdist)
library(stats4)

# read data
data = read.csv2('SWE_WS10_Budapest_2545.csv', header = T, sep = ';', dec = '.')


## ::::::::::::::::::::::::::::::::::
## MARGINAL FIT
## ::::::::::::::::::::::::::::::::::

df_marg_param = data.frame()
for (i in 1:2){
  obs      = data[,i]
  
  sd_max   = sd(obs)
  var_max  = sd_max^2
  mean_max = mean(obs)
  
  scale    = sqrt(6)/pi*sd_max
  location = mean_max - scale*-digamma(1)
  
  nLL_gumbel = function(scale, location) {
    -sum(log(dgumbel(obs, scale, location) ))
  }
  
  fit_gumbel = mle(nLL_gumbel, start = list(scale = scale, location = location),
                   method="Nelder-Mead")
  
  df = data.frame(t(fit_gumbel@coef))
  df_marg_param = rbind(df_marg_param, df)
}

row.names(df_marg_param) = c('SWE', 'WS10')
df_marg_param


## ::::::::::::::::::::::::::::::::::
## COPULA FUNCTION FIT
## ::::::::::::::::::::::::::::::::::
# Pearson correlation coefficient
ro = cor(data, method = 'pearson')[2]

mdata = as.matrix(data)
udata = pobs(mdata)
# plot(udata)

# selection of bivariate copula
# cc = BiCopSelect(udata[,'c'], udata[,'fi'], familyset = 1, selectioncrit = 'AIC')

# define copula function
plackett_cop    = plackettCopula()
frank_cop       = frankCopula(dim = 2)
gauss_cop       = normalCopula(dim = 2)
t_cop           = tCopula(dim = 2, df = 2, df.fixed = T)
clayton_cop     = claytonCopula(dim = 2)
gumbel_cop      = gumbelCopula(dim = 2)


# maximum likleihood fit
plackett_fit_ml = fitCopula(plackett_cop, udata, method="ml", start = 1)
frank_fit_ml    = fitCopula(frank_cop, udata, method="ml", start = 2)
gauss_fit_ml    = fitCopula(gauss_cop, udata, method="ml", start = 0)
t_fit_ml        = fitCopula(t_cop, udata, method="ml", start = 0)
clayton_fit_ml  = fitCopula(clayton_cop, udata, method="ml", start = 1)
gumbel_fit_ml   = fitCopula(gumbel_cop, udata, method="ml", start = 2)

## POST-PROCESSING
# parameters
plackett_theta  = plackett_fit_ml@copula@parameters
plackett_AIC    = -2*plackett_fit_ml@loglik + 2*1
plackett_BIC    = -2*plackett_fit_ml@loglik + log(dim(data)[1])*1

frank_theta     = frank_fit_ml@copula@parameters
frank_AIC       = -2*frank_fit_ml@loglik + 2*1
frank_BIC       = -2*frank_fit_ml@loglik + log(dim(data)[1])*1

gauss_theta     = gauss_fit_ml@copula@parameters
gauss_AIC       = -2*gauss_fit_ml@loglik + 2*1
gauss_BIC       = -2*gauss_fit_ml@loglik + log(dim(data)[1])*1

t_theta         = t_fit_ml@copula@parameters
t_AIC           = -2*t_fit_ml@loglik + 2*1
t_BIC           = -2*t_fit_ml@loglik + log(dim(data)[1])*1

clayton_theta   = clayton_fit_ml@copula@parameters
clayton_AIC     = -2*clayton_fit_ml@loglik + 2*1
clayton_BIC     = -2*clayton_fit_ml@loglik + log(dim(data)[1])*1

gumbel_theta    = gumbel_fit_ml@copula@parameters
gumbel_AIC      = -2*gumbel_fit_ml@loglik + 2*1
gumbel_BIC      = -2*gumbel_fit_ml@loglik + log(dim(data)[1])*1

df = data.frame(theta = c(gauss_theta, t_theta, clayton_theta, gumbel_theta, plackett_theta, frank_theta), 
               AIC = c(gauss_AIC, t_AIC, clayton_AIC, gumbel_AIC, plackett_AIC, frank_AIC), 
               BIC = c(gauss_BIC, t_BIC, clayton_BIC, gumbel_BIC, plackett_BIC, frank_BIC))

row.names(df) = c('Gauss', 't', 'Clayton', 'Gumbel','Plackett', 'Frank')
df
df$AIC-min(df$AIC)