## Tang et al. (2013). Impact of copula selection on geotechnical reliability under incomplete probability information. Computers and Geotechnics
##
## In the paper the copulas are fitted based on prescribing the equivalence of Pearson correaltion 
## c   - cohesion [kPa]
## fi  - friction angle [°]

rm(list=ls(all=TRUE))

# library(ggplot2)
library(VineCopula)
library(copula)

# read data
rawdata = read.csv2('c_fi_data_geotechnical_example.csv', header = T, sep = ';', dec = '.', skip = 1)

#CS-ET data
data = rawdata[,7:8]
data = na.omit(data)
names(data) = c('c', 'fi')
# plot(data)
# ggplot(data, aes(x = c, y = fi)) + geom_point()

# Pearson correlation coefficient
ro = cor(data, method = 'pearson')[2]

mdata = as.matrix(data)
udata = pobs(mdata)

# selection of bivariate copula
# cc = BiCopSelect(udata[,'c'], udata[,'fi'], familyset = 1, selectioncrit = 'AIC')

#
plackett_cop    = plackettCopula()
frank_cop       = frankCopula(dim = 2)
gauss_cop       = normalCopula(dim = 2)

plackett_fit_ml = fitCopula(plackett_cop, udata, method="ml", start = 0.5)
frank_fit_ml    = fitCopula(frank_cop, udata, method="ml", start = -2)
gauss_fit_ml    = fitCopula(gauss_cop, udata, method="ml", start = 0)

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

df = data.frame(theta = c(gauss_theta, plackett_theta, frank_theta), 
               AIC = c(gauss_AIC, plackett_AIC, frank_AIC), 
               BIC = c(gauss_BIC, plackett_BIC, frank_BIC))

row.names(df) = c('Gauss', 'Plackett', 'Frank')
df
