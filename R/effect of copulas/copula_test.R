rm(list=ls(all=TRUE))

library(copula)
library(ggplot2)
library(reshape2)
# library(gumbel)

ndim = 2

# cop <- new("acopula")
# cop@name <- "Gumbel"
# cop@theta <- 0.5

# cop <- archmCopula(family = "Gumbel", dim = ndim, param = 2)

# cop <- normalCopula(0.5, dim = ndim)
# cop <- tCopula(c(0.5, 0.3), dim = ndim, dispstr = "toep")

n = 100
x = seq(-3, 3, length.out = n)
y = x
X1 = matrix(rep(x, each=n), nrow=n)
X2 = matrix(rep(y, n), nrow=n)

k_tau = 0.5

# cop <- evCopula(family = "tev", param = c(0.8, 0.9), df = 2, dim = ndim, dispstr = "toep")
# cop <- evCopula(family = "huslerReiss", param = iTau(huslerReissCopula(), k_tau), dim = ndim, dispstr = "toep")
# # cop <- evCopula(family = "huslerReiss", param = 1e3, dim = ndim, dispstr = "toep")
# 
# 
# cop <- ellipCopula(family = "normal", dim = ndim, dispstr = "toep", param = sin(k_tau*pi/2))
# 
# # pCopula(rep(0.8, ndim), cop)
# 
X12 = matrix(c(c(X1),c(X2)), ncol = 2)
# pdf = matrix(dCopula(U12, cop), ncol = n)
# # 
# contour(z = pdf)
# # 


# xx = seq(0,1,length.out = 1e5)
# yy = iTau(huslerReissCopula(), xx)
# plot(xx, yy, type = "b")


# writeMat("hr_tau_delta.mat", hr_tau = xx, hr_delta = yy)
# u = matrix(rep(0.5, ndim), ncol = ndim)
# dgumbel(u, alpha = 2, dim = ndim)

myMvd1 <- mvdc(copula = ellipCopula(family = "normal", param = sin(k_tau*pi/2)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd2 <- mvdc(copula = ellipCopula(family = "t", df = 2, param = sin(k_tau*pi/2)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd3 <- mvdc(copula = archmCopula(family = "clayton", param = 2*k_tau/(1-k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd4 <- mvdc(copula = archmCopula(family = "gumbel", param = 1/(1-k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd5 <- mvdc(copula = evCopula(family = "tev", param = iTau(tevCopula(df = 2), k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd6 <- mvdc(copula = evCopula(family = "huslerReiss", param = iTau(huslerReissCopula(), k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

# par(mfrow = c(2, 3), mar = c(2, 2, 1, 1), oma = c(1, 1, 0, 0),
#     mgp = c(2, 1, 0))
# par(mfrow = c(2, 2), mar = c(2, 2, 1, 1), oma = c(1, 1, 0, 0),
#     mgp = c(2, 1, 0))
# 
# contour(myMvd1, dMvdc, xlim = c(-3, 3), ylim = c(-3, 3))
# text(-2,2, 'Gauss')
# contour(myMvd2, dMvdc, xlim = c(-3, 3), ylim = c(-3, 3))
# text(-2,2, 't (dof = 2)')
# contour(myMvd3, dMvdc, xlim = c(-3, 3), ylim = c(-3, 3))
# text(-2,2, 'Clayton')
# contour(myMvd4, dMvdc, xlim = c(-3, 3), ylim = c(-3, 3))
# text(-2,2, 'Gumbel')
# contour(myMvd5, dMvdc, xlim = c(-3, 3), ylim = c(-3, 3))
# text(-2,2, 'extreme t (dof = 2)')
# contour(myMvd6, dMvdc, xlim = c(-3, 3), ylim = c(-3, 3))
# text(-2,2, 'Hüsler-Reiss')

Z = dMvdc(X12, myMvd1)
df1 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Gauss")

Z = dMvdc(X12, myMvd2)
df2 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "t (dof = 2)")

Z = dMvdc(X12, myMvd3)
df3 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Clayton")

Z = dMvdc(X12, myMvd4)
df4 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Gumbel")

df = rbind(df1, df2, df3, df4)

g = ggplot(df, aes(x, y, z = z))
g = g + stat_contour(geom="polygon", aes(fill=..level.., group = type))
g = g + facet_wrap(~ type)

g = g + theme_bw(base_size = 10, base_family = 'Times New Roman') 
g = g + theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5))
# g = g + theme(legend.position="none")
g = g + theme(plot.margin = unit(c(0,0,0,0), "mm"))
g = g + ylab('') + xlab('')
print(g)


cop <- huslerReissCopula()
cop <- tevCopula(df=2)

# tau = c(seq(-1,0.90,length.out = 1e2), exp(seq(log(0.9), log(1), length.out = 1e3)))
tau = c(seq(0,1,length.out = 1e3))
tau = unique(tau)
delta = iTau(cop, tau)
plot(tau, delta)

# writeMat('d:/Working folder/Matlab working folder/add_function/statistics/distribution_functions/multivariate/husler_reiss/hr_tau_delta.mat', ktau = tau, delta = delta)
writeMat('d:/Working folder/Matlab working folder/add_function/statistics/distribution_functions/multivariate/t_extreme/tev_tau_rho.mat', ktau = tau, rho = delta)


cop <- claytonCopula(param = 10)
