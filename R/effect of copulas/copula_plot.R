rm(list=ls(all=TRUE))

# library(Cairo)
library(copula)
library(ggplot2)
# library(reshape2)
library(extrafont)
library(RColorBrewer)
library(grid)

Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.19/bin/gswin64c.exe")

save_plot = 1
# font_type   = "Times"
font_type   = "CM Roman"

ndim = 2

# ORDER!!! WARNING!
# for legend label
copulas_it     = c("Gauss"           = "Gauss", 
                   "t (dof=2)"       = expression(paste(italic("t"), " (dof=2)")),
                   "Gumbel"          = "Gumbel", 
                   "rotGumbel-180°"  = "rotGumbel-180°",
                   "rotClayton-180°" = "rotClayton-180°")

main_dir    = "d:/Working folder/R working folder/effect of copulas"
fig_dir     = file.path(main_dir, "figures")

n  = 100
x  = seq(-3, 3, length.out = n)
y  = x
X1 = matrix(rep(x, each=n), nrow=n)
X2 = matrix(rep(y, n), nrow=n)

k_tau = 0.5


X12 = matrix(c(c(X1),c(X2)), ncol = 2)

# Elliptical
myMvd1 <- mvdc(copula = ellipCopula(family = "normal", param = sin(k_tau*pi/2)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd2 <- mvdc(copula = ellipCopula(family = "t", df = 2, param = sin(k_tau*pi/2)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

# Archimedean
myMvd3 <- mvdc(copula = archmCopula(family = "clayton", param = 2*k_tau/(1-k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd4 <- mvdc(copula = archmCopula(family = "frank", param = iTau(frankCopula(), k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd5 <- mvdc(copula = archmCopula(family = "joe", param = iTau(joeCopula(), k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd6 <- mvdc(copula = archmCopula(family = "amh", param = iTau(amhCopula(), k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

# Extreme
myMvd7 <- mvdc(copula = archmCopula(family = "gumbel", param = 1/(1-k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd8 <- mvdc(copula = evCopula(family = "tev", param = iTau(tevCopula(df = 2), k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd9 <- mvdc(copula = evCopula(family = "huslerReiss", param = iTau(huslerReissCopula(), k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd10 <- mvdc(copula = evCopula(family = "galambos", param = iTau(galambosCopula(), k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd11 <- mvdc(copula = evCopula(family = "tawn", param = iTau(tawnCopula(), k_tau)),
                margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

# Other
myMvd12 <- mvdc(copula = plackettCopula(param = iTau(plackettCopula(), k_tau)),
               margins = c("norm", "norm"), paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))

myMvd13 <- mvdc(copula = fgmCopula(param = iTau(fgmCopula(), k_tau)),
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

# Elliptical
Z = dMvdc(X12, myMvd1)
df1 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Gauss", class = "elliptical")

Z = dMvdc(X12, myMvd2)
df2 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "t (dof = 2)", class = "elliptical")

# Archimedean
Z = dMvdc(X12, myMvd3)
df3 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Clayton", class = "Archimedean")

Z = dMvdc(X12, myMvd3)
df3r = data.frame(x = rev(X12[,2]), y = rev(X12[,1]), z = Z, type = "rotClayton-180°", class = "Archimedean")

Z = dMvdc(X12, myMvd4)
df4 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Frank", class = "Archimedean")

Z = dMvdc(X12, myMvd5)
df5 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Joe", class = "Archimedean")

Z = dMvdc(X12, myMvd6)
df6 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "AMH", class = "Archimedean")

# Extreme
Z = dMvdc(X12, myMvd7)
df7 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Gumbel", class = "Archimedean & extreme")

Z = dMvdc(X12, myMvd7)
df7r = data.frame(x = rev(X12[,2]), y = rev(X12[,1]), z = Z, type = "rotGumbel-180°", class = "Archimedean & extreme")

Z = dMvdc(X12, myMvd8)
df8 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "extreme t (dof = 2)", class = "extreme")

Z = dMvdc(X12, myMvd9)
df9 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Hüsler-Reiss", class = "extreme")

Z = dMvdc(X12, myMvd10)
df10 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Galambos", class = "extreme")

Z = dMvdc(X12, myMvd11)
df11 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Tawn", class = "extreme")

# Other
Z = dMvdc(X12, myMvd12)
df12 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "Plackett", class = "other")

Z = dMvdc(X12, myMvd13)
df13 = data.frame(x = X12[,1], y = X12[,2], z = Z, type = "FGM", class = "other")


df = rbind(df1, df2, df7, df7r, df3r)
# df = rbind(df1, df2, df3, df4, df5, df6)
# df = rbind(df1, df2, df3, df4, df5, df6, df7, df7r, df8, df9, df10, df11, df12, df13)

### coloring dataframe
## col_df = data.frame(levels(df$type))
## names(col_df) = "type"
## col_df$x = col_df$y = col_df$z = 1
## col_df$class = as.factor(levels(df$class)[1])

g = ggplot(df, aes(x, y, z = z))
g = g + stat_contour(geom="polygon", aes(fill=..level.., group = type), bins = 10)
g = g + stat_contour(aes(group = type), bins = 5, col = 'black')
g = g + scale_fill_gradient(low = "#ffff99", high = "#ff3300")
# g = g + facet_wrap(~ type)
g = g + facet_wrap(c("type", "class"), ncol = 5)
# g = g + facet_grid(type ~ class)
## g = g + geom_rect(data = col_df, aes(fill = type), xmin = -Inf, xmax = Inf,
##                   ymin = -Inf, ymax = Inf, alpha = 0.3)

g = g + theme_bw(base_size = 10, base_family = font_type)

g = g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
g = g + theme(panel.border = element_rect(fill = NULL, color = 'black', size = 0.3))

g = g + theme(legend.key = element_blank())
g = g + theme(legend.position = "none")

# g = g + theme(axis.ticks = element_line(size = 0.3), axis.line = element_line(size = 0.3, colour = "black"))
g = g + theme(axis.ticks = element_line(size = 0.3), axis.line = element_blank())

g = g + theme(strip.background = element_rect(colour = "black", size = 0.3, fill = "white"))
g = g + theme(strip.text.x = element_text(size = 7))
g = g + theme(strip.text = element_text(color = "black"))

# g = g + theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5))

g = g + theme(plot.margin = unit(c(0,1,0,0), "mm"))

g = g + ylab('') + xlab('')
g = g + coord_fixed()
print(g)

## SAVE PLOT
if (save_plot == 1){
  # setEPS()
  fig_name = "5_copulas.pdf" 
  # postscript(file = fig_name, family = font_type,
  #            width = 150/25.4, height = 50/25.4, pointsize = 9,
  #            colormodel = 'cmyk')
  pdf(file = file.path(fig_dir, fig_name), family = font_type,
          width = 160/25.4, height = 50/25.4,
          colormodel = 'cmyk')
  print(g)
  
  dev.off()
  
  embed_fonts(file.path(fig_dir, fig_name))
}