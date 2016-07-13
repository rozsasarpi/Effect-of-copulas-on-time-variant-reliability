# The effect of copulas on time-variant reliability 
# with continuous stochastic processes
#
# Sudret's corroding beam example
#

rm(list=ls(all=TRUE))

library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
# library(dplyr)
library(extrafont)
library(RColorBrewer)
library(grid)
library(XLConnect)
library(gridExtra)
library(latex2exp)

source("log10_minor_break.R")
source("log10space.R")

Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.19/bin/gswin64c.exe")

copulas     = c("Gauss-FORM", "Gauss-DI", "t (dof=2)-DI", "Gumbel-DI", "rotGumbel-180°-DI", "rotClayton-180°-DI")
copulas_excel = c("Gauss-FORM", "Gauss-DI", "t (dof=2)-DI", "Gumbel-DI", "rotClayton-180°-DI", "rotGumbel-180°-DI")

# for legend label
copulas_it     = c("Gauss-FORM", "Gauss-DI", expression(paste(italic("t"), " (dof=2)-DI")), "Gumbel-DI", "rotGumbel-180°-DI", "rotClayton-180°-DI")

nc          = length(copulas)
# set colors to copula types
cbSet1      = brewer.pal(nc, "Set1")
cbPalette   = c(cbSet1[1], cbSet1)
# set line types
ltPalette   = c(2,rep(1,nc))
shPalette   = c(1+15, 1+15, 2+15, 0+15, 3, 5+13) 

save_plot   = 1
log_plot    = 1
# font_type   = "Times"
font_type   = "CM Roman"

main_dir    = "d:/Working folder/R working folder/effect of copulas"
data_dir    = "d:/Working folder/Matlab working folder/time variant reliability/phi2 method/effect of copulas/Sudret's corroding beam"
fig_dir     = file.path(main_dir, "figures")

excel_file  = file.path(data_dir, "corroding beam_results2.xlsx")
sheet_name  = paste("results")

################################################################################
## LOAD DATA
################################################################################
wb          = loadWorkbook(excel_file)
nu_df       = readWorksheet(wb, sheet = sheet_name, region = "B4:G13", useCachedValues = TRUE)
# beta_df     = readWorksheet(wb, sheet = sheet_name, region = "B95:G101", useCachedValues = TRUE)

# WARNING!
beta_df     = nu_df

################################################################################
## PRE-PROCESS DATA
################################################################################
names(beta_df) = signif(beta_df[1,],2)
beta_df = beta_df[-c(1,2,6),]

tmp_df = beta_df
for (ii in 1:(dim(beta_df)[2])){
  tmp_df[ii,] = tmp_df[ii,]/beta_df[2,]
}

beta_df$copula = as.factor(copulas_excel)
beta_df$id     = as.factor("raw")
tmp_df$copula  = as.factor(copulas_excel)
tmp_df$id      = as.factor("normalized")
beta_df        = rbind(beta_df, tmp_df)

beta_df = melt(beta_df)
names(beta_df) = c("copula", "id", "t", "beta")
beta_df$t = as.numeric(levels(beta_df$t))[beta_df$t]

# convert ID-s to factors and change the order of specimens for plot
beta_df = within(beta_df, copula <- factor(copula, levels = as.factor(copulas)))

# ltID
beta_df$ltID = 1
beta_df$ltID[beta_df$copula == copulas[1]] = 2
beta_df$ltID = as.factor(beta_df$ltID)
################################################################################
## PLOT
################################################################################
basic_plot = function(df) {
  p = ggplot(df, aes(x = t, y = beta, color = copula, shape = copula))
  p = p + geom_path(aes(linetype = copula))
  p = p + geom_point(color = "white", size = 2)
  p = p + geom_point(size = 1.0)
  # p = p + facet_grid(id~., scales = "free_y")
  p = p + theme_bw(base_size = 10, base_family = font_type)
  p = p + theme(legend.key = element_blank())
  # p = p + theme(panel.grid.minor.x = element_line(colour = "grey90", size = 0.2))
  # p = p + theme(panel.grid.minor.y = element_line(colour = "grey90", size = 0.2))
  p = p + theme(panel.grid.minor.x = element_blank())
  p = p + theme(panel.grid.minor.y = element_blank())
  p = p + theme(legend.text.align = 0)
  p = p + theme(panel.border = element_rect(color = "black"))
  p = p + theme(legend.key.width = unit(9, "mm"))
  # p = p + theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5))
  # p = p + theme(plot.margin = unit(c(0,0,0,0), "mm"))
  p = p + ylab(TeX("$\\nu^+")) + xlab(expression(paste(italic(t), " [year]")))
  
  p = p + scale_colour_manual(values = cbPalette, labels = copulas_it)
  p = p + scale_linetype_manual(values = ltPalette, labels = copulas_it)
  p = p + scale_shape_manual(values = shPalette, labels = copulas_it)
}

# raw nu plot
p1 = basic_plot(beta_df[beta_df$id == "raw", ])
p1 = p1 + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), 
                axis.title.x=element_blank())
if (log_plot == 1){
p1 = p1 + scale_y_log10(breaks = c(1e-4, 1e-3, 1e-2), labels = trans_format("log10", math_format(10^.x)))#,
                      #minor_breaks = log10_minor_break())
# p1 = p1 + geom_hline(yintercept = log10space(c(1e-4, 1e-3, 1e-2)))
}

# normalized nu plot
p2 = basic_plot(beta_df[beta_df$id == "normalized", ])
p2 = p2 + ylab(TeX("$\\nu^+ / \\nu^+_\\mathrm{Gauss-DI}$"))

# put them toghether
grobs <- ggplotGrob(p1)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


pp = plot_grid(p1 + theme(legend.position = "none"),
          p2 + theme(legend.position = "none"),
          rel_heights = c(5,6),
          align = "v", nrow = 2)

p = plot_grid(pp, legend, ncol = 2, rel_widths = c(1, 0.4))

print(p)


################################################################################
## SAVE PLOT
################################################################################
if (save_plot == 1){
  setEPS()
  if (log_plot == 1){
    fig_name = "Sudrets_beam_log_nu.pdf" 
  } else {
    fig_name = "Sudrets_beam_nu.pdf" 
  }
    
  # postscript(file = file.path(fig_dir, fig_name), family = 'Times New Roman', title = fig_name,
  #            width = 130/25.4, height = 70/25.4,
  #            colormodel = 'cmyk')
  #   postscript(file = file.path(figure_dir, fig_name), family = 'Times New Roman',
  #            width = 100/25.4, height = 120/25.4)
  # print(g)
  pdf(file = file.path(fig_dir, fig_name), family = font_type,
      width = 130/25.4, height = 90/25.4,
      colormodel = 'cmyk')
  print(p)
  
  dev.off()
  
  embed_fonts(file.path(fig_dir, fig_name))
}