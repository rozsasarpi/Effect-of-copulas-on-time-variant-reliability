# The effect of copulas on time-variant reliability 
# with continuous stochastic processes
#
# Simple, time-variant example
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

copulas     = c("Gauss-DI", "t (dof=2)-DI", "Gumbel-DI", "rotGumbel-180°-DI", "rotClayton-180°-DI")
copulas_excel = c("Gauss-DI", "t (dof=2)-DI", "Gumbel-DI", "rotClayton-180°-DI", "rotGumbel-180°-DI")

# for legend label
copulas_it     = c("Gauss-DI", expression(paste(italic("t"), " (dof=2)-DI")), "Gumbel-DI", "rotGumbel-180°-DI", "rotClayton-180°-DI")

# set colors to copula types
cbSet1      = brewer.pal(5, "Set1")
cbPalette   = cbSet1
shPalette   = c(1+15, 2+15, 0+15, 3, 5+13) 

save_plot   = 1
# font_type   = "Times"
font_type   = "CM Roman"

main_dir    = "d:/Working folder/R working folder/effect of copulas"
data_dir    = "d:/Working folder/Matlab working folder/time variant reliability/phi2 method/effect of copulas/simple example"
fig_dir     = file.path(main_dir, "figures")

excel_file  = file.path(data_dir, "simple_example_results_gauss.xlsx")
sheet_name  = paste("summary")

################################################################################
## LOAD DATA
################################################################################
wb = loadWorkbook(excel_file)
beta_df = readWorksheet(wb, sheet = sheet_name, region = "K2:O10", useCachedValues = TRUE)

################################################################################
## PRE-PROCESS DATA
################################################################################
names(beta_df) = signif(beta_df[1,],2)
beta_df = beta_df[-c(1,2),]
# remove Clayton copula - not working - no results
beta_df = beta_df[-3,]

# @@@@@@@@@@@@@@@@@@@@@
# convert to failure probability
# col_names = names(beta_df)
# beta_df = data.frame(pnorm(-as.matrix(beta_df)))
# names(beta_df) = col_names
# @@@@@@@@@@@@@@@@@@@@@


tmp_df = beta_df
for (ii in 1:(dim(beta_df)[2])){
  tmp_df[ii,] = (tmp_df[ii,]/beta_df[1,])
  # @@@@@@@@@@@@@@@@@@@@@
  # tmp_df[ii,] = 1/(tmp_df[ii,]/beta_df[1,])
  # @@@@@@@@@@@@@@@@@@@@@
}

beta_df$copula = as.factor(copulas_excel)
beta_df$id     = as.factor("raw")
tmp_df$copula  = as.factor(copulas_excel)
tmp_df$id      = as.factor("normalized")
beta_df        = rbind(beta_df, tmp_df)

beta_df = melt(beta_df)
names(beta_df) = c("copula", "id", "P_f0", "beta")
beta_df$P_f0 = as.numeric(levels(beta_df$P_f0))[beta_df$P_f0]

# convert ID-s to factors and change the order of specimens for plot
beta_df = within(beta_df, copula <- factor(copula, levels = as.factor(copulas)))


################################################################################
## PLOT
################################################################################
basic_plot = function(df, log_y = F) {
  p = ggplot(df, aes(x = P_f0, y = beta, color = copula, shape = copula))
  p = p + geom_path()
  p = p + geom_point(color = "white", size = 2)
  p = p + geom_point(size = 1)
  p = p + scale_x_log10(breaks = unique(beta_df$P_f0), labels = trans_format("log10", math_format(10^.x)),
                        minor_breaks = log10_minor_break())
  if (log_y == T){
    p = p + scale_y_log10()
  }
  
  # p = p + facet_grid(id~., scales = "free_y")
  p = p + theme_bw(base_size = 10, base_family = font_type)
  p = p + theme(legend.key = element_blank())
  p = p + theme(panel.grid.minor.x = element_line(colour = "grey90", size = 0.2))
  p = p + theme(legend.text.align = 0)
  p = p + theme(legend.key.width = unit(9, "mm"))
  p = p + theme(panel.border = element_rect(color = "black"))
  # p = p + theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5))
  # p = p + theme(plot.margin = unit(c(0,0,0,0), "mm"))
  
  p = p + ylab(bquote(beta)) + xlab(TeX("$\\mathit{P}_{\\mathrm{f}}(0)$"))
  
  # @@@@@@@@@@@@@@@@@@@@@
  # p = p + ylab(bquote(italic(P)[f])) + xlab(bquote(italic(P)[f0]))
  # @@@@@@@@@@@@@@@@@@@@@
  
  p = p + scale_colour_manual(values = cbPalette, labels = copulas_it)
  p = p + scale_shape_manual(values = shPalette, labels = copulas_it)
}

# raw beta plot
p1 = basic_plot(beta_df[beta_df$id == "raw", ])
# @@@@@@@@@@@@@@@@@@@@@
# p1 = basic_plot(beta_df[beta_df$id == "raw", ], T)
# @@@@@@@@@@@@@@@@@@@@@
p1 = p1 + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), 
                axis.title.x=element_blank())

# normalized beta plot
p2 = basic_plot(beta_df[beta_df$id == "normalized", ])
p2 = p2 + ylab(expression(paste(italic(beta), "/", italic(beta)[Gauss])))
# @@@@@@@@@@@@@@@@@@@@@
# p2 = p2 + ylab(expression(paste(italic(P)[f.Gauss], "/", italic(P)[f])))
# @@@@@@@@@@@@@@@@@@@@@

# put them toghether
grobs <- ggplotGrob(p1)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


pp = plot_grid(p1 + theme(legend.position="none"),
          p2 + theme(legend.position="none"),
          rel_heights = c(5,6),
          align = "v", nrow = 2)

p = plot_grid(pp, legend, ncol = 2, rel_widths = c(1, 0.4))

print(p)


################################################################################
## SAVE PLOT
################################################################################
if (save_plot == 1){
  # setEPS()
  fig_name = "simple_example_beta.pdf"
  # @@@@@@@@@@@@@@@@@@@@@
  # fig_name = "simple_example_Pf.pdf"
  # @@@@@@@@@@@@@@@@@@@@@
  # postscript(file = file.path(fig_dir, fig_name), family = 'Times New Roman', title = fig_name,
  #            width = 130/25.4, height = 70/25.4,
  #            colormodel = 'cmyk')
  #   postscript(file = file.path(figure_dir, fig_name), family = 'Times New Roman',
  #            width = 100/25.4, height = 120/25.4)
  pdf(file = file.path(fig_dir, fig_name), family = font_type,
      width = 130/25.4, height = 90/25.4,
      colormodel = 'cmyk')
  print(p)
  
  dev.off()
  
  embed_fonts(file.path(fig_dir, fig_name))
}