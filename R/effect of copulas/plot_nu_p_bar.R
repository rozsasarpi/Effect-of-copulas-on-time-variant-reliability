rm(list=ls(all=TRUE))

library(Cairo)
library(ggplot2)
library(RColorBrewer)
library(XLConnect)
library(reshape2)
library(dplyr)
library(magrittr)
library(extrafont)
library(RColorBrewer)
library(latex2exp)

# loadfonts(device="postscript")
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.19/bin/gswin64c.exe")

main_dir    = "d:/Working folder/R working folder/effect of copulas"
data_dir    = "d:/Working folder/Matlab working folder/time variant reliability/phi2 method/effect of copulas/stochastic_snow_model"
fig_dir     = file.path(main_dir, "figures")

excel_file  = file.path(data_dir, "stoch_snow_results.xlsx")
sheet_name  = paste("results2")

save_plot   = 1
# font_type   = "Times"
font_type   = "CM Roman"

################################################################################
## LOAD & PREPROCESS DATA
################################################################################
# Gauss - autocorrelation function
wb          = loadWorkbook(excel_file)
nu_p_df     = readWorksheet(wb, sheet = sheet_name, region = "C7:E12", useCachedValues = TRUE)

marginals   = c("Gauss", "Lognormal", "Gumbel")
copulas     = c("Gauss-DI", "t (dof=2)-DI", "Gumbel-DI", "rotGumbel-180°-DI", "rotClayton-180°-DI", "tev-DI")

# for legend label
copulas_it     = c("Gauss-DI", expression(paste(italic("t"), " (dof=2)-DI")), "Gumbel-DI", "rotGumbel-180°-DI", "rotClayton-180°-DI", expression(paste(italic("tev"), "-DI")))

names(nu_p_df) = marginals
#!!!!!!!!!!!!!!!!!!!!!!!!!!!
# remove the tev copula that distrots data
ridx        = c(6)
nu_p_df     = nu_p_df[-ridx,]
copulas     = copulas[-ridx]
copulas_it  = copulas_it[-ridx]
#!!!!!!!!!!!!!!!!!!!!!!!!!!!

# normalized out-crossing rates
wg          = sapply(nu_p_df, function(x) x/x[1]) %>% data.frame()

wg$acorr    = factor("Gauss")
wg$copula   = as.factor(copulas)


# combine values into a single dataframe
w           = wg
w_df        = melt(w)
names(w_df) = c("acorr", "copula","marginal", "w")

# convert ID-s to factors and change the order of specimens for plot to match coloring
w_df = within(w_df, copula <- factor(copula, levels = as.factor(copulas)))
w_df = within(w_df, marginal <- factor(marginal, levels = as.factor(marginals)))

# set colors to copula types
nc          = length(copulas)
cbSet1      = brewer.pal(nc, "Set1")
cbPalette   = cbSet1

################################################################################
## PLOT
###############################################################################
# df          = filter(w_df, acorr == "Gauss")
df = w_df

p = ggplot(df, aes(x = copula, y = w, fill = copula)) + geom_bar(stat="identity")
p = p + facet_wrap("marginal")
p = p + scale_fill_manual(values = cbPalette, labels = copulas_it)
# p = p + scale_y_continuous(breaks = c(0, 0.5, 1.0))
# p = p + ylab(expression(paste("Normalized out-crossing rate")))
p = p + ylab(TeX("$\\nu^+ / \\nu^+_\\mathrm{Gauss}$"))
p = p + theme_classic(base_size = 10, base_family = font_type)
p = p + theme(axis.line.y = element_blank())
p = p + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
p = p + theme(strip.text = element_text(color = "black"))
p = p + theme(legend.text.align = 0)
# p = p + theme(panel.grid.major = element_line(color = "gray"))
# p = p + theme(panel.grid.major.x = element_blank())
p = p + geom_hline(yintercept=seq(1, 3, 1), col="white", lwd=1)
print(p)


################################################################################
## SAVE PLOT
################################################################################
if (save_plot == 1){
  # setEPS()

    fig_name = "normalized_nu_p.pdf"
    # fig_name = "normalized_nu_p.png"
    # postscript(file = file.path(fig_dir, fig_name), family = font_type, title = fig_name,
    #          width = 160/25.4, height = 50/25.4,
    #          colormodel = 'cmyk')
    
    # cairo_ps(filename = file.path(fig_dir, fig_name), family = font_type,
    #          width = 160/25.4, height = 50/25.4)
    pdf(file = file.path(fig_dir, fig_name), family = font_type,
        width = 160/25.4, height = 50/25.4,
        colormodel = 'cmyk')
    # png(filename = file.path(fig_dir, fig_name), family = font_type,
    #     width = 160/25.4, height = 50/25.4, units = "in",
    #     res = 1200)
  print(p)
  
  dev.off()
  
  embed_fonts(file.path(fig_dir, fig_name))
}
