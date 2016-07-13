
## Combine italic with upright font in legend
library(ggplot2)

dsamp <- diamonds[sample(nrow(diamonds), 1000), ] 

d <- qplot(carat, price, data=dsamp, colour=clarity)

lab <- levels(dsamp$clarity)

lab_italics <- c(expression(paste(italic("Text1"),"rr")), expression(italic("Text2")), lab[3:length(lab)])

d + scale_colour_hue(breaks=lab, labels=lab_italics)

## LaTeX to R
library(latex2exp)
x <- seq(0, 4, length.out=100)
alpha <- 1:5

plot(x, xlim=c(0, 4), ylim=c(0, 10), 
     xlab='x', ylab=TeX('$\\alpha  x^\\alpha$, where $\\alpha \\in 1\\ldots 5$'), 
     type='n', main=TeX('Using $\\LaTeX$ for plotting in base graphics!'))

invisible(sapply(alpha, function(a) lines(x, a*x^a, col=a)))

legend('topleft', legend=TeX(sprintf("$\\alpha = %d$", alpha)), 
       lwd=1, col=alpha)



## tikzDevice might be much better



##

mod_label = c("2seater" = "a",
              "compact" = "b",
              "midsize" = "c",
              "minivan" = "d",
              "pickup"  = "e",
              "subcompact" = "t",
              "suv"     = "g")

df = mpg
df$class = as.factor(df$class)
# tp = unique(df[,c('class')])
tp = as.data.frame(levels(df$class))
names(tp) = "class"
tp$displ = tp$hwy = 1

g = ggplot(df, aes(x = displ, y = hwy))
g = g + geom_point()
g = g + geom_rect(data = tp, aes(fill = class), xmin = -Inf, xmax = Inf,
          ymin = -Inf, ymax = Inf, alpha = 0.3)
# g = g + facet_wrap(~class, labeller = as_labeller(mod_label))
g = g + facet_wrap(~class)
g = g + theme_bw()
g = g + theme(strip.text = element_text(color = c("red", "blue")))
print(g)

##

tp <- unique(tips[,c('sex','day')])
tp$total_bill <- tp$tip <- 1

#Just Fri
ggplot(tips, aes(x = total_bill, y = tip/total_bill)) + 
  geom_rect(data = subset(tp,day == 'Fri'), aes(fill = day), xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3) +
  geom_point(shape=1) + 
  facet_grid(~day) 
