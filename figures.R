# Libraries
library(ggplot2)
library(gganimate)
library(gifski)
library(gapminder)
library(gridExtra)
library(grid)
library(MRColour)
library(plyr)
library(BayesFactor)
library(assortedRFunctions)

# Standard values
minN         <- 10
barWidth     <- 3
panMar1      <- 0.16
panMar2      <- 0.3 
bar_yAxis    <- 0.22
hist_yAxis   <- 100
middle_xlim  <- c(-5.5, 9.5)
crit1        <- 10
crit2        <- 1/6
breaksVal    <- c(-5, -2, 0, 2, 5, 9)
breaksLab    <- c('1/6', '1/3', '1', '3', '6', '10')
nIter        <- 10000
width        <- 20
heigth       <- width/1.111289

# Global settings
fontSize <- 15
theme_set(theme_gray(base_size = fontSize))


# load data
load("simulationResults.RData")

# /* 
# ----------------------------- Student's t ---------------------------
# */
set.seed(911225)
n    <- 60
data <- rnorm(n, 0.0, 1)

df             <- n - 1
standard.error <- sd(data)/sqrt(n)
t.value        <- (mean(data) - 0)/standard.error 

pl <- ggplot(data.frame(x = c(-8, 8)), aes(x = x)) +
  stat_function(fun = dt, args = list(df = df)) +
  geom_vline(xintercept = t.value) +
  geom_ribbon(data = data.frame(x = seq(t.value, 8, 0.01)), mapping = aes(ymin = 0, ymax = dt(seq(t.value, 8, 0.01), df))) +
  geom_ribbon(data = data.frame(x = seq(-8, -t.value, 0.01)), mapping = aes(ymin = 0, ymax = dt(seq(-8, -t.value, 0.01), df))) +
  labs(x = expression(italic("t")*"-value"), y = 'Density', title = expression("Student's "*italic("t")*"-distribution df = 59"))

ggsave(file = "figures/student.png", pl, width = width, height = heigth, units = "cm")

# /* 
# ----------------------------- Traditional design ---------------------------
# */
# Calculatge percentage
df1_2            <- subset(df1, df1$side == 2)
df1_2$overCrit1  <- df1_2$bf > 10
df1_2$belowCrit2 <- df1_2$bf < 1/6
df1_2[df1_2$d == 0.0, 'overCrit1']  <- NA
df1_2[df1_2$d == 0.5, 'belowCrit2'] <- NA
df1_2_agg         <- ddply(df1_2, c('n'), summarise, per_crit1 = mean(overCrit1, na.rm = TRUE), per_crit2 = mean(belowCrit2, na.rm = TRUE))

# Misleading evidence
df1_2_miss            <- df1_2
df1_2_miss$overCrit1  <- df1_2_miss$bf > 10
df1_2_miss$belowCrit2 <- df1_2_miss$bf < 1/6
df1_2_miss[df1_2_miss$d == 0.5, 'overCrit1']  <- NA
df1_2_miss[df1_2_miss$d == 0.0, 'belowCrit2'] <- NA
df1_2_miss_agg         <- ddply(df1_2_miss, c('n'), summarise, per_crit1 = mean(overCrit1, na.rm = TRUE), per_crit2 = mean(belowCrit2, na.rm = TRUE))

# Plot
pl <- ggplot(df1_2_agg, aes(x = n)) +
  geom_hline(yintercept = 0.80) +
  geom_line(aes(y = per_crit1, colour = "BF10 > 10")) + 
  geom_line(aes(y = per_crit2, colour = "BF10 < 1/6")) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 0), 
        legend.position = c(1, 0)) +
  labs(x = 'Sample size', y = 'Percentage')

ggsave(file = "figures/traditional.png", pl, width = width, height = heigth, units = "cm")

# /* 
# ----------------------------- Create animation ---------------------------
# */
set.seed(2)

# Parameters
minN      <- 10
batchSize <- 5
crit1     <- 10
crit2     <- 1/6
maxN      <- 300
nIter     <- 100
d         <- 0

for(i in 1:nIter){
  # First iteration
  data <- rnorm(minN, d, 1)
  n    <- minN
  bf   <- reportBF(ttestBF(data))
  
  # Loop
  while(bf[length(bf)] < crit1 & bf[length(bf)] > crit2){
    data <- c(data, rnorm(minN, d, 1))
    bf[length(bf) + 1] <- reportBF(ttestBF(data))
    n[length(n) + 1]   <- n[length(n)] + batchSize
  }
  
  if(i == 1){
    df <- data.frame(id = rep(i, length(bf)), n = n, bf = bf)
  } else {
    df <- rbind(df, data.frame(id =  rep(i, length(bf)), n = n, bf = bf))
  }
}

df$time <- 1:nrow(df)

# Transform BF
df$trans_bf <- NA
df$trans_bf[df$bf < 1] <- -1/df$bf[df$bf < 1] + 1
df$trans_bf[df$bf > 1] <- df$bf[df$bf > 1] - 1

# Create Agg
df_agg <- ddply(df, c('id'), summarise, n = n[length(n)], bf = bf[length(bf)])
df_agg$support <- 'undecided'
df_agg$support[df_agg$bf > crit1] <- 'H1'
df_agg$support[df_agg$bf < crit2] <- 'H0'

# Creates band
df_agg$band                                     <- '> 10'
df_agg$band[df_agg$bf < 10 & df_agg$bf > 6]     <- '> 6'
df_agg$band[df_agg$bf < 6 & df_agg$bf > 3]      <- '> 3'
df_agg$band[df_agg$bf < 3 & df_agg$bf > 1]      <- '> 1'
df_agg$band[df_agg$bf < 1 & df_agg$bf > 1/3]    <- '< 1'
df_agg$band[df_agg$bf < 1/3 & df_agg$bf > 1/6]  <- '< 1/3'
df_agg$band[df_agg$bf < 1/6 & df_agg$bf > 1/10] <- '< 1/6'
df_agg$band[df_agg$bf < 1/10]                   <- '< 1/10'

# Create factor band
df_agg$band <- factor(df_agg$band, levels = c('> 10', '> 6', '> 3', '> 1', '< 1', '< 1/3', '< 1/6', '< 1/10'))

# Get back to main DF
df$band <- rep(df_agg$band, table(df$id))

# Plot
linePlot <- ggplot(df, aes(x = n, y = trans_bf, group = id, colour = band)) + 
  geom_line(alpha = 0.5, show.legend = FALSE) + 
  scale_colour_BF() +
  geom_hline(yintercept = 9) +
  geom_hline(yintercept = -5) +
  scale_y_continuous(breaks = breaksVal,
                     labels = breaksLab) +
  coord_cartesian(xlim = c(minN, maxN), ylim = middle_xlim, expand = FALSE) +
  labs(y = expression(BF[10]), x = NULL) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(panMar, panMar, panMar, 0), "cm"))

anim <- linePlot + geom_point(show.legend = FALSE) + transition_reveal(time)

animate(anim, height = 300, width = 800)
anim_save("figures/anim.gif")

# /* 
# ----------------------------- Sequential design 1 --------------------------
# */
# For plotting
minN         <- 10
maxN         <- 1050
bar_yAxis    <- 0.25
nIter        <- 10000

tempDF <- subset(df2, side == 2 & d == 0.0)

tempDF$trans_bf <- NA
tempDF$trans_bf[tempDF$bf < 1] <- -1/tempDF$bf[tempDF$bf < 1] + 1
tempDF$trans_bf[tempDF$bf > 1] <- tempDF$bf[tempDF$bf > 1] - 1

tempDF_agg <- ddply(tempDF, c('id'), summarise, n = n[length(n)], bf = bf[length(bf)])
tempDF_agg$support <- 'undecided'
tempDF_agg$support[tempDF_agg$bf > crit1] <- 'H1'
tempDF_agg$support[tempDF_agg$bf < crit2] <- 'H0'

# Creates band
tempDF_agg$band                                             <- '> 10'
tempDF_agg$band[tempDF_agg$bf < 10 & tempDF_agg$bf > 6]     <- '> 6'
tempDF_agg$band[tempDF_agg$bf < 6 & tempDF_agg$bf > 3]      <- '> 3'
tempDF_agg$band[tempDF_agg$bf < 3 & tempDF_agg$bf > 1]      <- '> 1'
tempDF_agg$band[tempDF_agg$bf < 1 & tempDF_agg$bf > 1/3]    <- '< 1'
tempDF_agg$band[tempDF_agg$bf < 1/3 & tempDF_agg$bf > 1/6]  <- '< 1/3'
tempDF_agg$band[tempDF_agg$bf < 1/6 & tempDF_agg$bf > 1/10] <- '< 1/6'
tempDF_agg$band[tempDF_agg$bf < 1/10]                       <- '< 1/10'

# Create factor band
tempDF_agg$band <- factor(tempDF_agg$band, levels = c('> 10', '> 6', '> 3', '> 1', '< 1', '< 1/3', '< 1/6', '< 1/10'))


# Get back to main DF
tempDF$band <- rep(tempDF_agg$band, table(tempDF$id))


# Create upper histogram
tempDF_agg_supp_H1 <- ddply(subset(tempDF_agg, support == 'H1'),
                            c('n', 'band'),
                            summarise,
                            freq = length(bf)/nIter)

# Upper Histogram
upper_hist <- ggplot(tempDF_agg_supp_H1, aes(x = n, y = freq, fill = band)) + 
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.5, width = barWidth)  +
  scale_fill_BF() +
  coord_cartesian(ylim = c(0, bar_yAxis), xlim = c(minN - 0.5, maxN + 0.5), expand = FALSE) +
  labs(y = 'Frequency', x = NULL) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))

# Create line plot
linePlot <- ggplot(tempDF, aes(x = n, y = trans_bf, group = id, colour = band)) + 
  geom_line(alpha = 0.05, show.legend = FALSE) + 
  scale_colour_BF() +
  geom_hline(yintercept = 9) +
  geom_hline(yintercept = -5) +
  scale_y_continuous(breaks = breaksVal,
                     labels = breaksLab) +
  coord_cartesian(xlim = c(minN, maxN), ylim = middle_xlim, expand = FALSE) +
  labs(y = expression(BF[10]), x = NULL) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))

# Create lower histogram
tempDF_agg_supp_H0 <- ddply(subset(tempDF_agg, support == 'H0'),
                            c('n', 'band'),
                            summarise,
                            freq = length(bf)/nIter)

lower_hist <- ggplot(tempDF_agg_supp_H0, aes(x = n, y = freq, fill = band)) + 
  scale_y_reverse() +
  scale_fill_BF(drop = F) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.5, width = barWidth) + 
  labs(y = 'Frequency', x = 'Sample size') +
  coord_cartesian(ylim = c(bar_yAxis, 0), xlim = c(minN -0.5, maxN + 0.5), expand = FALSE) +
  theme(plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))

# Get blank plot
blank <- grid.rect(gp=gpar(col="white"))

# Get legend plot
legendPlot <- ggplot(tempDF_agg_supp_H1, aes(x = n, fill = band)) +
  geom_histogram(alpha = 0.5) +
  scale_fill_BF(drop = FALSE) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_legend(title = expression(BF[10])))

legend <- cowplot::get_legend(legendPlot)


# Combine everything into one plot
## Get Grobs
gplot1 <- ggplotGrob(upper_hist)
gplot2 <- ggplotGrob(linePlot)
gplot4 <- ggplotGrob(lower_hist)

# Align widths
maxWidth = grid::unit.pmax(gplot1$widths, gplot2$widths, gplot4$widths)
gplot1$widths <- as.list(maxWidth)
gplot2$widths <- as.list(maxWidth)
gplot4$widths <- as.list(maxWidth)

lay <- rbind(c(1,2,3),
             c(4,5,3),
             c(6,7,3))

g1 <- arrangeGrob(grobs = list(gplot1, blank, legend, gplot2, blank, gplot4, blank),
                  newpage = FALSE,
                  top =  textGrob(expression("Simulation of Bayesian sequential design analysis (d"[true]*" = 0.0)"), gp = gpar(fontsize = fontSize)),
                  layout_matrix = lay,
                  widths = c(12, 3, 2))

ggsave(file = "figures/pure1.png", g1, width = width, height = heigth, units = "cm")


# /* 
# ----------------------------- Sequential design 2--------------------------
# */
# For plotting
minN         <- 10
maxN         <- 170
bar_yAxis    <- 0.22
hist_yAxis   <- 100


tempDF <- subset(df2, side == 2 & d == 0.5)

tempDF$trans_bf <- NA
tempDF$trans_bf[tempDF$bf < 1] <- -1/tempDF$bf[tempDF$bf < 1] + 1
tempDF$trans_bf[tempDF$bf > 1] <- tempDF$bf[tempDF$bf > 1] - 1

tempDF_agg <- ddply(tempDF, c('id'), summarise, n = n[length(n)], bf = bf[length(bf)])
tempDF_agg$support <- 'undecided'
tempDF_agg$support[tempDF_agg$bf > crit1] <- 'H1'
tempDF_agg$support[tempDF_agg$bf < crit2] <- 'H0'


# Creates band
tempDF_agg$band                                             <- '> 10'
tempDF_agg$band[tempDF_agg$bf < 10 & tempDF_agg$bf > 6]     <- '> 6'
tempDF_agg$band[tempDF_agg$bf < 6 & tempDF_agg$bf > 3]      <- '> 3'
tempDF_agg$band[tempDF_agg$bf < 3 & tempDF_agg$bf > 1]      <- '> 1'
tempDF_agg$band[tempDF_agg$bf < 1 & tempDF_agg$bf > 1/3]    <- '< 1'
tempDF_agg$band[tempDF_agg$bf < 1/3 & tempDF_agg$bf > 1/6]  <- '< 1/3'
tempDF_agg$band[tempDF_agg$bf < 1/6 & tempDF_agg$bf > 1/10] <- '< 1/6'
tempDF_agg$band[tempDF_agg$bf < 1/10]                       <- '< 1/10'

# Create factor band
tempDF_agg$band <- factor(tempDF_agg$band, levels = c('> 10', '> 6', '> 3', '> 1', '< 1', '< 1/3', '< 1/6', '< 1/10'))


# Get back to main DF
tempDF$band <- rep(tempDF_agg$band, table(tempDF$id))


# Create upper histogram
tempDF_agg_supp_H1 <- ddply(subset(tempDF_agg, support == 'H1'),
                            c('n', 'band'),
                            summarise,
                            freq = length(bf)/nIter)

# Upper Histogram
upper_hist <- ggplot(tempDF_agg_supp_H1, aes(x = n, y = freq, fill = band)) + 
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.5, width = barWidth)  +
  scale_fill_BF() +
  coord_cartesian(ylim = c(0, bar_yAxis), xlim = c(minN - 0.5, maxN + 0.5), expand = FALSE) +
  labs(y = 'Frequency', x = NULL) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))

# Create line plot
linePlot <- ggplot(tempDF, aes(x = n, y = trans_bf, group = id, colour = band)) + 
  geom_line(alpha = 0.05, show.legend = FALSE) + 
  scale_colour_BF() +
  geom_hline(yintercept = 9) +
  geom_hline(yintercept = -5) +
  scale_y_continuous(breaks = breaksVal,
                     labels = breaksLab) +
  coord_cartesian(xlim = c(minN, maxN), ylim = middle_xlim, expand = FALSE) +
  labs(y = expression(BF[10]), x = NULL) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))


# Create lower histogram
tempDF_agg_supp_H0 <- ddply(subset(tempDF_agg, support == 'H0'),
                            c('n', 'band'),
                            summarise,
                            freq = length(bf)/nIter)



# Creates empty df for plotting if no values
if(nrow(tempDF_agg_supp_H0) == 0){
  tempDF_agg_supp_H0 <- data.frame(n = seq(1, 5, 1),
                                   band = rep('< 1/10', 5),
                                   freq = rep(0, 5))
}


lower_hist <- ggplot(tempDF_agg_supp_H0, aes(x = n, y = freq, fill = band)) + 
  scale_y_reverse() +
  scale_fill_BF(drop = F) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.5, width = barWidth) + 
  labs(y = 'Frequency', x = 'Sample size') +
  coord_cartesian(ylim = c(bar_yAxis, 0), xlim = c(minN -0.5, maxN + 0.5), expand = FALSE) +
  theme(plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))

# Get blank plot
blank <- grid.rect(gp=gpar(col="white"))

# Get legend plot
legendPlot <- ggplot(tempDF_agg_supp_H1, aes(x = n, fill = band)) +
  geom_histogram(alpha = 0.5) +
  scale_fill_BF(drop = FALSE) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_legend(title = expression(BF[10])))

legend <- cowplot::get_legend(legendPlot)


# Combine everything into one plot
## Get Grobs
gplot1 <- ggplotGrob(upper_hist)
gplot2 <- ggplotGrob(linePlot)
#gplot3 <- ggplotGrob(right_hist) # not adjusted
gplot4 <- ggplotGrob(lower_hist)

# Align widths
maxWidth = grid::unit.pmax(gplot1$widths, gplot2$widths, gplot4$widths)
gplot1$widths <- as.list(maxWidth)
gplot2$widths <- as.list(maxWidth)
gplot4$widths <- as.list(maxWidth)

# Align heights
#maxHeight = grid::unit.pmax(gplot2$heights, gplot3$heights)
#gplot2$heights <- as.list(maxHeight)
#gplot3$heights <- as.list(maxHeight)

lay <- rbind(c(1,2,3),
             c(4,5,3),
             c(6,7,3))

g2 <- arrangeGrob(grobs = list(gplot1, blank, legend, gplot2, blank, gplot4, blank),
                  newpage = FALSE,
                  top =  textGrob(expression("Simulation of Bayesian sequential design analysis (d"[true]*" = 0.5)"), gp = gpar(fontsize = fontSize)),
                  layout_matrix = lay,
                  widths = c(12, 3, 2))

ggsave(file = "figures/pure2.png", g2, width = width, height = heigth, units = "cm")
# /* 
# ----------------------------- Sequential design with limit 1---------------------------
# */
# For plotting
minN         <- 10
maxN         <- 100
bar_yAxis    <- 0.22
hist_yAxis   <- bar_yAxis*nIter

tempDF <- subset(df3, d == 0.0)

tempDF$trans_bf <- NA
tempDF$trans_bf[tempDF$bf < 1] <- -1/tempDF$bf[tempDF$bf < 1] + 1
tempDF$trans_bf[tempDF$bf > 1] <- tempDF$bf[tempDF$bf > 1] - 1

tempDF_agg <- ddply(tempDF, c('id'), summarise, n = n[length(n)], bf = bf[length(bf)])
tempDF_agg$support <- 'undecided'
tempDF_agg$support[tempDF_agg$bf > crit1] <- 'H1'
tempDF_agg$support[tempDF_agg$bf < crit2] <- 'H0'

# Creates band
tempDF_agg$band                                             <- '> 10'
tempDF_agg$band[tempDF_agg$bf < 10 & tempDF_agg$bf > 6]     <- '> 6'
tempDF_agg$band[tempDF_agg$bf < 6 & tempDF_agg$bf > 3]      <- '> 3'
tempDF_agg$band[tempDF_agg$bf < 3 & tempDF_agg$bf > 1]      <- '> 1'
tempDF_agg$band[tempDF_agg$bf < 1 & tempDF_agg$bf > 1/3]    <- '< 1'
tempDF_agg$band[tempDF_agg$bf < 1/3 & tempDF_agg$bf > 1/6]  <- '< 1/3'
tempDF_agg$band[tempDF_agg$bf < 1/6 & tempDF_agg$bf > 1/10] <- '< 1/6'
tempDF_agg$band[tempDF_agg$bf < 1/10]                       <- '< 1/10'

# Create factor band
tempDF_agg$band <- factor(tempDF_agg$band, levels = c('> 10', '> 6', '> 3', '> 1', '< 1', '< 1/3', '< 1/6', '< 1/10'))


# Get back to main DF
tempDF$band <- rep(tempDF_agg$band, table(tempDF$id))


# Create upper histogram
tempDF_agg_supp_H1 <- ddply(subset(tempDF_agg, support == 'H1'),
                            c('n', 'band'),
                            summarise,
                            freq = length(bf)/nIter)

# Upper Histogram
upper_hist <- ggplot(tempDF_agg_supp_H1, aes(x = n, y = freq, fill = band)) + 
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.5, width = barWidth)  +
  scale_fill_BF() +
  coord_cartesian(ylim = c(0, bar_yAxis), xlim = c(minN - 0.5, maxN + 0.5), expand = FALSE) +
  labs(y = 'Frequency', x = NULL) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))

# Create line plot
linePlot <- ggplot(tempDF, aes(x = n, y = trans_bf, group = id, colour = band)) + 
  geom_line(alpha = 0.05, show.legend = FALSE) + 
  scale_colour_BF() +
  geom_hline(yintercept = 9) +
  geom_hline(yintercept = -5) +
  scale_y_continuous(breaks = breaksVal,
                     labels = breaksLab) +
  coord_cartesian(xlim = c(minN, maxN), ylim = middle_xlim, expand = FALSE) +
  labs(y = expression(BF[10]), x = NULL) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))

# Create right histogram
tempDF_agg_undecided <- subset(tempDF_agg, tempDF_agg$n == maxN & tempDF_agg$bf < crit1 & tempDF_agg$bf > crit2)
tempDF_agg_undecided$trans_bf <- NA
tempDF_agg_undecided$trans_bf[tempDF_agg_undecided$bf < 1] <- -1/tempDF_agg_undecided$bf[tempDF_agg_undecided$bf < 1] + 1
tempDF_agg_undecided$trans_bf[tempDF_agg_undecided$bf > 1] <- tempDF_agg_undecided$bf[tempDF_agg_undecided$bf > 1] - 1

right_hist <- ggplot(tempDF_agg_undecided, aes(x = trans_bf, fill = band)) +
  coord_flip(ylim = c(0, hist_yAxis), xlim = middle_xlim, expand = FALSE) +
  geom_histogram(alpha = 0.5, show.legend = FALSE) +
  scale_fill_BF(drop = FALSE) +
  labs(y = 'Count', x = NULL) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(panMar1, panMar1, panMar1, 0), "cm"))

# Create lower histogram
tempDF_agg_supp_H0 <- ddply(subset(tempDF_agg, support == 'H0'),
                            c('n', 'band'),
                            summarise,
                            freq = length(bf)/nIter)



# Creates empty df for plotting if no values
if(nrow(tempDF_agg_supp_H0) == 0){
  tempDF_agg_supp_H0 <- data.frame(n = seq(1, 5, 1),
                                   band = rep('< 1/10', 5),
                                   freq = rep(0, 5))
}


lower_hist <- ggplot(tempDF_agg_supp_H0, aes(x = n, y = freq, fill = band)) + 
  scale_y_reverse() +
  scale_fill_BF(drop = F) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.5, width = barWidth) + 
  labs(y = 'Frequency', x = 'Sample size') +
  coord_cartesian(ylim = c(bar_yAxis, 0), xlim = c(minN -0.5, maxN + 0.5), expand = FALSE) +
  theme(plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))

# Get blank plot
blank <- grid.rect(gp=gpar(col="white"))

# Get legend plot
legendPlot <- ggplot(tempDF_agg_supp_H1, aes(x = n, fill = band)) +
  geom_histogram(alpha = 0.5) +
  scale_fill_BF(drop = FALSE) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_legend(title = expression(BF[10])))

legend <- cowplot::get_legend(legendPlot)


# Combine everything into one plot
## Get Grobs
gplot1 <- ggplotGrob(upper_hist)
gplot2 <- ggplotGrob(linePlot)
gplot3 <- ggplotGrob(right_hist) # not adjusted
gplot4 <- ggplotGrob(lower_hist)

# Align widths
maxWidth = grid::unit.pmax(gplot1$widths, gplot2$widths, gplot4$widths)
gplot1$widths <- as.list(maxWidth)
gplot2$widths <- as.list(maxWidth)
gplot4$widths <- as.list(maxWidth)

# Align heights
maxHeight = grid::unit.pmax(gplot2$heights, gplot3$heights)
gplot2$heights <- as.list(maxHeight)
gplot3$heights <- as.list(maxHeight)

lay <- rbind(c(1,2,3),
             c(4,5,3),
             c(6,7,3))

g1 <- arrangeGrob(grobs = list(gplot1, blank, legend, gplot2, gplot3, gplot4, blank),
             newpage = FALSE,
             top =  textGrob(expression("Simulation of Bayesian sequential design analysis (d"[true]*" = 0.0)"), gp = gpar(fontsize = fontSize)),
             layout_matrix = lay,
             widths = c(12, 3, 2))


ggsave(file = "figures/limit1.png", g1, width = width, height = heigth, units = "cm")

# /* 
# ----------------------------- Sequential design with limit 2---------------------------
# */
# For plotting
minN         <- 10
bar_yAxis    <- 0.22
hist_yAxis   <- bar_yAxis*nIter

tempDF <- subset(df3, d == 0.5)

tempDF$trans_bf <- NA
tempDF$trans_bf[tempDF$bf < 1] <- -1/tempDF$bf[tempDF$bf < 1] + 1
tempDF$trans_bf[tempDF$bf > 1] <- tempDF$bf[tempDF$bf > 1] - 1

tempDF_agg <- ddply(tempDF, c('id'), summarise, n = n[length(n)], bf = bf[length(bf)])
tempDF_agg$support <- 'undecided'
tempDF_agg$support[tempDF_agg$bf > crit1] <- 'H1'
tempDF_agg$support[tempDF_agg$bf < crit2] <- 'H0'

# Creates band
tempDF_agg$band                                             <- '> 10'
tempDF_agg$band[tempDF_agg$bf < 10 & tempDF_agg$bf > 6]     <- '> 6'
tempDF_agg$band[tempDF_agg$bf < 6 & tempDF_agg$bf > 3]      <- '> 3'
tempDF_agg$band[tempDF_agg$bf < 3 & tempDF_agg$bf > 1]      <- '> 1'
tempDF_agg$band[tempDF_agg$bf < 1 & tempDF_agg$bf > 1/3]    <- '< 1'
tempDF_agg$band[tempDF_agg$bf < 1/3 & tempDF_agg$bf > 1/6]  <- '< 1/3'
tempDF_agg$band[tempDF_agg$bf < 1/6 & tempDF_agg$bf > 1/10] <- '< 1/6'
tempDF_agg$band[tempDF_agg$bf < 1/10]                       <- '< 1/10'

# Create factor band
tempDF_agg$band <- factor(tempDF_agg$band, levels = c('> 10', '> 6', '> 3', '> 1', '< 1', '< 1/3', '< 1/6', '< 1/10'))


# Get back to main DF
tempDF$band <- rep(tempDF_agg$band, table(tempDF$id))


# Create upper histogram
tempDF_agg_supp_H1 <- ddply(subset(tempDF_agg, support == 'H1'),
                            c('n', 'band'),
                            summarise,
                            freq = length(bf)/nIter)

# Upper Histogram
upper_hist <- ggplot(tempDF_agg_supp_H1, aes(x = n, y = freq, fill = band)) + 
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.5, width = barWidth)  +
  scale_fill_BF() +
  coord_cartesian(ylim = c(0, bar_yAxis), xlim = c(minN - 0.5, maxN + 0.5), expand = FALSE) +
  labs(y = 'Frequency', x = NULL) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))

# Create line plot
linePlot <- ggplot(tempDF, aes(x = n, y = trans_bf, group = id, colour = band)) + 
  geom_line(alpha = 0.05, show.legend = FALSE) + 
  scale_colour_BF() +
  geom_hline(yintercept = 9) +
  geom_hline(yintercept = -5) +
  scale_y_continuous(breaks = breaksVal,
                     labels = breaksLab) +
  coord_cartesian(xlim = c(minN, maxN), ylim = middle_xlim, expand = FALSE) +
  labs(y = expression(BF[10]), x = NULL) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))

# Create right histogram
tempDF_agg_undecided <- subset(tempDF_agg, tempDF_agg$n == maxN & tempDF_agg$bf < crit1 & tempDF_agg$bf > crit2)
tempDF_agg_undecided$trans_bf <- NA
tempDF_agg_undecided$trans_bf[tempDF_agg_undecided$bf < 1] <- -1/tempDF_agg_undecided$bf[tempDF_agg_undecided$bf < 1] + 1
tempDF_agg_undecided$trans_bf[tempDF_agg_undecided$bf > 1] <- tempDF_agg_undecided$bf[tempDF_agg_undecided$bf > 1] - 1

right_hist <- ggplot(tempDF_agg_undecided, aes(x = trans_bf, fill = band)) +
  coord_flip(ylim = c(0, hist_yAxis), xlim = middle_xlim, expand = FALSE) +
  geom_histogram(alpha = 0.5, show.legend = FALSE) +
  scale_fill_BF(drop = FALSE) +
  labs(y = 'Count', x = NULL) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(panMar1, panMar1, panMar1, panMar1), "cm"))

# Create lower histogram
tempDF_agg_supp_H0 <- ddply(subset(tempDF_agg, support == 'H0'),
                            c('n', 'band'),
                            summarise,
                            freq = length(bf)/nIter)



# Creates empty df for plotting if no values
if(nrow(tempDF_agg_supp_H0) == 0){
  tempDF_agg_supp_H0 <- data.frame(n = seq(1, 5, 1),
                                   band = rep('< 1/10', 5),
                                   freq = rep(0, 5))
}


lower_hist <- ggplot(tempDF_agg_supp_H0, aes(x = n, y = freq, fill = band)) + 
  scale_y_reverse() +
  scale_fill_BF(drop = F) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.5, width = barWidth) + 
  labs(y = 'Frequency', x = 'Sample size') +
  coord_cartesian(ylim = c(bar_yAxis, 0), xlim = c(minN -0.5, maxN + 0.5), expand = FALSE) +
  theme(plot.margin = unit(c(panMar1, panMar2, panMar1, 0), "cm"))

# Get blank plot
blank <- grid.rect(gp=gpar(col="white"))

# Get legend plot
legendPlot <- ggplot(tempDF_agg_supp_H1, aes(x = n, fill = band)) +
  geom_histogram(alpha = 0.5) +
  scale_fill_BF(drop = FALSE) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_legend(title = expression(BF[10])))

legend <- cowplot::get_legend(legendPlot)


# Combine everything into one plot
## Get Grobs
gplot1 <- ggplotGrob(upper_hist)
gplot2 <- ggplotGrob(linePlot)
gplot3 <- ggplotGrob(right_hist) # not adjusted
gplot4 <- ggplotGrob(lower_hist)

# Align widths
maxWidth = grid::unit.pmax(gplot1$widths, gplot2$widths, gplot4$widths)
gplot1$widths <- as.list(maxWidth)
gplot2$widths <- as.list(maxWidth)
gplot4$widths <- as.list(maxWidth)

# Align heights
maxHeight = grid::unit.pmax(gplot2$heights, gplot3$heights)
gplot2$heights <- as.list(maxHeight)
gplot3$heights <- as.list(maxHeight)

lay <- rbind(c(1,2,3),
             c(4,5,3),
             c(6,7,3))

g2 <- arrangeGrob(grobs = list(gplot1, blank, legend, gplot2, gplot3, gplot4, blank),
                  newpage = FALSE,
                  top =  textGrob(expression("Simulation of Bayesian sequential design analysis (d"[true]*" = 0.5)"), gp = gpar(fontsize = fontSize)),
                  layout_matrix = lay,
                  widths = c(12, 3, 2))


ggsave(file = "figures/limit2.png", g2, width = width, height = heigth, units = "cm")
