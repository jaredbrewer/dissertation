# Jared Brewer
# Created: 05 November 2020
# Last Edited: 10 November 2022
# RT-qPCR Analysis 

library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(reshape2)

# Read in a file formatted in a wide form - you're best off doing this manually in Excel or in a text editor.
# Row names: Samples
# Col names: Gene Targets

today <- format(Sys.time(), "%m%d%y")
thp1 <- read.csv("./gmtb_inca6_qpcr_all.csv", header = T)

# Do some simple math - it appends them to the end (very important).
thp1$dVEGFA <- thp1$VEGFA - thp1$GAPDH
thp1.agg <- aggregate(dVEGFA ~ cond + rep + exp, thp1, mean)

for (exp in thp1.agg$exp[!duplicated(thp1.agg$exp)]) {
  ref <- mean(thp1.agg[thp1.agg$exp == exp & thp1.agg$cond == "dmso_ctrl",]$dVEGFA)
  thp1.agg$ref[thp1.agg$exp == exp] <- ref
}

thp1.agg$ddCt <- thp1.agg$dVEGFA - thp1.agg$ref
thp1.agg$rq <- 2**-(thp1.agg$ddCt)

thp1.agg <- separate(thp1.agg, cond, into = c("drug", "bacteria"), sep = "_", remove = T)

labs <- c(expression(paste(bold("Control"))), 
          expression(paste(bolditalic("γMtb"))))

for (exp in thp1.agg$exp[!duplicated(thp1.agg$exp)]) {
  subs <- thp1.agg[thp1.agg$exp == exp,]
  thp1.aov <- aov(subs$rq ~ subs$drug*subs$bacteria)
  print(TukeyHSD(thp1.aov))
  if (exp == 1) {
    lims <- c(13, 14.5)
  }
  if (exp == 2) {
    lims <- c(11, 12.5)
  }
  if (exp == 3) {
    lims <- c(13, 14.5)
  }
  thp1.plot <- ggplot(subs, aes(x = drug, y = rq, fill = bacteria)) + scale_y_continuous(limits = c(0, 15)) +
    geom_beeswarm(aes(shape = drug, size = 5), dodge.width = 0.9, show.legend = F) +
    scale_shape_manual(values = c(21, 22, 24)) +
    geom_boxplot(aes(fill = bacteria), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9)) +
    stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
    xlab("Drug Treatment") + ylab("RQ VEGFA") +
    geom_signif(y_position = lims, xmin = c(0.8, 1.2), xmax = c(1.2, 2.2), annotations = c("p < 0.001", "p < 0.001"), textsize = 5.5, color = "black") +
    scale_x_discrete(limits = c("dmso", "inca"), labels = c("DMSO", "INCA-6")) +
    scale_fill_manual(name = "Bacteria", labels = labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") +
    guides(color = "none") + theme_minimal() +
    theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom")
  
  ggsave(paste("THP-1_VEGFA_qPCR_exp", exp, "_", today, ".png", sep = ""), thp1.plot, width = 5.15, height = 6.2, dpi = 300)

  # This code generates bar plots -- you need a single point to plot that represents all of the various images and then data to populate the error bars.
  
  sem <- ddply(subs, c("drug", "bacteria"), summarize, mean = mean(rq), sem = sd(rq)/sqrt(length(rq)))
  sem <- transform(sem, lower = mean-sem, upper = mean+sem)
  colnames(sem) <- c("drug", "bacteria", "rq", "sem", "lower", "upper")

  thp1.col <- ggplot(sem, aes(x = drug, y = rq, fill = bacteria)) + scale_y_continuous(limits = c(0, 15)) +
    geom_col(position = position_dodge(1)) +
    geom_errorbar(data = sem, aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(1), size = 1.5, color = "black") +
    xlab("Treatment") + ylab("RQ VEGFA") + # scale_y_continuous(trans='log10') +
    scale_x_discrete(limits = c("dmso", "inca"), labels = c("DMSO", "INCA-6")) +
    # c("-INCA-6 \n" mmtb, "+INCA-6, \n -Mtb", "-INCA-6, \n +Mtb", "+INCA-6, \n +Mtb")
    scale_fill_manual(name = "Bacteria", labels = labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") +
    geom_signif(y_position = lims, xmin = c(0.8, 1.2), xmax = c(1.2, 2.2), annotations = c("p < 0.001", "p < 0.001"), textsize = 5.5, color = "black") +
    # guides(color = "none") +
    # dark_mode() +
    theme_minimal() +
    theme(text = element_text(size = 20, face = "bold")) + theme(legend.position="bottom")
  
  ggsave(paste("THP-1_VEGFA_qPCR_bar_exp", exp, "_", today, ".png", sep = ""), thp1.col, width = 5.15, height = 6.2, dpi = 300)
}
 
thp1.dagg <- aggregate(rq ~ drug + bacteria + exp, thp1.agg, mean)

thp1.aov <- aov(rq ~ drug * bacteria, data = thp1.dagg)
TukeyHSD(thp1.aov)

lims = c(12, 14, 10)

thp1.plot <- ggplot(thp1.dagg, aes(x = drug, y = rq, fill = bacteria)) + # scale_y_continuous(limits = c(0,15)) +
  geom_beeswarm(aes(shape = drug, size = 5), dodge.width = 0.9, show.legend = F) +
  scale_shape_manual(values = c(21, 22, 24)) +
  geom_boxplot(aes(fill = bacteria), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Drug Treatment") + ylab("RQ VEGFA") + 
  geom_signif(y_position = lims, xmin = c(0.8, 1.2, 1.8), xmax = c(1.2, 2.2, 2.2), annotations = c("p < 0.001", "p < 0.001", "p = 0.04"), textsize = 5.5, color = "black") +
  scale_x_discrete(limits = c("dmso", "inca"), labels = c("DMSO", "INCA-6")) +
  scale_fill_manual(name = "Bacteria", labels = labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") +
  guides(color = "none") + theme_minimal() +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom")

ggsave(paste("THP-1_VEGFA_qPCR_agg_all", "_", today, ".png", sep = ""), thp1.plot, width = 5.15, height = 6.2, dpi = 300)
  
sem <- ddply(thp1.dagg, c("drug", "bacteria"), summarize, mean = mean(rq), sem = sd(rq)/sqrt(length(rq)))
sem <- transform(sem, lower = mean-sem, upper = mean+sem)
colnames(sem) <- c("drug", "bacteria", "rq", "sem", "lower", "upper")

thp1.col <- ggplot(sem, aes(x = drug, y = rq, fill = bacteria)) + 
  geom_col(position = position_dodge(1)) + 
  geom_errorbar(data = sem, aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(1), size = 1.5, color = "black") + 
  xlab("Treatment") + ylab("RQ VEGFA") + # scale_y_continuous(trans='log10') +
  scale_x_discrete(limits = c("dmso", "inca"), labels = c("DMSO", "INCA-6")) +
  # c("-INCA-6 \n" mmtb, "+INCA-6, \n -Mtb", "-INCA-6, \n +Mtb", "+INCA-6, \n +Mtb")
  scale_fill_manual(name = "Bacteria", labels = labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") +
  geom_signif(y_position = c(12, 15, 8), xmin = c(0.8, 1.2, 1.8), xmax = c(1.2, 2.2, 2.2), annotations = c("p < 0.001", "p < 0.001", "p = 0.04"), textsize = 5.5, color = "black") +
  # guides(color = "none") + 
  # dark_mode() + 
  theme_minimal() +
  theme(text = element_text(size = 20, face = "bold")) + theme(legend.position="bottom") #, plot.background = element_rect(fill = "grey15"), 
# legend.background = element_rect(fill = "grey15"), panel.background = element_rect(fill = "grey10"))

ggsave(paste("THP-1_VEGFA_qPCR_agg_bar", "_", today, ".png", sep = ""), thp1.col, width = 5.15, height = 6.2, dpi = 300)