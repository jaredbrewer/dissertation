library(ggplot2)
library(scales)
library(extrafont)
library(reshape)
library(plyr)
library(ggsignif)
library(RColorBrewer)
library(FSA)
library(gghighlight)

setwd("~/Documents")
today <- format(Sys.time(), "%m%d%y")

cfu <- read.csv("./tdTomato_VIVIT_CFU_median_021622.csv", header = T)

t.test(cfu$rel_count~cfu$gen)

labs <- c(expression(paste(bolditalic("irg1:tdTomato"))), 
          expression(paste(bolditalic("irg1:VIVIT"))))

ttl.1 <- expression(atop(paste(bolditalic("irg1:tdTomato "), bold("v.")), paste(bolditalic("irg1:VIVIT "), bold("CFU"))))

cfu.plot <- ggplot(cfu, aes(x = gen, y = rel_count, color = gen, fill = gen)) + scale_y_continuous(limits = c(0,1.25)) +
  geom_dotplot(position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.9), binaxis = "y", binwidth = 0.05, method = "histodot", stackdir = "centerwhole", stackratio = 0.1, dotsize = 1, binpositions = "all") +
  geom_boxplot(aes(fill = gen), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Genotype") + ylab("Relative CFU") + 
  geom_signif(y_position = c(1.1), xmin = c(1), xmax = c(2), annotations = c("p = 0.018"), textsize = 5.5, color = "black") +
  scale_x_discrete(limits = c("tdTomato", "VIVIT"), labels = labs) +
  scale_fill_manual(name = "Genotype", labels = labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") + 
  guides(color = "none") + theme_minimal() +
  labs(title = ttl.1) +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(size = 24, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 24, face = "bold", hjust = 0.5)) +
  theme(legend.position = "none")

cfu.plot

fn <- paste("JB30[234]_cfu_", today, ".png", sep = "")
ggsave(fn, cfu.plot, width = 4, height = 5.875, units = "in", dpi = 300)