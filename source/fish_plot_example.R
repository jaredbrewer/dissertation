library(ggplot2)
library(ggbeeswarm)
library(scales)
library(extrafont)
library(reshape)
library(plyr)
library(ggsignif)
library(RColorBrewer)
library(FSA)
library(gghighlight)
library(mdthemes)

setwd("~/Documents")
today <- format(Sys.time(), "%m%d%y")

ttl <- expression(paste(bolditalic("nfatc2a"^"xt69")))

nfatc2a <- read.csv("./Nfatc2a_Infection/JB292/nfatc2a/JB292_quant.csv", header = T)
nfatc2a$exp <- "JB292"
nfatc2a <- nfatc2a[, c("Label", "Length")]
key <- read.csv("./Nfatc2a_Infection/JB292/JB292_genotypes.csv", header = T)
key$rep <- paste(key$rep, ".czi", sep = "")

nfatc2a.vasc <- aggregate(Length ~ Label, nfatc2a, sum)
colnames(nfatc2a.vasc) <- c("rep", "vasc")

merged <- merge(nfatc2a.vasc, key, by="rep", all=TRUE)
nfatc2a.merged <- na.omit(merged)

results <- dunnTest(nfatc2a.merged$vasc ~ as.factor(nfatc2a.merged$gen), two.sided = F, method = "bh")
results

nfatc2a.plot <- ggplot(nfatc2a.merged, aes(x = gen, y = vasc, color, fill = gen)) + # scale_y_continuous(limits = c(0,1100)) +
  geom_beeswarm(aes(shape = gen, size = 5), dodge.width = 0.9) +
  scale_shape_manual(values = c(21, 22, 24)) +
  geom_boxplot(aes(fill = gen), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Genotype") + ylab("Length of Abnormal Vasculature (μm)") + 
  geom_signif(y_position = c(450, 525, 600), xmin = c(1, 1, 2), xmax = c(2, 3, 3), annotations = c("p < 0.001", "p < 0.001", "p < 0.001"), textsize = 5.5, color = "black") +
  scale_x_discrete(limits = c("ctrl", "het", "mut"), labels = c("Wild-type", "Heterozygous", "Mutant")) +
  scale_fill_manual(name = "Genotype", labels = c("ctrl", "het", "mut"), values = c("firebrick3", "springgreen3", "deepskyblue")) + theme(legend.position = "none") + 
  guides(color = "none") + theme_minimal() +
  ggtitle(ttl) +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

fn <- paste("JB292_nfatc2a_", today, ".png", sep = "")
ggsave(fn, nfatc2a.plot, width = 6.625, height = 5.875, units = "in", dpi = 300)