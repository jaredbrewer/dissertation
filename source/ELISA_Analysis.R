# Jared Brewer
# Created: 25 September 2021
# Last Edited: 10 November 2022
# ELISA Analysis

library(drc)
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

stdcrv <- read.csv("./THP1_ELISA_StdCrv_092721.csv", header = T)

model1 <- drm(OD~Conc,
              fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
              data=stdcrv)
plot(model1)

##

exp3 <- read.csv("./THP1_ELISA_Rep3-0924_092721.csv", header = T)

exp3 <- aggregate(exp3$od,
                  by = list(rep = exp3$dup, treat = exp3$treat, cond = exp3$cond),
                  FUN = function(x) c(mean = mean(x)))

DOSEx <- ED(model1, exp3$x, type = "relative", display = F)
exp3.data <- as.data.frame(DOSEx)
exp3.data$adjust <- exp3.data$Estimate*3

exp3$conc <- exp3.data$adjust
exp3[is.na(exp3)] <- 0
names(exp3)[names(exp3) == "x"] <- "raw"


exp3.aov <- aov(conc ~ treat*cond, data = exp3)

##

exp2 <- read.csv("./THP1_ELISA_Rep2-0923_092721.csv", header = T)

DOSEx.2 <- ED(model1, exp2$od, type = "relative", display = F)
exp2.data <- as.data.frame(DOSEx.2)
exp2.data$adjust <- exp2.data$Estimate*3

exp2$conc <- exp2.data$adjust
exp2[is.na(exp2)] <- 0

exp2.aov <- aov(conc ~ treat*cond, data = exp2)

##

exp1 <- read.csv("./R_VEGF_ELISA_gMtb_INCA-6_8+24hr_091721.csv", header = T)
exp1.aov <- aov(blanked_od ~ treat*cond, data = exp1)

labs <- c(expression(paste(bold("Control"))), 
          expression(paste(bolditalic("γMtb"))))
 
#

elisa.plot <- ggplot(exp3, aes(x = cond, y = conc, color, fill = treat)) + scale_y_continuous(limits = c(0,22)) +
  geom_quasirandom(aes(shape = cond, size = 5), dodge.width = 0.9, show.legend = F) +
  scale_shape_manual(values = c(21, 22)) +
  geom_boxplot(aes(fill = treat), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Exposure Condition") + ylab("[VEGFA] pg/mL") + 
  geom_signif(y_position = c(15, 17, 19), xmin = c(0.8, 0.8, 1.8), xmax = c(1.8, 2.2, 2.2), annotations = c("p < 0.0001", "p = 0.101", "p < 0.0001"), textsize = 8, color = "black", vjust = -0.35) +
  scale_x_discrete(limits = c("control", "mtb"), labels = c("Control", "γMtb")) +
  scale_fill_manual(name = "Treatment", labels = c("DMSO", "40μM INCA-6 "), values = c("firebrick3", "deepskyblue")) +
  theme_minimal() +
  theme(text = element_text(size = 20, face = "bold")) +
  theme(legend.position="bottom")

fn <- paste("exp3_ELISA_", today, ".png", sep = "")
ggsave(fn, elisa.plot, width = 5, height = 7, units = "in", dpi = 300)

#

elisa.plot <- ggplot(exp2, aes(x = cond, y = conc, color, fill = treat)) + scale_y_continuous(limits = c(0,40)) +
  geom_quasirandom(aes(shape = cond, size = 5), dodge.width = 0.9, show.legend = F) +
  scale_shape_manual(values = c(21, 22)) +
  geom_boxplot(aes(fill = treat), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Exposure Condition") + ylab("[VEGFA] pg/mL") + 
  geom_signif(y_position = c(30, 34, 38), xmin = c(0.8, 0.8, 1.8), xmax = c(1.8, 2.2, 2.2), annotations = c("p < 0.0001", "p = 0.193", "p = 0.0001"), textsize = 8, color = "black", vjust = -0.35) +
  scale_x_discrete(limits = c("control", "mtb"), labels = c("Control", "γMtb")) +
  scale_fill_manual(name = "Treatment", labels = c("DMSO", "40μM INCA-6 "), values = c("firebrick3", "deepskyblue")) +
  theme_minimal() +
  theme(text = element_text(size = 20, face = "bold")) +
  theme(legend.position="bottom")

fn <- paste("exp2_ELISA_", today, ".png", sep = "")
ggsave(fn, elisa.plot, width = 5, height = 7, units = "in", dpi = 300)

#

elisa.plot <- ggplot(exp1, aes(x = cond, y = blanked_od, color, fill = treat)) + scale_y_continuous(limits = c(0,1.5)) +
  geom_quasirandom(aes(shape = cond, size = 5), dodge.width = 0.9, show.legend = F) +
  scale_shape_manual(values = c(21, 22)) +
  geom_boxplot(aes(fill = treat), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Exposure Condition") + ylab("Blanked OD") + 
  geom_signif(y_position = c(0.75, 1, 1.25), xmin = c(0.8, 0.8, 1.8), xmax = c(1.8, 2.2, 2.2), annotations = c("p < 0.0001", "p = 0.205", "p < 0.0001"), textsize = 8, color = "black", vjust = -0.35) +
  scale_x_discrete(limits = c("control", "mtb"), labels = c("Control", "γMtb")) +
  scale_fill_manual(name = "Treatment", labels = c("DMSO", "40μM INCA-6 "), values = c("firebrick3", "deepskyblue")) +
  theme_minimal() +
  theme(text = element_text(size = 20, face = "bold")) +
  theme(legend.position="bottom")

fn <- paste("exp1_ELISA_", today, ".png", sep = "")
ggsave(fn, elisa.plot, width = 5, height = 7, units = "in", dpi = 300)

