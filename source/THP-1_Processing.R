# Jared Brewer
# Created: 19 June 2022
# Last Edited: 15 November 2022
# Cell Culture Analysis Pipeline

library(DataCombine)
library(tidyr)
library(ggplot2)
  
today <- format(Sys.time(), "%m%d%y")

# The key is a master key with all of the included experiments in one file.
key <- read.csv("./THP1_quant/keyfile.csv")
colnames(key) <- c("original", "file.path")

replacer <-  function(x){
  str_replace(x, ".tif", "")
}

key <- key |> mutate_all(funs(replacer))

treat.labs <- c(expression(paste(bold("Control"))), 
                expression(paste(bolditalic("γMtb"))))

# Inhibitor Analysis

counts.dir <- "./THP1_quant/inhibitor"
counts <- list.files(counts.dir, full.names = T)

inhibitor <- data.frame()

for (file in counts) {
  if (endsWith(file, ".csv")) {
    csv <- read.csv(file)
    base <- basename(file)
    fn.col <- mutate(csv, file.path = paste0(base))
    inhibitor <- bind_rows(inhibitor, fn.col)
  }
}

inhibitor = subset(inhibitor, Slice %in% "Total")
inhibitor <- dplyr::select(inhibitor, c("Type.1", "Type.2", "Type.3", "Type.4", "file.path"))
colnames(inhibitor) <- c("num_cells", "vegf_pos", "nfat_nuc", "intersect", "file.path")
inhibitor <- data.frame(lapply(inhibitor, function(x) {gsub(".csv", "", x) }))
inhibitor <- merge(inhibitor, key, by = "file.path", all = F)

inhibitor <- data.frame(lapply(inhibitor, function(x) {gsub("mtb_inca", "mtbinca", x) }))
inhibitor <- separate(inhibitor, original, into = c("slide", "condition", "rep"), sep = "_", remove = T)

# Not sure at what point in this pipeline it got confused and thought there were characters?
inhibitor[, c(2:5)] <- sapply(inhibitor[, c(2:5)], as.numeric)

inhibitor$pct_vegf <- inhibitor$vegf_pos/inhibitor$num_cells
inhibitor$pct_int <- inhibitor$intersect/inhibitor$vegf_pos
inhibitor$pct_nuc <- inhibitor$nfat_nuc/inhibitor$num_cells
inhibitor$vegf_nuc <- inhibitor$nfat_nuc/inhibitor$vegf_pos
inhibitor$vegf_int <- inhibitor$intersect/inhibitor$vegf_pos
inhibitor$int_tot <- inhibitor$intersect/inhibitor$num_cells
inhibitor$nuc_vegf <- inhibitor$intersect/inhibitor$nfat_nuc
inhibitor$vegf_nuc_norm <- (inhibitor$nfat_nuc/inhibitor$vegf_pos)/inhibitor$num_cells

inhibitor[is.na(inhibitor)] <- 0

replaces = data.frame(from = c("ctrl", "inca", "mtb", "mtbinca"), to = c("dmso_ctrl", "inca_ctrl", "dmso_mtb", "inca_mtb"))
inhibitor <- FindReplace(data = inhibitor, Var = "condition", 
                         replaceData = replaces, from = "from", to = "to", exact = T)
inhibitor <- separate(inhibitor, condition, into = c("drug", "bacteria"), sep = "_", remove = T)

main.labs <- c(expression(paste(bold("DMSO"))), 
          expression(paste(bold("INCA-6"))))

inh.aov <- aov(inhibitor$pct_vegf ~ inhibitor$treatment*inhibitor$isoform)
TukeyHSD(inh.aov)

inhibitor.plot <- ggplot(inhibitor, aes(x = drug, y = pct_vegf*100, color = bacteria, fill = bacteria)) +
  scale_y_continuous(limits = c(0,50)) +
  geom_beeswarm(aes(shape = drug, size = 5), dodge.width = 0.9, show.legend = F) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  geom_boxplot(aes(fill = bacteria), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9), fatten = 5) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Drug Treatment") + ylab("% VEGFA+ Cells") +
  geom_signif(y_position = c(40, 45, 30), xmin = c(0.8, 1.2, 1.8), xmax = c(1.2, 2.2, 2.2), annotations = c("p < 0.001", "p < 0.001", "p = 0.885"), textsize = 5.5, color = "black") +
  scale_x_discrete(limits = c("dmso", "inca"), labels = main.labs) +
  scale_fill_manual(name = "Treatment", labels = treat.labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") + 
  guides(color = "none") + theme_minimal() +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom")

fn <- paste("VEGFA_inhibitor_", today, ".png", sep = "")
ggsave(fn, inhibitor.plot, height = 6, width = 5, units = "in", dpi = 300)

inh.aov <- aov(inhibitor$pct_nuc ~ inhibitor$treatment*inhibitor$isoform)
TukeyHSD(inh.aov)

inhibitor.plot <- ggplot(inhibitor, aes(x = drug, y = pct_nuc*100, color = bacteria, fill = bacteria)) +
  scale_y_continuous(limits = c(0,75)) +
  geom_beeswarm(aes(shape = drug, size = 5), dodge.width = 0.9, show.legend = F) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  geom_boxplot(aes(fill = bacteria), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9), fatten = 5) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Drug Treatment") + ylab("% Nuclear NFAT") +
  geom_signif(y_position = c(65, 70, 40), xmin = c(0.8, 1.2, 1.8), xmax = c(1.2, 2.2, 2.2), annotations = c("p < 0.001", "p < 0.001", "p = 0.998"), textsize = 5.5, color = "black") +
  scale_x_discrete(limits = c("dmso", "inca"), labels = main.labs) +
  scale_fill_manual(name = "Treatment", labels = treat.labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") + 
  guides(color = "none") + theme_minimal() +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom")

fn <- paste("NFAT_nuc_inhibitor_", today, ".png", sep = "")
ggsave(fn, inhibitor.plot, height = 6, width = 5, units = "in", dpi = 300)

inh.aov <- aov(inhibitor$int_tot ~ inhibitor$treatment*inhibitor$isoform)
TukeyHSD(inh.aov)

inhibitor.plot <- ggplot(inhibitor, aes(x = drug, y = int_tot*100, color = bacteria, fill = bacteria)) +
  scale_y_continuous(limits = c(0,50)) +
  geom_beeswarm(aes(shape = drug, size = 5), dodge.width = 0.9, show.legend = F) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  geom_boxplot(aes(fill = bacteria), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9), fatten = 5) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Drug Treatment") + ylab("% VEGFA+ and NFAT+") +
  geom_signif(y_position = c(35, 40, 25), xmin = c(0.8, 1.2, 1.8), xmax = c(1.2, 2.2, 2.2), annotations = c("p < 0.001", "p < 0.001", "p = 1.00"), textsize = 5.5, color = "black") +
  scale_x_discrete(limits = c("dmso", "inca"), labels = main.labs) +
  scale_fill_manual(name = "Treatment", labels = treat.labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") + 
  guides(color = "none") + theme_minimal() +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom")

fn <- paste("NFAT_VEGFA_inhibitor_", today, ".png", sep = "")
ggsave(fn, inhibitor.plot, height = 6, width = 5, units = "in", dpi = 300)

# Isoforms Analysis

counts.dir <- "./THP1_quant/isoforms"
counts <- list.files(counts.dir, full.names = T)

isoforms <- data.frame()

for (file in counts) {
  if (endsWith(file, ".csv")) {
    csv <- read.csv(file)
    base <- basename(file)
    fn.col <- mutate(csv, file.path = paste0(base))
    isoforms <- bind_rows(isoforms, fn.col)
  }
}

isoforms = subset(isoforms, Slice %in% "Total")
isoforms <- dplyr::select(isoforms, c("Type.1", "Type.2", "Type.3", "Type.4", "file.path"))
colnames(isoforms) <- c("num_cells", "vegf_pos", "nfat_nuc", "intersect", "file.path")
isoforms <- data.frame(lapply(isoforms, function(x) {gsub(".csv", "", x) }))
isoforms <- merge(isoforms, key, by = "file.path", all = F)
isoforms <- separate(isoforms, original, into = c("treatment", "MAX", "isoform", "rep"), sep = "_", remove = T)
isoforms <- data.frame(lapply(isoforms, function(x) {gsub("uninfected", "control", x) }))

isoforms[, c(2:5)] <- sapply(isoforms[, c(2:5)], as.numeric)

isoforms$pct_vegf <- isoforms$vegf_pos/isoforms$num_cells
isoforms$pct_int <- isoforms$intersect/isoforms$vegf_pos
isoforms$pct_nuc <- isoforms$nfat_nuc/isoforms$num_cells
isoforms$vegf_nuc <- isoforms$nfat_nuc/isoforms$vegf_pos
isoforms$vegf_int <- isoforms$intersect/isoforms$vegf_pos
isoforms$int_tot <- isoforms$intersect/isoforms$num_cells
isoforms$nuc_vegf <- isoforms$intersect/isoforms$nfat_nuc
isoforms$vegf_nuc_norm <- (isoforms$nfat_nuc/isoforms$vegf_pos)/isoforms$num_cells

isoforms[is.na(isoforms)] <- 0

iso.aov <- aov(isoforms$pct_int ~ isoforms$treatment*isoforms$isoform)
TukeyHSD(iso.aov)

isoforms.plot <- ggplot(isoforms, aes(x = isoform, y = pct_int*100, color = treatment, fill = treatment)) +
  scale_y_continuous(limits = c(0,120)) +
  geom_beeswarm(aes(shape = isoform, size = 5), dodge.width = 0.9, show.legend = F) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  geom_boxplot(aes(fill = treatment), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9), fatten = 5) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Isoform") + ylab("% Double Positive Cells / VEGFA+ Cells") + 
  geom_signif(y_position = c(70, 115, 105, 60), xmin = c(0.8, 1.2, 1.8, 2.8), xmax = c(1.2, 2.2, 2.2, 4.2), annotations = c("p < 0.001", "p < 0.001", "p < 0.001", "p > 0.5"), textsize = 5.5, color = "black") +
  scale_x_discrete(limits = c("nfatc1", "nfatc2", "nfatc3", "nfatc4"), labels = c("NFATC1", "NFATC2", "NFATC3", "NFATC4")) +
  scale_fill_manual(name = "Treatment", labels = treat.labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") + 
  guides(color = "none") + theme_minimal() +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom")

fn <- paste("NFAT_isoforms_", today, ".png", sep = "")
ggsave(fn, inhibitor.plot, height = 7, width = 7, units = "in", dpi = 300)

# Lentivirus Analysis

counts.dir <- "./THP1_quant/lentivirus"
counts <- list.files(counts.dir, full.names = T)

lentivirus <- data.frame()

for (file in counts) {
  if (endsWith(file, ".csv")) {
    csv <- read.csv(file)
    base <- basename(file)
    fn.col <- mutate(csv, file.path = paste0(base))
    lentivirus <- rbind(lentivirus, fn.col)
  }
}

lentivirus <- subset(lentivirus, Slice %in% "Total")
lentivirus <- dplyr::select(lentivirus, c("Type.1", "Type.2", "Type.3", "Type.4", "Type.5", "Type.6", "file.path"))
colnames(lentivirus) <- c("num_cells", "vegf_pos", "nfat_nuc", "cas9_pos", "cas9_vegf", "intersect", "file.path")

lentivirus <- data.frame(lapply(lentivirus, function(x) {gsub(".csv", "", x) }))
lentivirus <- merge(lentivirus, key, by = "file.path", all = F)
lentivirus <- separate(lentivirus, original, into = c("slide", "genotype", "treatment", "rep"), sep = "_", remove = T)
lentivirus <- data.frame(lapply(lentivirus, function(x) {gsub("ST", "CT", x) }))

lentivirus[, c(2:7)] <- sapply(lentivirus[, c(2:7)], as.numeric)

lentivirus$pct_vegf <- lentivirus$vegf_pos/lentivirus$num_cells
lentivirus$pct_int <- lentivirus$intersect/lentivirus$vegf_pos
lentivirus$pct_nuc <- lentivirus$nfat_nuc/lentivirus$num_cells
lentivirus$nuc_vegf <- lentivirus$intersect/lentivirus$nfat_nuc
lentivirus$int_tot <- lentivirus$intersect/lentivirus$num_cells
lentivirus$cas9_pct <- lentivirus$cas9_vegf/lentivirus$vegf_pos

lentivirus[is.na(lentivirus)] <- 0

main.labs <- c(expression(paste(bold("Safe Targeting"))), 
          expression(paste(bolditalic("NFATC2"))))

lenti.aov <- aov(lentivirus$pct_vegf ~ lentivirus$treatment*lentivirus$genotype)
TukeyHSD(lenti.aov)

lenti.plot <- ggplot(lentivirus, aes(x = genotype, y = pct_vegf*100, color = treatment, fill = treatment)) +
  # scale_y_continuous(limits = c(0,50)) +
  geom_beeswarm(aes(shape = genotype, size = 5), dodge.width = 0.9, show.legend = F) +
  scale_shape_manual(values = c(21, 22)) +
  geom_boxplot(aes(fill = treatment), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9), fatten = 5) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Genotype") + ylab("% VEGF+ Cells") +
  geom_signif(y_position = c(40, 45), xmin = c(0.8, 1.2), xmax = c(1.2, 2.2), annotations = c("p = 0.002", "p = 0.026"), textsize = 5.5, color = "black") +
  scale_x_discrete(limits = c("CT", "NFATC2"), labels = main.labs) +
  scale_fill_manual(name = "Treatment", labels = treat.labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") + 
  guides(color = "none") + theme_minimal() +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom")

fn <- paste("lenti_VEGF_", today, ".png", sep = "")
ggsave(fn, lenti.plot, width = 6.375, height = 5.875, units = "in", dpi = 300)

lenti.aov <- aov(lentivirus$int_tot ~ lentivirus$treatment*lentivirus$genotype)
TukeyHSD(lenti.aov)

lenti.plot <- ggplot(lentivirus, aes(x = genotype, y = int_tot*100, color = treatment, fill = treatment)) +
  # scale_y_continuous(limits = c(0,50)) +
  geom_beeswarm(aes(shape = genotype, size = 5), dodge.width = 0.9, show.legend = F) +
  scale_shape_manual(values = c(21, 22)) +
  geom_boxplot(aes(fill = treatment), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9), fatten = 5) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Genotype") + ylab("% VEGFA+ and NFAT+ of Total") +
  geom_signif(y_position = c(35, 40), xmin = c(0.8, 1.2), xmax = c(1.2, 2.2), annotations = c("p = 0.002", "p = 0.003"), textsize = 5.5, color = "black") +
  scale_x_discrete(limits = c("CT", "NFATC2"), labels = main.labs) +
  scale_fill_manual(name = "Treatment", labels = treat.labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") + 
  guides(color = "none") + theme_minimal() +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom")

fn <- paste("lenti_int_", today, ".png", sep = "")
ggsave(fn, lenti.plot, width = 6.375, height = 5.875, units = "in", dpi = 300)

lenti.aov <- aov(lentivirus$nuc_vegf ~ lentivirus$treatment*lentivirus$genotype)
TukeyHSD(lenti.aov)

lenti.plot <- ggplot(lentivirus, aes(x = genotype, y = nuc_vegf*100, color = treatment, fill = treatment)) +
  # scale_y_continuous(limits = c(0,50)) +
  geom_beeswarm(aes(shape = genotype, size = 5), dodge.width = 0.9, show.legend = F) +
  scale_shape_manual(values = c(21, 22)) +
  geom_boxplot(aes(fill = treatment), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9), fatten = 5) +
  stat_boxplot(geom='errorbar', width = 0.25, position = position_dodge(width = 0.9)) +
  xlab("Genotype") + ylab("% VEGF+ and NFAT+ \n of NFAT+ Cells") +
  geom_signif(y_position = c(70, 80), xmin = c(0.8, 1.2), xmax = c(1.2, 2.2), annotations = c("p = 0.003", "p = 0.021"), textsize = 5.5, color = "black") +
  scale_x_discrete(limits = c("CT", "NFATC2"), labels = main.labs) +
  scale_fill_manual(name = "Treatment", labels = treat.labs, values = c("firebrick3", "deepskyblue")) + theme(legend.position = "none") + 
  guides(color = "none") + theme_minimal() +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom")

fn <- paste("lenti_VEGFnuc_", today, ".png", sep = "")
ggsave(fn, lenti.plot, width = 6.375, height = 5.875, units = "in", dpi = 300)
