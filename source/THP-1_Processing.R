# Jared Brewer
# Created: 19 June 2022
# Last Edited: 10 November 2022
# Cell Culture Analysis Pipeline

library(DataCombine)
library(tidyr)
library(ggplot2)
  
today <- format(Sys.time(), "%m%d%y")

# Inhibitor Analysis

counts.dir <- "/Volumes/JB_SD3/THP1_quant/inhibitor/"
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

key <- read.csv("/Volumes/JB_SD3/THP1_quant/keyfile.csv")
colnames(key) <- c("original", "file.path")

inhibitor = subset(inhibitor, Slice %in% "Total")
inhibitor <- select(inhibitor, c("Type.1", "Type.2", "Type.3", "Type.4", "file.path"))
colnames(inhibitor) <- c("num_cells", "vegf_pos", "nfat_nuc", "intersect", "file.path")
replaces = data.frame(from = ".csv", to = ".tif")
inhibitor <- FindReplace(data = inhibitor, Var = "file.path", 
                         replaceData = replaces, from = "from", to = "to", exact = F)
inhibitor <- merge(inhibitor, key, by = "file.path", all = F)

# replaces = data.frame(from = ".csv", to = ".tif")
# inhibitor <- FindReplace(data = inhibitor, Var = "file.path", 
#                          replaceData = replaces, from = "from", to = "to", exact = F)
replaces = data.frame(from = "mtb_inca", to = "mtbinca")
inhibitor <- FindReplace(data = inhibitor, Var = "original", 
                         replaceData = replaces, from = "from", to = "to", exact = F)
inhibitor <- separate(inhibitor, original, into = c("slide", "condition", "rep"), sep = "_", remove = T)

replaces = data.frame(from = ".tif", to = "")
inhibitor <- FindReplace(data = inhibitor, Var = "rep", 
                         replaceData = replaces, from = "from", to = "to", exact = F)

inhibitor$pct_vegf <- inhibitor$vegf_pos/inhibitor$num_cells
inhibitor$pct_int <- inhibitor$intersect/inhibitor$vegf_pos
inhibitor$pct_nuc <- inhibitor$nfat_nuc/inhibitor$num_cells
inhibitor$vegf_nuc <- inhibitor$nfat_nuc/inhibitor$vegf_pos
inhibitor$vegf_int <- inhibitor$intersect/inhibitor$vegf_pos
inhibitor$int_tot <- inhibitor$intersect/inhibitor$num_cells
inhibitor$nuc_vegf <- inhibitor$intersect/inhibitor$nfat_nuc
inhibitor$vegf_nuc_norm <- (inhibitor$nfat_nuc/inhibitor$vegf_pos)/inhibitor$num_cells
inh.aov <- aov(inhibitor$nuc_vegf ~ inhibitor$isoform * inhibitor$inf)
inhibitor[is.na(inhibitor)] <- 0

replaces = data.frame(from = c("ctrl", "inca", "mtb", "mtbinca"), to = c("dmso_ctrl", "inca_ctrl", "dmso_mtb", "inca_mtb"))
inhibitor <- FindReplace(data = inhibitor, Var = "condition", 
                         replaceData = replaces, from = "from", to = "to", exact = T)
inhibitor <- separate(inhibitor, condition, into = c("drug", "bacteria"), sep = "_", remove = T)

inhibitor.2 <- inhibitor[inhibitor$vegf_pos != "1", ]
inhibitor.2[is.na(inhibitor.2)] <- 0

replaces <- data.frame(from = "uninfected", to = "control")
inhibitor <- FindReplace(data = inhibitor, Var = "treatment", 
                        replaceData = replaces, from = "from", to = "to", exact = T)

inh.aov <- aov(inhibitor$pct_int ~ inhibitor$treatment*inhibitor$isoform)
TukeyHSD(inh.aov)

main.labs <- c(expression(paste(bold("DMSO"))), 
          expression(paste(bold("INCA-6"))))

treat.labs <- c(expression(paste(bold("Control"))), 
                expression(paste(bolditalic("γMtb"))))

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
  # ggtitle("") +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom")
  # theme(legend.position = "none")

ggsave("VEGFA_inhibitor_101122.png", plot = inhibitor.plot, dpi = 300, height = 6, width = 5)

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
  # ggtitle("") +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom")
# theme(legend.position = "none")

ggsave("NFAT_nuc_inhibitor_101122.png", plot = inhibitor.plot, dpi = 300, height = 6, width = 5)

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
  # ggtitle("") +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom")

ggsave("NFAT_VEGFA_inhibitor_101122.png", plot = inhibitor.plot, dpi = 300, height = 6, width = 5)

# Isoforms Analysis

counts.dir <- "/Volumes/JB_SD2/THP1_NFAT_031422/Isoforms/Blinded/quant"
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

key <- read.csv("/Volumes/JB_SD2/THP1_NFAT_031422/Isoforms/Blinded/keyfile.csv", header = T)
colnames(key) <- c("original", "file.path")

isoforms = subset(isoforms, Slice %in% "Total")
isoforms <- select(isoforms, c("Type.1", "Type.2", "Type.3", "Type.4", "file.path"))
colnames(isoforms) <- c("num_cells", "vegf_pos", "nfat_nuc", "intersect", "file.path")
replaces = data.frame(from = ".csv", to = ".tif")
isoforms <- FindReplace(data = isoforms, Var = "file.path", 
                        replaceData = replaces, from = "from", to = "to", exact = F)
isoforms <- merge(isoforms, key, by = "file.path", all = F)

# replaces = data.frame(from = ".csv", to = ".tif")
# isoforms <- FindReplace(data = isoforms, Var = "file.path", 
#                          replaceData = replaces, from = "from", to = "to", exact = F)
replaces = data.frame(from = "mtb_inca", to = "mtbinca")
isoforms <- FindReplace(data = isoforms, Var = "original", 
                        replaceData = replaces, from = "from", to = "to", exact = F)
isoforms <- separate(isoforms, original, into = c("treatment", "MAX", "isoform", "rep"), sep = "_", remove = T)

replaces = data.frame(from = ".tif", to = "")
isoforms <- FindReplace(data = isoforms, Var = "rep", 
                        replaceData = replaces, from = "from", to = "to", exact = F)
isoforms$pct_vegf <- isoforms$vegf_pos/isoforms$num_cells
isoforms$pct_int <- isoforms$intersect/isoforms$vegf_pos
isoforms$pct_nuc <- isoforms$nfat_nuc/isoforms$num_cells
isoforms$vegf_nuc <- isoforms$nfat_nuc/isoforms$vegf_pos
isoforms$vegf_int <- isoforms$intersect/isoforms$vegf_pos
isoforms$int_tot <- isoforms$intersect/isoforms$num_cells
isoforms$nuc_vegf <- isoforms$intersect/isoforms$nfat_nuc
isoforms$vegf_nuc_norm <- (isoforms$nfat_nuc/isoforms$vegf_pos)/isoforms$num_cells
inh.aov <- aov(isoforms$nuc_vegf ~ isoforms$isoform * isoforms$inf)
isoforms[is.na(isoforms)] <- 0

replaces = data.frame(from = c("ctrl", "inca", "mtb", "mtbinca"), to = c("dmso_ctrl", "inca_ctrl", "dmso_mtb", "inca_mtb"))
isoforms <- FindReplace(data = isoforms, Var = "treat", 
                        replaceData = replaces, from = "from", to = "to", exact = T)
isoforms <- separate(isoforms, treat, into = c("drug", "bacteria"), sep = "_", remove = T)

isoforms.2 <- isoforms[isoforms$vegf_pos != "1", ]
isoforms.2[is.na(isoforms.2)] <- 0

replaces <- data.frame(from = "uninfected", to = "control")
isoforms <- FindReplace(data = isoforms, Var = "treatment", 
                        replaceData = replaces, from = "from", to = "to", exact = T)

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


# Lentivirus Analysis

counts.dir <- "/Volumes/JB_SD3/THP1_quant/lentivirus"
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

key <- read.csv("/Volumes/JB_SD3/THP1_quant/keyfile.csv")
colnames(key) <- c("original", "file.path")

lentivirus <- subset(lentivirus, Slice %in% "Total")
lentivirus <- select(lentivirus, c("Type.1", "Type.2", "Type.3", "Type.4", "Type.5", "Type.6", "file.path"))
colnames(lentivirus) <- c("num_cells", "vegf_pos", "nfat_nuc", "cas9_pos", "cas9_vegf", "intersect", "file.path")

replaces = data.frame(from = ".csv", to = ".tif")
lentivirus <- FindReplace(data = lentivirus, Var = "file.path", 
                         replaceData = replaces, from = "from", to = "to", exact = F)

lentivirus <- merge(lentivirus, key, by = "file.path", all = F)

replaces = data.frame(from = ".csv", to = ".tif")
lentivirus <- FindReplace(data = lentivirus, Var = "file.path", 
                         replaceData = replaces, from = "from", to = "to", exact = F)

replaces = data.frame(from = "mtb_inca", to = "mtbinca")
lentivirus <- FindReplace(data = lentivirus, Var = "original", 
                         replaceData = replaces, from = "from", to = "to", exact = F)

lentivirus <- separate(lentivirus, original, into = c("slide", "genotype", "treatment", "rep"), sep = "_", remove = T)

replaces = data.frame(from = "ST", to = "CT")
lentivirus <- FindReplace(data = lentivirus, Var = "genotype", 
                          replaceData = replaces, from = "from", to = "to", exact = F)

replaces = data.frame(from = ".tif", to = "")
lentivirus <- FindReplace(data = lentivirus, Var = "rep", 
                         replaceData = replaces, from = "from", to = "to", exact = F)

lentivirus$pct_vegf <- lentivirus$vegf_pos/lentivirus$num_cells
lentivirus$pct_int <- lentivirus$intersect/lentivirus$vegf_pos
lentivirus$pct_nuc <- lentivirus$nfat_nuc/lentivirus$num_cells
lentivirus$nuc_vegf <- lentivirus$intersect/lentivirus$nfat_nuc
lentivirus$int_tot <- lentivirus$intersect/lentivirus$num_cells
lentivirus$cas9_pct <- lentivirus$cas9_vegf/lentivirus$vegf_pos
lenti.aov <- aov(lentivirus$pct_int ~ lentivirus$treatment*lentivirus$genotype)
lentivirus[is.na(lentivirus)] <- 0

replaces = data.frame(from = "ST", to = "CT")
lentivirus <- FindReplace(data = lentivirus, Var = "genotype", 
                          replaceData = replaces, from = "from", to = "to", exact = T)


main.labs <- c(expression(paste(bold("Safe Targeting"))), 
          expression(paste(bolditalic("NFATC2"))))

treat.labs <- c(expression(paste(bold("Control"))), 
                expression(paste(bolditalic("γMtb"))))

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
  # ggtitle("Lentivirus THP-1 \n Immunofluorescence") +
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
  # ggtitle("Lentivirus THP-1 \n Immunofluorescence") +
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
  # ggtitle("Lentivirus THP-1 \n Immunofluorescence") +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom")

fn <- paste("lenti_VEGFnuc_", today, ".png", sep = "")
ggsave(fn, lenti.plot, width = 6.375, height = 5.875, units = "in", dpi = 300)
