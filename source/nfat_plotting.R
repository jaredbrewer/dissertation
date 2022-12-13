# Jared Brewer
# Created: 1 December 2022
# Last Edited: 13 December 2022
# NFAT Plotting

labs <- c(expression(paste(bolditalic("nfatc1"))), 
          expression(paste(bolditalic("nfatc2a"))),
          expression(paste(bolditalic("nfatc2b"))),
          expression(paste(bolditalic("nfatc3a"))),
          expression(paste(bolditalic("nfatc3b"))),
          expression(paste(bolditalic("nfatc4"))))


ggplot(nfat, aes(x = as.factor(source), y = fpkm, color = gene, fill = gene)) + 
  geom_beeswarm(aes(shape = gene, size = 5), dodge.width = 0.9, show.legend = F) + 
  geom_boxplot(aes(fill = gene), alpha = 0.25, outlier.shape = NA, position = position_dodge(width = 0.9), fatten = 5) + 
  xlab("Source") + ylab("FPKM") + 
  scale_fill_manual(name = "NFAT Isoforms", limits = c("nfatc1", "nfatc2a", "nfatc2b", "nfatc3a", "nfatc3b", "nfatc4"), 
                    labels = labs, 
                    values = c("firebrick3", "deepskyblue", "springgreen3", "darkorchid4", "darkorange1", "darkgoldenrod1")) + 
  guides(color = "none") + 
  scale_color_manual(name = "NFAT Isoforms", limits = c("nfatc1", "nfatc2a", "nfatc2b", "nfatc3a", "nfatc3b", "nfatc4"), 
                     labels = labs, 
                     values = c("firebrick3", "deepskyblue", "springgreen3", "darkorchid4", "darkorange1", "darkgoldenrod1")) + 
  guides(color = "none") + 
  scale_x_discrete(limits = c("granuloma", "macrophage"), labels = c("Granulomas", "Kidney Macrophages")) + 
  theme_minimal() +
  theme(text = element_text(size = 20, face = "bold")) +
  theme(legend.position="bottom")

ggsave("nfat_expression.png", width = 11, height = 6, units = "in", dpi = 300)

