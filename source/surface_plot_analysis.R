# Jared Brewer
# Created: 21 July 2022
# Last Edited: 21 July 2022
# 3D Surface Profile Plot Analysis

library(ggplot2)
library(ggbeeswarm)
library(ggsignif)
library(dplyr)
library(DescTools)
library(stringr)

# Plug in the path were the .csv files are - make sure they are the only CSVs in the directory -
 # unless you want to get more creative with some regex.

counts.dir <- ""
counts <- list.files(counts.dir, full.names = T)

cells <- data.frame()

# Our objective here is to combine all the individual CSV files into a
 # single file, but to do a series of operations on it first.
# We read in the file, then aggregate into the minimum value
 # (the data is inverted - need min, not max.)
# We then do some basic math operations to normalize the data across the length of the cell.

for (file in counts) {
  if (endsWith(file, ".csv")) {
    cells <- read.csv(file, header = F) |>
      aggregate(V3 ~ V1 + V2, data = _, min) |>
      mutate(pct_x = (V2 - min(V2))/(max(V2) - min(V2)),
             pct_y = (V3 - min(V3))/(max(V3) - min(V3))) |>
      mutate(x_pct_bin = ntile(pct_x, n=100),
             x_pct_corr = pct_x*100,
             y_pct_corr = -pct_y*100+100,
             treat = str_extract(V1, "WTCAFA|WTMAFA"),
             rep = as.numeric(str_extract(V1, "\\d+$"))) |>
      rbind(cells)
  }
}

aucs <- data.frame()

# Now we want to measure the area under the curve for each replicate.
# We use nested for loops to read each treatment and each replicate for each treatment and filter.
# Then, we normalize across the total Y range, calculate the AUCs, and select only the columns we need.

for (t in unique(cells$treat)) {
  for (r in unique(cells$rep)) {
    aucs <- filter(cells, treat == t) |>
      filter(rep == r) |>
      mutate(y_pct_corr = y_pct_corr/sum(y_pct_corr)) |>
      mutate(per = AUC(x_pct_corr, y_pct_corr, from = min(0), to = max(25))/
                        AUC(x_pct_corr, y_pct_corr, from = min(25), to = max(75)),
             mid = AUC(x_pct_corr, y_pct_corr, from = min(37.5), to = max(62.5))) |>
      select(per, mid, treat, rep) |>
      unique() |>
      na.omit() |>
      rbind(aucs)
  }
}

t.test(aucs$per ~ aucs$treat)

# Plot the AUCs - must add your own labels.

auc.plot <- ggplot(aucs, aes(x = treat, y = per)) +
  geom_beeswarm(aes(shape = treat, size = 3)) +
  geom_boxplot(aes(fill = treat), alpha = 0.25, outlier.shape = NA, fatten = 5) +
  stat_boxplot(geom='errorbar', width = 0.25) +
  xlab("Genotype") + ylab("AUC") +
  geom_signif(y_position = c(0.6), xmin = c(1), xmax = c(2), annotations = c("p = 0.002"), textsize = 5.5, color = "black") +
  scale_x_discrete(limits = c("", ""), labels = c("", "")) +
  scale_fill_manual(name = "", labels = c("", ""), values = c("black", "black")) + theme(legend.position = "none") +
  guides(color = "none") + theme_minimal() +
  theme(text = element_text(size = 20, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  theme_prism() +
  theme(legend.position = "none")
