# Sleuth Analysis of Kallisto Output
# Created: 29 September 2017
# Last Edited: 5 August 2022
# Jared Brewer

# You have to provide exactly two things:
  # The directory where you did everything and set the working directory to that location using setwd()
  # A sample table formatted with two columns: sample and condition and named sample_table.txt
# After running, you'll get an interactive window where you can play with your data and look for effects.
# There is merit in breaking out comparisons if you have a bunch of conditions...

setwd("~/esxM")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install(c("rhdf5", "biomaRt"))

install.packages("devtools")
devtools::install_github("pachterlab/sleuth")

library(sleuth)

library(tidyverse)
library(biomaRt)
library(data.table)
library(pheatmap)
library(tidyr)

# This analysis was optimized for animals - biomaRt access to bacterial genomes has been discontinued and is a mess.
mart <- read_tsv("./gene_list.txt")
t2g <- as.data.frame(mart$Name)
names(t2g)[1] <- "target_id"
t2g.2 <- t2g[!duplicated(t2g$target_id),]
t2g.2 <- as.data.frame(t2g.2)
names(t2g.2)[1] <- "target_id"
t2g.2$gene_id <- t2g.2$target_id

# Eukaryotes only:
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
t2g <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- rename(t2g, target_id = ensembl_transcript_id_version, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
#

files.kal <- list.files(pattern = "./*_quant")
files.dirs <- sapply(files.kal, function(id) file.path(".", id))
sample.table <- read.table("./sample_table.txt", header = T)
sample.table <- dplyr::mutate(sample.table, path = files.dirs)

so <- sleuth_prep(sample.table,
                  full_model = ~condition,
                  target_mapping = t2g.2,
                  extra_bootstrap_summary = T,
                  aggregation_column = "target_id")
so <- sleuth_fit(so, ~condition, "full")
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
# so <- sleuth_wt(so, which_beta = 'condition', which_model = 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, pval <= 0.05)

sleuth_live(so)

de <- subset(sleuth_table, pval < 0.05)
paths <- enrichPathway(gene = trimws(de$entrez), pvalueCutoff = 0.2, readable = T, organism = "Mycobacterium_tuberculosis_h37rv")

tpm.mat <- sleuth_to_matrix(so, which_df = "obs_raw", which_units = "tpm")
tpm <- as.data.frame(tpm.mat)

# Read straight from Sleuth matrix:
setDT(tpm, keep.rownames = "genes")
esx <- tpm |> filter(endsWith(genes, "_5")) |> column_to_rownames(var = "genes")
esx <- as.matrix(esx)
write.csv(x = tpm.mat, file = "esx_5_tpm.csv")
pheatmap(esx, show_colnames = F, annotation_legend = T, filename = "rnaseq_heatmap_070622.png", width = 6, height = 2)

# Read from saved .csv file (this is the method used for the actual figure):
esx <- read.csv("./esx_5_tpm.csv", header = T)
row.names(esx) <- esx$target_id
esx <- esx |> filter(endsWith(target_id, "_5"))
esx <- esx[-c(1:2)]
esx <- as.matrix(esx)
pheatmap(esx, show_colnames = F, annotation_legend = T, filename = "rnaseq_heatmap_070622.png", width = 6, height = 2)
