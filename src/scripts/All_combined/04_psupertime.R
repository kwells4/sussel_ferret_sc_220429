library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
library(psupertime)

source("src/scripts/functions.R")

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "All_combined"

cell_types <- "RNA_combined_celltype"
clusters <- "RNA_cluster"
celltype_two <- "RNA_krentz_celltype"

pval <- 0.05
logfc <- 0.5
cell_cutoff <- 20

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

# Just WT first

wt_samples <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14")
wt_seurat <- subset(merged_seurat, subset = sample %in% wt_samples)

wt_sce <- as.SingleCellExperiment(wt_seurat)

day <- wt_seurat$day

day <- factor(day, levels = c("D2", "D5", "D7", "D9", "D14"))

psuper_obj  <- psupertime(wt_sce, day,
                          sel_genes="hvg")

saveRDS(psuper_obj, file.path(sample_dir, "rda_obj", "wt_psupertime.rds"))


 
# Just CFKO

cfko_samples <- c("CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")
cfko_seurat <- subset(merged_seurat, subset = sample %in% cfko_samples)

ckfo_sce <- as.SingleCellExperiment(cfko_seurat)

day <- cfko_seurat$day

day <- factor(day, levels = c("D2", "D5", "D7", "D9", "D14"))

psuper_obj_cfko  <- psupertime(ckfo_sce, day,
                          sel_genes="hvg")

saveRDS(psuper_obj_cfko, file.path(sample_dir, "rda_obj",
                                   "cfko_psupertime.rds"))

