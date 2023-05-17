library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
library(ggridges)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "All_combined"

wt_samples <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14")

cfko_samples <- c("CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")

sample_colors <- as.character(LaCroixColoR::lacroix_palette("Coconut", 10))
sample_colors[5] <- "#F4E3C7"
names(sample_colors) <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                          "CFKO_D14", "CFKO_D9", "CFKO_D7", "CFKO_D5", "CFKO_D2")
sample_colors <- sample_colors[c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                                 "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")]



all_colors <- c("acinar" = "#D4405B",
                "ductal" = "#A5903E",
                "Prolif_acinar" = "#55A470",
                "Prolif_ductal" = "#767FC9",
                "progenitor_like_cells" = "#297878",
                "transitional_to_acinar" = "#874652",
                "centroacinar" = "#78295D")

phase_colors <- LaCroixColoR::lacroix_palette("CranRaspberry", n = 3)
names(phase_colors) <- c("S", "G1", "G2M")

sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

pseudobulk <- AverageExpression(merged_seurat,
                                features = VariableFeatures(merged_seurat),
                                group.by = "RNA_combined_celltype")

pseudobulk <- pseudobulk$RNA


all_corr <- cor(pseudobulk)



pseudobulk <- AverageExpression(merged_seurat,
                                features = VariableFeatures(merged_seurat),
                                group.by = "celltype_sample")

pseudobulk <- pseudobulk$RNA


all_corr <- cor(pseudobulk)
all_corr <- round(all_corr, digits = 2)

pheatmap::pheatmap(all_corr, display_numbers = all_corr)


all_corr <- all_corr[grepl("CFKO", rownames(all_corr)),
                     grepl("CFKO", colnames(all_corr))]

pheatmap::pheatmap(all_corr)

pseudobulk <- AverageExpression(merged_seurat,
                                features = VariableFeatures(merged_seurat),
                                group.by = "celltype_sample")

pseudobulk <- pseudobulk$RNA


all_corr <- cor(pseudobulk)


all_corr <- all_corr[ , grepl("transitional_to", colnames(all_corr))]

pheatmap::pheatmap(all_corr)

all_corr <- round(all_corr, digits = 2)

pheatmap::pheatmap(all_corr, display_numbers = all_corr)
