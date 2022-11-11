library(slingshot)
library(scAnalysisR)
library(Seurat)
library(here)
library(tidyverse)
library(tradeSeq)
library(pheatmap)
library(MetBrewer)

source(here("src/scripts/functions.R"))

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

all_samples <- "All_combined"

all_sample_dir <- here("results", all_samples, "R_analysis")

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))


wtd2 <- subset(merged_seurat, subset = orig.ident == "WT_D2")

wtd2 <- subset(wtd2, subset = RNA_combined_celltype == "ductal")

wtd2$type <- ifelse(wtd2$original_cluster == 4, "ductal2", "ductal1")

Idents(wtd2) <- "type"
markers <- FindMarkers(wtd2, ident.1 = "ductal1", ident.2 = "ductal2")

markers <- markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::arrange(desc(abs(avg_log2FC)))


v_plots <- featDistPlot(wtd2, geneset = c("percent.mt", "nCount_RNA", "nFeature_RNA"),
             sep_by = "type", combine = FALSE)

# TODO run GO analysis on these genes


