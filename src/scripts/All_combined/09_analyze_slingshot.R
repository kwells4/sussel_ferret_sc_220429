library(slingshot)
library(scAnalysisR)
library(Seurat)
library(here)
library(tidyverse)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

all_samples <- "All_combined"

all_sample_dir <- here("results", all_samples, "R_analysis")

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))


# cfko slingshot
cfko_slingshot <- readRDS(file.path(all_sample_dir, "files",
                                   "slingshot", "CFKO_res.rds"))


all_pseudotime <- slingPseudotime(cfko_slingshot)
colnames(all_pseudotime) <- paste0("cfko_", colnames(all_pseudotime))

all_pseudotime <- data.frame(all_pseudotime)

merged_seurat <- AddMetaData(merged_seurat, metadata = all_pseudotime)

pseudotime_plots <- plotDimRed(merged_seurat, colnames(all_pseudotime),
                               plot_type= "rna.umap")



# wt slingshot
wt_slingshot <- readRDS(file.path(all_sample_dir, "files",
                                    "slingshot", "WT_res.rds"))


wt_pseudotime <- slingPseudotime(wt_slingshot)
colnames(wt_pseudotime) <- paste0("wt_", colnames(wt_pseudotime))

wt_pseudotime <- data.frame(wt_pseudotime)

# Rename barcodes to match
meta_mapper <- merged_seurat[[]] %>%
  dplyr::filter(genotype == "WT")

wt_pseudotime <- wt_pseudotime[order(match(rownames(wt_pseudotime),
                                           meta_mapper$WT_barcodes)),]

rownames(wt_pseudotime) <- rownames(meta_mapper)

merged_seurat <- AddMetaData(merged_seurat, metadata = wt_pseudotime)

pseudotime_plots <- plotDimRed(merged_seurat, colnames(wt_pseudotime),
                               plot_type= "rna.umap")


# To do 
# Try the WT with 2 starting populations---> remind yourself how to do this