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


cfko_pseudotime <- slingPseudotime(cfko_slingshot)
colnames(cfko_pseudotime) <- paste0("cfko_", colnames(cfko_pseudotime))

cfko_pseudotime <- data.frame(cfko_pseudotime)

# Check barcodes
dim(cfko_pseudotime)
length(intersect(rownames(cfko_pseudotime), colnames(merged_seurat)))

rownames(cfko_pseudotime) <- rownames(meta_mapper)

merged_seurat <- AddMetaData(merged_seurat, metadata = cfko_pseudotime)

pseudotime_plots <- plotDimRed(merged_seurat, colnames(cfko_pseudotime),
                               plot_type= "rna.umap")



# wt slingshot
wt_slingshot <- readRDS(file.path(all_sample_dir, "files",
                                    "slingshot", "WT_res.rds"))


wt_pseudotime <- slingPseudotime(wt_slingshot)
colnames(wt_pseudotime) <- paste0("wt_", colnames(wt_pseudotime))

wt_pseudotime <- data.frame(wt_pseudotime)

# check barcodes
dim(wt_pseudotime)
length(intersect(rownames(wt_pseudotime), colnames(merged_seurat)))

merged_seurat <- AddMetaData(merged_seurat, metadata = wt_pseudotime)

pseudotime_plots <- plotDimRed(merged_seurat, colnames(wt_pseudotime),
                               plot_type= "rna.umap")
