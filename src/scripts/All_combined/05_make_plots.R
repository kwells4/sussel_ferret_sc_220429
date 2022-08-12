library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)

source("src/scripts/functions.R")

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "All_combined"

wt_samples <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14")

cfko_samples <- c("CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")

sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

# Add psuedotime to object
psuper_obj_wt <- readRDS(file.path(sample_dir, "rda_obj", "wt_psupertime.rds"))

psuper_obj_cfko <- readRDS(file.path(sample_dir, "rda_obj",
                                   "cfko_psupertime.rds"))

project_info_wt <- psuper_obj_wt$proj_dt %>%
  tibble::column_to_rownames("cell_id")

project_info_cfko <- psuper_obj_cfko$proj_dt %>%
  tibble::column_to_rownames("cell_id")

project_info <- rbind(project_info_cfko, project_info_wt)

merged_seurat <- AddMetaData(merged_seurat, metadata = project_info)

all_samples <- c(wt_samples, cfko_samples)

cluster_plots <- lapply(all_samples, function(x){
  plotDimRed(merged_seurat, col_by = "RNA_cluster", plot_type = "rna.umap",
             highlight_group = TRUE, group = x, meta_data_col = "sample")[[1]]
})


names(cluster_plots) <- all_samples

pdf(file.path(sample_dir, "images", "cluster_plots_separate.pdf"),
    width = 8, height = 8)
print(cluster_plots)
dev.off()

pdf(file.path(sample_dir, "images", "cluster_plots_combined.pdf"),
    width = 15, height = 12)

print(cowplot::plot_grid(cluster_plots$WT_D2, cluster_plots$WT_D5,
                   cluster_plots$WT_D7, cluster_plots$CFKO_D2,
                   cluster_plots$CFKO_D5, cluster_plots$CFKO_D7,
                   cluster_plots$WT_D9, cluster_plots$WT_D14,
                   NULL, cluster_plots$CFKO_D9, cluster_plots$CFKO_D14,
                   NULL, nrow = 4, ncol = 3))

dev.off()


pdf(file.path(sample_dir, "images", "pseudotime_plots.pdf"),
    width = 8, height = 8)

psupertime_plot_wt <- plotDimRed(merged_seurat, col_by = "psuper",
                                 plot_type = "rna.umap",
                                 highlight_group = TRUE, group = "WT",
                                 meta_data_col = "genotype")


psupertime_plot_cfko <- plotDimRed(merged_seurat, col_by = "psuper",
                                 plot_type = "rna.umap",
                                 highlight_group = TRUE, group = "CFKO",
                                 meta_data_col = "genotype")

print(psupertime_plot_wt)
print(psupertime_plot_cfko)

dev.off()


saveRDS(merged_seurat, file.path(sample_dir, "rda_obj/seurat_processed.rds"))
