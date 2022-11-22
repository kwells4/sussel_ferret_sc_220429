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

new_colors <- c(all_colors,
                "activated_ductal" = "#D67229")

celltype_levels <- c("acinar", "activated_ductal", "ductal", "Prolif_acinar",
                     "Prolif_ductal",
                     "progenitor_like_cells", "centroacinar",
                     "transitional_to_acinar")

sample_levels <- names(sample_colors)

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

best_markers <- c(
  "CTCF", # proliferation
  "HMGB2", # proliferation
  "ID3", # ductal
  "IGFBP3", # activate ductal
  "KRT19", # ductal
  "cytokeritins1", # ductal
  "STAT3", # progenitor
  "FOS", # progenitor
  "TSPAN8", # progenitor
  "SOX4", # early progenitor
  "ALDH1A1" # centroacinar
)


# Make a new UMAP --------------------------------------------------------------
all_samples <- names(sample_colors)
all_plots <- lapply(all_samples, function(x){
  plotDimRed(merged_seurat, "new_celltype", plot_type = "rna.umap",
             highlight_group = TRUE, group = x, meta_data_col = "sample",
             color = new_colors)[[1]]
})

names(all_plots) <- all_samples


all_cell_types <- cowplot::plot_grid(all_plots$WT_D2, all_plots$WT_D5,
                                     all_plots$WT_D7, all_plots$WT_D9,
                                     all_plots$WT_D14, NULL, all_plots$CFKO_D2,
                                     all_plots$CFKO_D5, all_plots$CFKO_D7,
                                     all_plots$CFKO_D9, all_plots$CFKO_D14,
                                     NULL, nrow = 4, ncol = 3)


pdf(file.path(all_sample_dir, "images", "combined_cell_type_umap_updated.pdf"),
    width = 15, height = 10)
print(all_cell_types)

dev.off()

# Pseudotime density plots -----------------------------------------------------
# Make density plots for the two lineages colored by cell type
# Get meta data, select only cell type, sample, and lineage,
# keep only values that aren't NA make plot as below:
wt_plotting_df <- merged_seurat[[]] %>%
  dplyr::select(new_celltype, sample, wt_Lineage4) %>%
  dplyr::filter(!is.na(wt_Lineage4))

wt_plotting_df$new_celltype <- factor(wt_plotting_df$new_celltype,
                                               levels = c("activated_ductal",
                                                          "Prolif_ductal",
                                                          "ductal",
                                                          "acinar",
                                                          "transitional_to_acinar",
                                                          "centroacinar",
                                                          "progenitor_like_cells"))

pdf(file.path(all_sample_dir, "images", "wt_pseudotime_ridges_updated.pdf"))
print(ggplot2::ggplot(wt_plotting_df, ggplot2::aes(x = wt_Lineage4,
                                                   y = new_celltype,
                                                   group = new_celltype,
                                                   fill = new_celltype)) + 
        ggridges::geom_density_ridges() +
        ggplot2::scale_fill_manual(values = new_colors) +
        ggplot2::ggtitle("WT lineage 4"))

dev.off()

cfko_plotting_df <- merged_seurat[[]] %>%
  dplyr::select(new_celltype, sample, cfko_Lineage2) %>%
  dplyr::filter(!is.na(cfko_Lineage2))

cfko_plotting_df$new_celltype <- factor(cfko_plotting_df$new_celltype,
                                                 levels = c("transitional_to_acinar",
                                                            "activated_ductal",
                                                            "ductal",
                                                            "acinar",
                                                            "centroacinar",
                                                            "progenitor_like_cells"))


pdf(file.path(all_sample_dir, "images", "cfko_pseudotime_ridges_updated.pdf"))
print(ggplot2::ggplot(cfko_plotting_df, ggplot2::aes(x = cfko_Lineage2,
                                                     y = new_celltype,
                                                     group = new_celltype,
                                                     fill = new_celltype)) + 
        ggridges::geom_density_ridges() +
        ggplot2::scale_fill_manual(values = new_colors) +
        ggplot2::ggtitle("CFKO lineage 2"))

dev.off()


# Violin of genes --------------------------------------------------------------
violin_plots <- featDistPlot(merged_seurat, unique(best_markers),
                              sep_by = "new_celltype",
                              combine = FALSE, color = new_colors)

pdf(file.path(all_sample_dir, "images", "cell_types.pdf"))
print(violin_plots)
dev.off()

merged_seurat$celltype_sample <- paste(merged_seurat$new_celltype,
                                       merged_seurat$sample, sep = "_")

merged_seurat$new_celltype <- factor(merged_seurat$new_celltype,
                                              levels = celltype_levels)

merged_seurat$sample <- factor(merged_seurat$sample, levels = sample_levels)

combined_levels <- merged_seurat[[c("sample", "new_celltype",
                                    "celltype_sample")]] %>%
  distinct() 

combined_levels <- combined_levels[order(combined_levels$new_celltype, 
                                         combined_levels$sample),]



merged_seurat$celltype_sample <- factor(merged_seurat$celltype_sample,
                                        levels = combined_levels$celltype_sample)


pdf(file.path(all_sample_dir, "images",
              "celltype_heatmaps_sample_celltype_updated.pdf"))


meta_df <- merged_seurat[[c("sample", "new_celltype",
                            "celltype_sample")]]


meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$celltype_sample

color_list <- list("sample" = sample_colors,
                   "new_celltype" = new_colors)

print(plot_heatmap(merged_seurat, gene_list = best_markers,
                   meta_col = "celltype_sample", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   meta_df = meta_ave, color_list = color_list,
                   plot_meta_col = FALSE, cluster_rows = TRUE,
                   cluster_cols = FALSE))

dev.off()

pdf(file.path(all_sample_dir, "images",
              "celltype_heatmaps_celltype_updated.pdf"))

meta_df <- merged_seurat[[c("new_celltype")]]


meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$new_celltype

color_list <- list("new_celltype" = new_colors)

print(plot_heatmap(merged_seurat, gene_list = best_markers,
                   meta_col = "new_celltype", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   meta_df = meta_ave, color_list = color_list,
                   plot_meta_col = TRUE, cluster_rows = TRUE,
                   cluster_cols = FALSE))


dev.off()
