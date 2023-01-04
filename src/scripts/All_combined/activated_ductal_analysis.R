library(scAnalysisR)
library(Seurat)
library(here)
library(tidyverse)
library(pheatmap)
library(MetBrewer)
library(openxlsx)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

source("src/scripts/functions.R")


normalization_method <- "log" # can be SCT or log


pval <- 0.05
logfc <- 0.5
cell_cutoff <- 20

all_samples <- "All_combined"

all_sample_dir <- here("results", all_samples, "R_analysis")

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

mapping_file <- read.csv(here("files/species_mapping_file.csv"))


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

celltype_levels <- c("acinar", "ductal", "Prolif_acinar", "Prolif_ductal",
                     "progenitor_like_cells", "centroacinar",
                     "transitional_to_acinar")

pdf(file.path(all_sample_dir, "images", "activated_ductal",
              "celltype_barplot.pdf"))

print(stacked_barplots(merged_seurat, meta_col = "new_celltype",
                       color = new_colors, split_by = "sample") +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)))

dev.off()

# Find conserved markers -------------------------------------------------------
Idents(merged_seurat) <- "new_celltype"
all_markers <- FindConservedMarkers(object = merged_seurat,
                                    ident.1 = "activated_ductal",
                                    grouping.var = "sample")


sig_markers <- all_markers %>%
  dplyr::filter(max_pval < 0.05)

sig_markers$gene <- rownames(sig_markers)

# Find all markers -------------------------------------------------------------

# marker_list <- find_write_markers_orthologs(seurat_object = merged_seurat,
#                                             meta_col = "new_celltype",
#                                             pval = pval,
#                                             logfc = logfc,
#                                             assay = "RNA",
#                                             save_dir = all_sample_dir,
#                                             mapping_file = mapping_file,
#                                             mapping_gene_col = "gene_id",
#                                             mapping_ortholog_col = c("Mouse.gene.name",
#                                                                      "Human.gene.name",
#                                                                      "Dog.gene.name",
#                                                                      "Pig.gene.name"))
# 
# 
# 
# # Logfc 0 ----------------------------------------------------------------------
# 
# activated_ductal_marker <- marker_list %>%
#   dplyr::filter(cluster == "activated_ductal") %>%
#   dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0)
# 
# 
# plot_seurat <- subset(merged_seurat, 
#                       subset = new_celltype %in% c("ductal", "activated_ductal"))
# 
# 
# print(plot_heatmap(plot_seurat, gene_list = activated_ductal_marker$gene,
#                    meta_col = "sample", average_expression = TRUE,
#                    colors = sample_colors, plot_rownames = FALSE,
#                    cluster_rows = TRUE,
#                    cluster_cols = FALSE))
# 
# dev.off()
# 
# plot_seurat$celltype_sample <- paste0(plot_seurat$new_celltype,
#                                       "_", plot_seurat$sample)
# 
# meta_df <- plot_seurat[[c("sample", "new_celltype",
#                             "celltype_sample")]]
# 
# 
# meta_ave <- meta_df
# rownames(meta_ave) <- NULL
# meta_ave <- distinct(meta_ave)
# rownames(meta_ave) <- meta_ave$celltype_sample
# 
# color_list <- list("sample" = sample_colors,
#                    "new_celltype" = new_colors)
# 
# print(plot_heatmap(plot_seurat, gene_list = activated_ductal_marker$gene,
#                    meta_col = "celltype_sample", average_expression = TRUE,
#                    colors = sample_colors, plot_rownames = FALSE,
#                    meta_df = meta_ave, color_list = color_list,
#                    plot_meta_col = FALSE, cluster_rows = TRUE,
#                    cluster_cols = FALSE))
# 
# dev.off()
# 
# meta_ave <- meta_ave %>%
#   arrange(new_celltype)
# 
# meta_ave$celltype_sample <- factor(meta_ave$celltype_sample,
#                                    levels = meta_ave$celltype_sample)
# 
# print(plot_heatmap(plot_seurat, gene_list = activated_ductal_marker$gene,
#                    meta_col = "celltype_sample", average_expression = TRUE,
#                    colors = sample_colors, plot_rownames = FALSE,
#                    meta_df = meta_ave, color_list = color_list,
#                    plot_meta_col = FALSE, cluster_rows = TRUE,
#                    cluster_cols = FALSE))
# 
# dev.off()
# 
# 
# # Logfc 0.5 --------------------------------------------------------------------
# 
# activated_ductal_marker <- marker_list %>%
#   dplyr::filter(cluster == "activated_ductal") %>%
#   dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.5)
# 
# 
# plot_seurat <- subset(merged_seurat, 
#                       subset = new_celltype %in% c("ductal", "activated_ductal"))
# 
# 
# print(plot_heatmap(plot_seurat, gene_list = activated_ductal_marker$gene,
#                    meta_col = "sample", average_expression = TRUE,
#                    colors = sample_colors, plot_rownames = FALSE,
#                    cluster_rows = TRUE,
#                    cluster_cols = FALSE))
# 
# dev.off()
# 
# plot_seurat$celltype_sample <- paste0(plot_seurat$new_celltype,
#                                       "_", plot_seurat$sample)
# 
# meta_df <- plot_seurat[[c("sample", "new_celltype",
#                           "celltype_sample")]]
# 
# 
# meta_ave <- meta_df
# rownames(meta_ave) <- NULL
# meta_ave <- distinct(meta_ave)
# rownames(meta_ave) <- meta_ave$celltype_sample
# 
# color_list <- list("sample" = sample_colors,
#                    "new_celltype" = new_colors)
# 
# print(plot_heatmap(plot_seurat, gene_list = activated_ductal_marker$gene,
#                    meta_col = "celltype_sample", average_expression = TRUE,
#                    colors = sample_colors, plot_rownames = FALSE,
#                    meta_df = meta_ave, color_list = color_list,
#                    plot_meta_col = FALSE, cluster_rows = TRUE,
#                    cluster_cols = FALSE))
# 
# dev.off()
# 
# meta_ave <- meta_ave %>%
#   arrange(new_celltype)
# 
# meta_ave$celltype_sample <- factor(meta_ave$celltype_sample,
#                                    levels = meta_ave$celltype_sample)
# 
# print(plot_heatmap(plot_seurat, gene_list = activated_ductal_marker$gene,
#                    meta_col = "celltype_sample", average_expression = TRUE,
#                    colors = sample_colors, plot_rownames = FALSE,
#                    meta_df = meta_ave, color_list = color_list,
#                    plot_meta_col = FALSE, cluster_rows = TRUE,
#                    cluster_cols = FALSE))
# 
# dev.off()
# 
# # Logfc 1 ----------------------------------------------------------------------
# 
# activated_ductal_marker <- marker_list %>%
#   dplyr::filter(cluster == "activated_ductal") %>%
#   dplyr::filter(p_val_adj < 0.05, avg_log2FC > 1)
# 
# 
# plot_seurat <- subset(merged_seurat, 
#                       subset = new_celltype %in% c("ductal", "activated_ductal"))
# 
# 
# print(plot_heatmap(plot_seurat, gene_list = activated_ductal_marker$gene,
#                    meta_col = "sample", average_expression = TRUE,
#                    colors = sample_colors, plot_rownames = FALSE,
#                    cluster_rows = TRUE,
#                    cluster_cols = FALSE))
# 
# dev.off()
# 
# plot_seurat$celltype_sample <- paste0(plot_seurat$new_celltype,
#                                       "_", plot_seurat$sample)
# 
# meta_df <- plot_seurat[[c("sample", "new_celltype",
#                           "celltype_sample")]]
# 
# 
# meta_ave <- meta_df
# rownames(meta_ave) <- NULL
# meta_ave <- distinct(meta_ave)
# rownames(meta_ave) <- meta_ave$celltype_sample
# 
# color_list <- list("sample" = sample_colors,
#                    "new_celltype" = new_colors)
# 
# print(plot_heatmap(plot_seurat, gene_list = activated_ductal_marker$gene,
#                    meta_col = "celltype_sample", average_expression = TRUE,
#                    colors = sample_colors, plot_rownames = TRUE,
#                    meta_df = meta_ave, color_list = color_list,
#                    plot_meta_col = FALSE, cluster_rows = TRUE,
#                    cluster_cols = FALSE))
# 
# dev.off()
# 
# meta_ave <- meta_ave %>%
#   arrange(new_celltype)
# 
# meta_ave$celltype_sample <- factor(meta_ave$celltype_sample,
#                                    levels = meta_ave$celltype_sample)
# 
# print(plot_heatmap(plot_seurat, gene_list = activated_ductal_marker$gene,
#                    meta_col = "celltype_sample", average_expression = TRUE,
#                    colors = sample_colors, plot_rownames = TRUE,
#                    meta_df = meta_ave, color_list = color_list,
#                    plot_meta_col = FALSE, cluster_rows = TRUE,
#                    cluster_cols = FALSE))
# 
# dev.off()
# 
# 
# conserved markers ------------------------------------------------------------

mean_logfc <- sig_markers %>%
  select(contains("avg_log2FC")) %>%
  dplyr::mutate("logfc_mean" = rowMeans(.), gene = rownames(.)) %>%
  dplyr::filter(logfc_mean > 0)



plot_seurat <- subset(merged_seurat,
                      subset = new_celltype %in% c("ductal", "activated_ductal"))

# 
# print(plot_heatmap(plot_seurat, gene_list = mean_logfc$gene,
#                    meta_col = "sample", average_expression = TRUE,
#                    colors = sample_colors, plot_rownames = FALSE,
#                    cluster_rows = TRUE,
#                    cluster_cols = FALSE))
# 
# dev.off()

plot_seurat$celltype_sample <- paste0(plot_seurat$new_celltype,
                                      "_", plot_seurat$sample)

meta_df <- plot_seurat[[c("sample", "new_celltype",
                          "celltype_sample")]]


meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$celltype_sample

color_list <- list("sample" = sample_colors,
                   "new_celltype" = new_colors)


pdf(file.path(all_sample_dir, "images", "activated_ductal",
              "activated_ductal_de_sample_group.pdf"))

print(plot_heatmap(plot_seurat, gene_list = mean_logfc$gene,
                   meta_col = "celltype_sample", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   meta_df = meta_ave, color_list = color_list,
                   plot_meta_col = FALSE, cluster_rows = TRUE,
                   cluster_cols = FALSE))

dev.off()

meta_ave <- meta_ave %>%
  arrange(new_celltype)

meta_ave$celltype_sample <- factor(meta_ave$celltype_sample,
                                   levels = meta_ave$celltype_sample)

pdf(file.path(all_sample_dir, "images", "activated_ductal",
              "activated_ductal_de_cell_type_group.pdf"))

print(plot_heatmap(plot_seurat, gene_list = mean_logfc$gene,
                   meta_col = "celltype_sample", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   meta_df = meta_ave, color_list = color_list,
                   plot_meta_col = FALSE, cluster_rows = TRUE,
                   cluster_cols = FALSE))

dev.off()

all_logfc <- all_markers %>%
  select(contains("avg_log2FC")) %>%
  dplyr::mutate("logfc_mean" = rowMeans(.), gene = rownames(.)) %>%
  dplyr::filter(logfc_mean > 0)

write_markers <- all_markers %>%
  dplyr::mutate(gene = rownames(.)) %>%
  dplyr::filter(gene %in% all_logfc$gene)

write.csv(write_markers, file.path(all_sample_dir, "images", "activated_ductal",
                               "activated_ductal_de.csv"))

pdf(file.path(all_sample_dir, "images", "activated_ductal",
              "IGFBP7_violin.pdf"))

print(featDistPlot(merged_seurat, geneset = "IGFBP7", combine = FALSE,
                   sep_by = "new_celltype", color = new_colors,
                   col_by = "new_celltype"))

dev.off()

pdf(file.path(all_sample_dir, "images", "activated_ductal",
              "IGFBP3_violin.pdf"))

print(featDistPlot(merged_seurat, geneset = "IGFBP3", combine = FALSE,
                   sep_by = "new_celltype", color = new_colors,
                   col_by = "new_celltype"))

dev.off()


pdf(file.path(all_sample_dir, "images", "activated_ductal",
              "IGFBP7_umap.pdf"))

print(plotDimRed(merged_seurat, col_by = "IGFBP7", plot_type = "rna.umap"))

dev.off()

pdf(file.path(all_sample_dir, "images", "activated_ductal",
              "IGFBP3_umap.pdf"))

print(plotDimRed(merged_seurat, col_by = "IGFBP3", plot_type = "rna.umap"))


dev.off()

