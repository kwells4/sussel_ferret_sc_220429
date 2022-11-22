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


merged_seurat$sample <- merged_seurat$orig.ident


# Gene sets
gene_path <- here("files/GSEA_signaling_pathways_with_orthologs.xlsx")

all_sheets <- openxlsx::getSheetNames(gene_path)

gene_lists <- lapply(all_sheets, function(x){
  gene_df <- openxlsx::readWorkbook(gene_path, sheet = x)
  #all_genes <- unique(gene_df$gene_id)
  return(unique(gene_df$gene_id))
})

all_sheets <- sub(" ", "_", all_sheets)

names(gene_lists) <- all_sheets


# Keep only the one time point
subset_seurat <- subset(merged_seurat,
                        subset = day == "D9")

# Remove genes that aren't expressed
expression_counts <- GetAssayData(subset_seurat, slot = "counts") %>%
  rowSums() 

expressed_genes <- expression_counts[expression_counts > 3]

subset_seurat <- subset(subset_seurat, features = names(expressed_genes))




make_heatmaps <- function(seurat_object,
                          plot_list,
                          meta_col = "sample_celltype",
                          name_addition = ""){
  plot_genes <- gene_lists[[plot_list]]
  
  fig_height <- round(length(plot_genes) / 10)
  
  # Make meta data and colors
  meta_df <- seurat_object[[c("sample", "RNA_combined_celltype",
                              meta_col)]]
  
  
  meta_ave <- meta_df
  rownames(meta_ave) <- NULL
  meta_ave <- distinct(meta_ave)
  rownames(meta_ave) <- meta_ave[[meta_col]]
  
  
  color_list <- list("sample" = sample_colors,
                     "RNA_combined_celltype" = all_colors)

  
  # Plot heatmap across all cell types
  pdf(file.path(sample_dir, "images", "gene_list_heatmap",
                paste0(plot_list, "_heatmap_average.pdf")))
  print(plot_heatmap(seurat_object, gene_list = unique(plot_genes),
                     meta_col = meta_col, average_expression = TRUE,
                     colors = sample_colors, plot_rownames = FALSE,
                     meta_df = meta_ave, color_list = color_list,
                     plot_meta_col = FALSE, cluster_rows = TRUE,
                     cluster_cols = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "gene_list_heatmap",
                paste0(plot_list, "_heatmap_average_labeled.pdf")),
      width = 8, height = fig_height)
  print(plot_heatmap(seurat_object, gene_list = unique(plot_genes),
                     meta_col = meta_col, average_expression = TRUE,
                     colors = sample_colors, plot_rownames = TRUE,
                     meta_df = meta_ave, color_list = color_list,
                     plot_meta_col = FALSE,
                     cluster_rows = TRUE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "gene_list_heatmap",
                paste0(plot_list, "_heatmap_cells.pdf")))
  
  print(plot_heatmap(seurat_object, gene_list = unique(plot_genes),
                     meta_col = meta_col, average_expression = FALSE,
                     colors = sample_colors, plot_rownames = FALSE,
                     meta_df = meta_df, color_list = color_list,
                     plot_meta_col = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "gene_list_heatmap",
                paste0(plot_list, "_heatmap_cells_labeled.pdf")),
      width = 8, height = fig_height)
  
  print(plot_heatmap(seurat_object, gene_list = unique(plot_genes),
                     meta_col = meta_col, average_expression = FALSE,
                     colors = sample_colors, plot_rownames = TRUE,
                     meta_df = meta_df, color_list = color_list,
                     plot_meta_col = FALSE))
  
  dev.off()
  
}


subset_seurat$sample_celltype <- factor(subset_seurat$sample_celltype)

for (plot_name in names(gene_lists)){
  print(plot_name)
  
  make_heatmaps(subset_seurat, plot_name)
}

subset_seurat$celltype_sample <- paste(subset_seurat$RNA_combined_celltype,
                                       subset_seurat$sample, sep = "_")

subset_seurat$celltype_sample <- factor(subset_seurat$celltype_sample)

for (plot_name in names(gene_lists)){
  print(plot_name)
  
  make_heatmaps(subset_seurat, plot_name, meta_col = "celltype_sample",
                name_addition = "new_order")
}
