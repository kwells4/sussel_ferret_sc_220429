library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)

source("src/scripts/functions.R")

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "All_combined"

sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

de_dir <- file.path(sample_dir, "files/DE")

samples <- c("D2", "D5", "D7", "D9", "D14")

gene_path <- here("files/GSEA_signaling_pathways_with_orthologs.xlsx")

all_sheets <- openxlsx::getSheetNames(gene_path)

gene_lists <- lapply(all_sheets, function(x){
  gene_df <- openxlsx::readWorkbook(gene_path, sheet = x)
  return(unique(gene_df$gene_id))
})

names(gene_lists) <- all_sheets

list_of_interest <- gene_lists$`EP genes`


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
## Sample DE -------------------------------------------------------------------

de_tests <- c("D2_all", "D5_all", "D7_all", "D9_all", "D14_all")

# Make plots for all tests 
invisible(lapply(de_tests, function(de_test){
  print(de_test)
  
  excel_file <- file.path(de_dir, paste0(de_test, ".xlsx"))
  
  excel_sheets <- openxlsx::getSheetNames(excel_file)
  
  excel_sheets <- excel_sheets[!(grepl("_gse", excel_sheets))]
  
  de_genes <- lapply(excel_sheets, function(x){
    de_df <- openxlsx::readWorkbook(excel_file, sheet = x)
    de_df$up_in <- x
    return(de_df)
  })
  
  names(de_genes) <- excel_sheets
  
  de_genes <- do.call(rbind, de_genes)
  
  de_genes <- unique(de_genes$gene)
  
  de_genes <- intersect(de_genes, list_of_interest)
  
  
  fig_height <- round(length(de_genes) / 5)
  
  # Plot heatmap across all samples
  pdf(file.path(sample_dir, "images", "sample_ep_heatmaps",
                paste0(de_test, "_heatmap_average.pdf")))
  
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes),
                     meta_col = "sample", average_expression = TRUE,
                     colors = sample_colors, plot_rownames = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "sample_ep_heatmaps",
                paste0(de_test, "_heatmap_average_labeled.pdf")),
      width = 8, height = fig_height)
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes),
                     meta_col = "sample", average_expression = TRUE,
                     colors = sample_colors, plot_rownames = TRUE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "sample_ep_heatmaps",
                paste0(de_test, "_heatmap_cells.pdf")))
  
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes),
                     meta_col = "sample", average_expression = FALSE,
                     colors = sample_colors, plot_rownames = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "sample_ep_heatmaps",
                paste0(de_test, "_heatmap_cells_labeled.pdf")),
      width = 8, height = fig_height)
  
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes),
                     meta_col = "sample", average_expression = FALSE,
                     colors = sample_colors, plot_rownames = TRUE))
  
  dev.off()
  
}))

## Celltype DE -----------------------------------------------------------------
de_tests <- c("D2_combined_celltype", "D5_combined_celltype",
              "D7_combined_celltype", "D9_combined_celltype",
              "D14_combined_celltype")


merged_seurat$sample_celltype <- paste(merged_seurat$sample,
                                       merged_seurat$RNA_combined_celltype,
                                       sep = "_")

# Make plots for all tests 
invisible(lapply(de_tests, function(de_test){
  print(de_test)
  
  excel_file <- file.path(de_dir, paste0(de_test, ".xlsx"))
  
  excel_sheets <- openxlsx::getSheetNames(excel_file)
  
  excel_sheets <- excel_sheets[!(grepl("_gse", excel_sheets))]
  
  de_genes <- lapply(excel_sheets, function(x){
    de_df <- openxlsx::readWorkbook(excel_file, sheet = x)
    de_df$up_in <- x
    return(de_df)
  })
  
  names(de_genes) <- excel_sheets
  
  de_genes <- do.call(rbind, de_genes)
  
  de_genes_intersect <- intersect(unique(de_genes$gene), list_of_interest)
  
  de_genes <- de_genes[de_genes$gene %in% de_genes_intersect,]
  
  fig_height <- round(length(de_genes_intersect) / 10)
  
  # Keep only the one time point
  subset_seurat <- subset(merged_seurat,
                          subset = day == gsub("_.*", "", de_test))
  
  # Make meta data and colors
  meta_df <- subset_seurat[[c("sample", "RNA_combined_celltype",
                              "sample_celltype")]]
  meta_ave <- data.frame(table(paste0(subset_seurat$sample, ";",
                                      subset_seurat$RNA_combined_celltype)))
  
  meta_ave$RNA_combined_celltype <- gsub(".*_D[0-9]+;", "", meta_ave$Var1)
  meta_ave$sample <- gsub(";.*$", "", meta_ave$Var1)
  meta_ave$Var1 <- NULL
  meta_ave$Freq <- NULL
  meta_ave$sample_celltype <- paste(meta_ave$sample,
                                    meta_ave$RNA_combined_celltype,
                                    sep = "_")
  
  meta_ave$sample_celltype <- factor(meta_ave$sample_celltype)
  
  rownames(meta_ave) <- meta_ave$sample_celltype
  color_list <- list("sample" = sample_colors,
                     "RNA_combined_celltype" = all_colors)
  
  
  # Plot heatmap across all cell types
  pdf(file.path(sample_dir, "images", "cell_type_ep_heatmaps",
                paste0(de_test, "_heatmap_average.pdf")))
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample_celltype", average_expression = TRUE,
                     colors = sample_colors, plot_rownames = FALSE,
                     meta_df = meta_ave, color_list = color_list,
                     plot_meta_col = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "cell_type_ep_heatmaps",
                paste0(de_test, "_heatmap_average_labeled.pdf")),
      width = 8, height = fig_height)
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample_celltype", average_expression = TRUE,
                     colors = sample_colors, plot_rownames = TRUE,
                     meta_df = meta_ave, color_list = color_list,
                     plot_meta_col = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "cell_type_ep_heatmaps",
                paste0(de_test, "_heatmap_cells.pdf")))
  
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample_celltype", average_expression = FALSE,
                     colors = sample_colors, plot_rownames = FALSE,
                     meta_df = meta_df, color_list = color_list,
                     plot_meta_col = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "cell_type_ep_heatmaps",
                paste0(de_test, "_heatmap_cells_labeled.pdf")),
      width = 8, height = fig_height)
  
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample_celltype", average_expression = FALSE,
                     colors = sample_colors, plot_rownames = TRUE,
                     meta_df = meta_df, color_list = color_list,
                     plot_meta_col = FALSE))
  
  dev.off()
  
  
  # Plot heatmaps across one cell type all samples
  de_genes$celltype <- gsub(".*_D[0-9]+_", "", de_genes$up_in)
  
  invisible(lapply(unique(de_genes$celltype), function(y){
    celltype_seurat <- subset(merged_seurat,
                              subset = RNA_combined_celltype == y)
    
    celltype_seurat$sample <- droplevels(celltype_seurat$sample)
    
    new_de_genes <- de_genes %>%
      dplyr::filter(celltype == y)
    
    new_de_genes <- unique(new_de_genes$gene)
    
    new_de_genes <- intersect(new_de_genes, list_of_interest)
    
    
    fig_height <- round(length(new_de_genes) / 10)
    
    # Plot heatmap across all samples
    pdf(file.path(sample_dir, "images", "cell_type_ep_heatmaps",
                  paste0(de_test, "_", y,
                         "_heatmap_average.pdf")))
    print(plot_heatmap(celltype_seurat, gene_list = unique(new_de_genes),
                       meta_col = "sample", average_expression = TRUE,
                       colors = sample_colors, plot_rownames = FALSE))
    
    dev.off()
    
    pdf(file.path(sample_dir, "images", "cell_type_ep_heatmaps",
                  paste0(de_test, "_", y, "_heatmap_average_labeled.pdf")),
        width = 8, height = fig_height)
    print(plot_heatmap(celltype_seurat, gene_list = unique(new_de_genes),
                       meta_col = "sample", average_expression = TRUE,
                       colors = sample_colors, plot_rownames = TRUE))
    
    dev.off()
    
    pdf(file.path(sample_dir, "images", "cell_type_ep_heatmaps",
                  paste0(de_test, "_", y,
                         "_heatmap_cells.pdf")))
    
    print(plot_heatmap(celltype_seurat, gene_list = unique(new_de_genes),
                       meta_col = "sample", average_expression = FALSE,
                       colors = sample_colors, plot_rownames = FALSE))
    
    dev.off()
    
    pdf(file.path(sample_dir, "images", "cell_type_ep_heatmaps",
                  paste0(de_test, "_", y, "_heatmap_cells_labeled.pdf")),
        width = 8, height = fig_height)
    
    print(plot_heatmap(celltype_seurat, gene_list = unique(new_de_genes),
                       meta_col = "sample", average_expression = FALSE,
                       colors = sample_colors, plot_rownames = TRUE))
    
    dev.off()
    
  }))
  
  
}))

