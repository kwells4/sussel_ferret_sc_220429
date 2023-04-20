# Save values from heatmaps made in 06_make_plots.R and
# EP_heatmaps.R

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

# DE heatmaps ------------------------------------------------------------------

# First some building

merged_seurat$sample <- factor(merged_seurat$sample,
                               levels = names(sample_colors))

# Read in DE genes
de_directory <- file.path(sample_dir, "files", "DE")

## Sample DE -------------------------------------------------------------------

de_tests <- c("D2_all", "D5_all", "D7_all", "D9_all", "D14_all")

# Make plots for all tests 
all_results <- lapply(de_tests, function(de_test){
  print(de_test)
  
  excel_file <- file.path(de_directory, paste0(de_test, ".xlsx"))
  
  excel_sheets <- openxlsx::getSheetNames(excel_file)
  
  excel_sheets <- excel_sheets[!grepl("gse", excel_sheets)]
  
  de_genes <- lapply(excel_sheets, function(x){
    de_df <- openxlsx::readWorkbook(excel_file, sheet = x)
    de_df$up_in <- x
    return(de_df)
  })
  
  names(de_genes) <- excel_sheets
  
  de_genes <- do.call(rbind, de_genes)
  
  
  res1 <- plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                       meta_col = "sample", average_expression = TRUE,
                       colors = sample_colors, plot_rownames = FALSE,
                       return_data = TRUE)
  
  return(res1)
})

names(all_results) <- de_tests

all_results <- lapply(all_results, function(x){
  x$counts <- log1p(x$counts)
  return(x)
})

sample_de <- openxlsx::createWorkbook()

invisible(lapply(names(all_results), function(sample_info){
  save_data <- all_results[[sample_info]]
  # Add worksheets
  openxlsx::addWorksheet(wb = sample_de,
                         sheetName = paste0(sample_info, "_zscore"))
  
  openxlsx::addWorksheet(wb = sample_de,
                         sheetName = paste0(sample_info, "_avg_counts"))
  
  openxlsx::writeData(wb = sample_de, 
                      sheet = paste0(sample_info, "_zscore"),
                      x = save_data$z_score,
                      rowNames = TRUE)
  
  openxlsx::writeData(wb = sample_de, 
                      sheet = paste0(sample_info, "_avg_counts"),
                      x = save_data$counts,
                      rowNames = TRUE)
  
}))

openxlsx::saveWorkbook(wb = sample_de, 
                       file = file.path(sample_dir, "files",
                                        "z_score_counts_sample_de.xlsx"),
                       overwrite = TRUE)

## EP Celltype DE --------------------------------------------------------------
# EP genes in pancreatic progenitors and CACs in D9 of WT and CF differentiation

de_tests <- c("D9_combined_celltype")

plot_celltypes <- c("centroacinar", "progenitor_like_cells")

list_of_interest <- gene_lists$`EP genes`

merged_seurat$sample_celltype <- paste(merged_seurat$sample,
                                       merged_seurat$RNA_combined_celltype,
                                       sep = "_")

# Make plots for all tests 
final_res <- lapply(de_tests, function(de_test){
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
  

  # Plot heatmaps across one cell type all samples
  de_genes$celltype <- gsub(".*_D[0-9]+_", "", de_genes$up_in)
  
  celltype_plots <- unique(de_genes$celltype)[unique(de_genes$celltype
                                                     %in% plot_celltypes)]
  
  all_results <- lapply(celltype_plots, function(y){
    celltype_seurat <- subset(merged_seurat,
                              subset = RNA_combined_celltype == y)
    
    celltype_seurat$sample <- droplevels(celltype_seurat$sample)
    
    new_de_genes <- de_genes %>%
      dplyr::filter(celltype == y)
    
    new_de_genes <- unique(new_de_genes$gene)
    
    new_de_genes <- intersect(new_de_genes, list_of_interest)
    
    # Plot heatmap across all samples
    
    res1 <- plot_heatmap(celltype_seurat, gene_list = unique(new_de_genes),
                         meta_col = "sample", average_expression = TRUE,
                         colors = sample_colors, plot_rownames = FALSE,
                         return_data = TRUE)
    
    return(res1)
    
  })
  
  names(all_results) <- celltype_plots
  
  return(all_results)
  
})

names(final_res) <- de_tests

final_res <- lapply(final_res$D9_combined_celltype, function(x){
  x$counts <- log1p(x$counts)
  return(x)
})

names(final_res) <- c("CAC", "Pancreatic_prog")

ep_de <- openxlsx::createWorkbook()

invisible(lapply(names(final_res),
                 function(sample_info){
  save_data <- final_res[[sample_info]]
  # Add worksheets
  openxlsx::addWorksheet(wb = ep_de,
                         sheetName = paste0(sample_info, "_zscore"))
  
  openxlsx::addWorksheet(wb = ep_de,
                         sheetName = paste0(sample_info, "_avg_counts"))
  
  openxlsx::writeData(wb = ep_de, 
                      sheet = paste0(sample_info, "_zscore"),
                      x = save_data$z_score,
                      rowNames = TRUE)
  
  openxlsx::writeData(wb = ep_de, 
                      sheet = paste0(sample_info, "_avg_counts"),
                      x = save_data$counts,
                      rowNames = TRUE)
  
}))

openxlsx::saveWorkbook(wb = ep_de, 
                       file = file.path(sample_dir, "files",
                                        "z_score_counts_ep_de.xlsx"),
                       overwrite = TRUE)
