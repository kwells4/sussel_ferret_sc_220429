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

cell_types <- "RNA_combined_celltype"
clusters <- "RNA_cluster"
celltype_two <- "RNA_krentz_celltype"

pval <- 0.05
logfc <- 0.5
cell_cutoff <- 20

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))



save_dir <- file.path(sample_dir, "files/DE")

# Seurat object is the starting object.
# Sample 1 is the first sample
# Sample 2 is the second sample
# Split by is if there should be any subsetting by cell type or cluster
# this should be a column from metadata. If it is NULL (the default)
# than DE between the two samples as a whole is found
# sample_col is the name of the metadata column with sample info.
# Assumed to be "orig.ident"

find_sample_de_genes <- function(seurat_object, sample_1, sample_2,
                                 excel_file_path,
                                 split_by = NULL,
                                 sample_col = "orig.ident",
                                 p_value = 0.05,
                                 logfc = 0.5,
                                 cell_cutoff = 20,
                                 mapping_file = NULL,
                                 mapping_gene_col = NULL,
                                 mapping_ortholog_col = NULL){
  
  excel_workbook <- createWorkbook()
  
  Idents(seurat_object) <- sample_col
  
  seurat_subset <- subset(seurat_object, idents =
                           c(sample_1, sample_2))
  
  if(is.null(split_by)){
    # Just do de between the two samples
    all_markers <- FindMarkers(seurat_subset, ident.1 = sample_1,
                               ident.2 = sample_2)
    
    all_markers$gene <- rownames(all_markers)
    
    if(!is.null(mapping_file)){
      # Add in ortholog gene
      new_mapping <- mapping_file %>%
        dplyr::select(all_of(c(mapping_gene_col,
                               mapping_ortholog_col)))
      
      all_markers <- merge(all_markers, new_mapping, by.x = "gene",
                            by.y = mapping_gene_col, all.x = TRUE,
                            all.y = FALSE) %>%
        distinct() %>%
        arrange(p_val_adj)
    }
    
    sample1_markers <- all_markers %>%
      filter(p_val_adj < p_value & avg_log2FC > logfc)
    
    sample2_markers <- all_markers %>%
      filter(p_val_adj < p_value & avg_log2FC < -logfc) %>%
      mutate(avg_log2FC = abs(avg_log2FC))
    
    addWorksheet(wb = excel_workbook, sheetName = sample_1)
    
    writeData(wb = excel_workbook, sheet = sample_1,
              x = sample1_markers)
    
    addWorksheet(wb = excel_workbook, sheetName = sample_2)
    
    writeData(wb = excel_workbook, sheet = sample_2,
              x = sample2_markers)
  } else {
    # take only cases that are in both samples
    cell_info <- table(seurat_subset[[split_by]][[1]],
                       seurat_subset[[sample_col]][[1]])
    
    cell_info <- rowSums(cell_info > cell_cutoff)
    
    cell_info <- cell_info[cell_info > 1]
    
    cluster_names <- names(cell_info)
    
    # Find markers for each cluster
    invisible(lapply(cluster_names, function(y){
      return_markers <- cluster_sample_de(cluster_name = y,
                                          sample_column = sample_col,
                                          cluster_column = split_by,
                                          seurat_object = seurat_subset,
                                          sample_1 = sample_1,
                                          sample_2 = sample_2)
      
      return_markers$gene <- rownames(return_markers)
      
      
      if(!is.null(mapping_file)){
        # Add in ortholog gene
        new_mapping <- mapping_file %>%
          dplyr::select(all_of(c(mapping_gene_col,
                                 mapping_ortholog_col)))
        
        return_markers <- merge(return_markers, new_mapping,
                                by.x = "gene",
                             by.y = mapping_gene_col, all.x = TRUE,
                             all.y = FALSE) %>%
          distinct() %>%
          arrange(p_val_adj)
      }
      
      sample1_markers <- return_markers %>%
        filter(p_val_adj < p_value & avg_log2FC > logfc)
      
      sample2_markers <- return_markers %>%
        filter(p_val_adj < p_value & avg_log2FC < -logfc) %>%
        mutate(avg_log2FC = abs(avg_log2FC))
      
      sheet_name_1 <- paste0(sample_1, "_", y)
      sheet_name_2 <- paste0(sample_2, "_", y)
      
      addWorksheet(wb = excel_workbook, sheetName = sheet_name_1)
      
      writeData(wb = excel_workbook, sheet = sheet_name_1,
                x = sample1_markers)
      
      addWorksheet(wb = excel_workbook, sheetName = sheet_name_2)
      
      writeData(wb = excel_workbook, sheet = sheet_name_2,
                x = sample2_markers)
    }))
  }
  
  saveWorkbook(wb = excel_workbook, file = excel_file_path,
               overwrite = TRUE)
}

cluster_sample_de <- function(cluster_name, cluster_column,
                              sample_column,
                              seurat_object, sample_1, sample_2){
  
  Idents(seurat_object) <- cluster_column
  
  seurat_subset <- subset(seurat_object, idents = cluster_name)
  
  Idents(seurat_subset) <- sample_column
  
  all_markers <- FindMarkers(seurat_subset, ident.1 = sample_1,
                             ident.2 = sample_2)
  
  return(all_markers)
}


test_samples <- list("D2" = c("WT_D2", "CFKO_D2"),
                     "D5" = c("WT_D5", "CFKO_D5"),
                     "D7" = c("WT_D7", "CFKO_D7"),
                     "D9" = c("WT_D9", "CFKO_D9"),
                     "D14" = c("WT_D14", "CFKO_D14"))

invisible(lapply(names(test_samples), function(x){
  print(x)
  
  sample_1 <- test_samples[[x]][1]
  sample_2 <- test_samples[[x]][2]
  print(sample_1)
  print(sample_2)
  
  excel_file_path <- file.path(save_dir,
                               paste0(x, "_all.xlsx"))
  
  find_sample_de_genes(seurat_object = merged_seurat,
                       sample_1 = sample_1,
                       sample_2 = sample_2,
                       excel_file_path = excel_file_path,
                       split_by = NULL,
                       sample_col = "orig.ident",
                       p_value = pval, logfc = logfc,
                       cell_cutoff = cell_cutoff,
                       mapping_file = mapping_file,
                       mapping_gene_col = "gene_id",
                       mapping_ortholog_col = c("Mouse.gene.name",
                                                "Human.gene.name",
                                                "Dog.gene.name",
                                                "Pig.gene.name"))
  
  excel_file_path <- file.path(save_dir,
                               paste0(x, "_combined_celltype.xlsx"))
  
  find_sample_de_genes(seurat_object = merged_seurat,
                       sample_1 = sample_1,
                       sample_2 = sample_2,
                       excel_file_path = excel_file_path,
                       split_by = cell_types,
                       sample_col = "orig.ident",
                       p_value = pval, logfc = logfc,
                       cell_cutoff = cell_cutoff,
                       mapping_file = mapping_file,
                       mapping_gene_col = "gene_id",
                       mapping_ortholog_col = c("Mouse.gene.name",
                                                "Human.gene.name",
                                                "Dog.gene.name",
                                                "Pig.gene.name"))
  
  excel_file_path <- file.path(save_dir,
                               paste0(x, "_krentz_celltype.xlsx"))
  
  find_sample_de_genes(seurat_object = merged_seurat,
                       sample_1 = sample_1,
                       sample_2 = sample_2,
                       excel_file_path = excel_file_path,
                       split_by = celltype_two,
                       sample_col = "orig.ident",
                       p_value = pval, logfc = logfc,
                       cell_cutoff = cell_cutoff,
                       mapping_file = mapping_file,
                       mapping_gene_col = "gene_id",
                       mapping_ortholog_col = c("Mouse.gene.name",
                                                "Human.gene.name",
                                                "Dog.gene.name",
                                                "Pig.gene.name"))
  
  
}))

