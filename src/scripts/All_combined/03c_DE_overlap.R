library(openxlsx)
library(tidyverse)
library(here)

sample <- "All_combined"

sample_dir <- here("results", sample, "R_analysis")

save_dir <- file.path(sample_dir, "files")

de_dirs <- c("DE", "DE_mast_re")

samples <- c("D2", "D5", "D7", "D9", "D14")

# Check overlaps

get_de_genes <- function(excel_file){
  sheetname <- openxlsx::getSheetNames(excel_file)
  sheetname <- sheetname[!grepl("gse", sheetname)]
  gene_list <- lapply(sheetname, function(sheet){
    results <- openxlsx::readWorkbook(xlsxFile = excel_file, sheet = sheet)
    return(unique(results$gene))  
  })
  names(gene_list) <- sheetname
  return(gene_list)
}

get_overlaps <- function(samples, de_dirs, excel_suffix){
  overlapping_vals <- lapply(samples, function(x){
    individual_de_sets <- lapply(de_dirs, function(dir){
      excel_file <- file.path(save_dir, dir, paste0(x, excel_suffix))
      de_list <- get_de_genes(excel_file)
    })
    names(individual_de_sets) <- de_dirs
    
    all_comparisions <- lapply(names(individual_de_sets[[1]]), function(comparison){
      list_one <- individual_de_sets[[1]][[comparison]]
      list_two <- individual_de_sets[[2]][[comparison]]
      intersecting_genes <- intersect(list_one, list_two)
      return_info <- data.frame("comparison" = comparison,
                                "DE_length" = length(list_one),
                                "DE_mast_re_len" = length(list_two),
                                "intersection" = length(intersecting_genes))
      return(return_info)
      
    })
    
    all_comparisions <- do.call(rbind, all_comparisions)
    
    all_comparisions$sample <- x
    
    return(all_comparisions)
  })
  
  overlapping_vals <- do.call(rbind, overlapping_vals)
  
  overlapping_vals$percent_DE <- overlapping_vals$intersection / 
    overlapping_vals$DE_length
  
  overlapping_vals$percent_mast <- overlapping_vals$intersection /
    overlapping_vals$DE_mast_re_len
  
  return(overlapping_vals)
}

all_overlaps <- get_overlaps(samples = samples,
                             de_dirs = de_dirs,
                             excel_suffix = "_all.xlsx")

combined_celltype_overlaps <- get_overlaps(samples = samples,
                                           de_dirs = de_dirs,
                                           excel_suffix = "_combined_celltype.xlsx")
