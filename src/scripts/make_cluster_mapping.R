library(Seurat)
library(tidyverse)
library(openxlsx)
library(here)

samples <- c("CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14", 
             "WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14")

celltypes <- c("RNA_celltype_byrnes", "RNA_baron_celltype",
               "RNA_tabula_muris_celltype", "RNA_baron_human_celltype",
               "RNA_qadir_celltype", "RNA_muraro_celltype")

base_dir <- here()

cluster_workbook <- openxlsx::createWorkbook()

lapply(samples, function(x){
  base_dir_proj <- file.path(base_dir, "results", x)
  save_dir <- file.path(base_dir_proj, "R_analysis")
  
  # Read in the data
  seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))
  
  all_celltypes <- lapply(celltypes, function(y){
    new_data <- data.frame(table(paste(seurat_data$RNA_cluster,
                                       seurat_data[[y]][[1]], sep = "_")))
    
    new_data$Var1 <- factor(new_data$Var1)
    new_data <- new_data[order(new_data$Var1),]
    new_data$Freq <- NULL
    new_data$cluster <- sub("_.*", "", new_data$Var1)
    new_data$Var1 <- sub("[0-9]*_", "", new_data$Var1)
    colnames(new_data) <- c(y, "cluster")
    return(new_data)
  })
  
  all_celltypes <- do.call(cbind, all_celltypes)
  
  column_order <- c("cluster", celltypes)
  all_celltypes <- all_celltypes[ , column_order]
  
  openxlsx::addWorksheet(wb = cluster_workbook, sheetName = x)
  openxlsx::writeData(wb = cluster_workbook, sheet = x, x = all_celltypes)
  
})

openxlsx::saveWorkbook(wb = cluster_workbook,
                       file = here("results/all_cluster_mapping.xlsx"),
                       overwrite = TRUE)