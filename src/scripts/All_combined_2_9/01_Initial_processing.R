library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "All_combined_2_9"

normalization_method <- "log" # can be SCT or log

HTO <- FALSE
ADT <- FALSE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
base_dir <- here()
base_dir_proj <- file.path(base_dir, "results", sample)

save_dir <- file.path(base_dir_proj, "R_analysis")

# Make output directories
ifelse(!dir.exists(file.path(base_dir_proj)),
       dir.create(file.path(base_dir_proj)), FALSE)

ifelse(!dir.exists(file.path(save_dir)),
       dir.create(file.path(save_dir)), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

samples <- c("CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", 
             "WT_D2", "WT_D5", "WT_D7", "WT_D9")


seurat_objects <- lapply(samples, function(x){
  base_dir_proj <- file.path(base_dir, "results", x)
  save_dir <- file.path(base_dir_proj, "R_analysis")
  
  # Read in the data
  seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))
  
  seurat_data$original_cluster <- seurat_data$RNA_cluster
  
  return(seurat_data)
  
})


# Merge Seurat objects
seurat_data <- merge(seurat_objects[[1]], seurat_objects[2:length(seurat_objects)])

seurat_data$sample <- seurat_data$orig.ident

# Remove doublet finding columns
remove_cols <- c(colnames(seurat_data[[]])[grepl("DF.classifications_",
                                                 colnames(seurat_data[[]]))],
                 colnames(seurat_data[[]])[grepl("pANN_",
                                                 colnames(seurat_data[[]]))])

for (i in remove_cols){
  seurat_data[[i]] <- NULL
}


seurat_data <- seurat_data %>%
  FindVariableFeatures() %>%
  ScaleData()

saveRDS(seurat_data, file = file.path(save_dir, "rda_obj",
                                      "seurat_start.rds"))

