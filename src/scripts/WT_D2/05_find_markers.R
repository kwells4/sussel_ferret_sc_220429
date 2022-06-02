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

sample <- "WT_D2"

cell_types <- "RNA_combined_celltype"
clusters <- "RNA_cluster"
celltype_two <- "RNA_tabula_muris_celltype"

pval <- 0.05
logfc <- 0.5

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

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

# Combined ---------------------------------------------------------------------

## Cell type DE ----------------------------------------------------------------

marker_list <- find_write_markers_orthologs(seurat_object = seurat_data,
                                  meta_col = cell_types,
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir,
                                  mapping_file = mapping_file,
                                  mapping_gene_col = "gene_id",
                                  mapping_ortholog_col = c("Mouse.gene.name",
                                                           "Human.gene.name",
                                                           "Dog.gene.name",
                                                           "Pig.gene.name"))

if(ADT){
  marker_list <- find_write_markers_orthologs(seurat_object = seurat_data,
                                    meta_col = cell_types,
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

## RNA cluster DE --------------------------------------------------------------

seurat_data$cluster_celltype <- paste0(seurat_data[[clusters]][[1]], "_",
                                       seurat_data[[cell_types]][[1]])

marker_list <- find_write_markers_orthologs(seurat_object = seurat_data,
                                  meta_col = "cluster_celltype",
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir,
                                  mapping_file = mapping_file,
                                  mapping_gene_col = "gene_id",
                                  mapping_ortholog_col = c("Mouse.gene.name",
                                                           "Human.gene.name",
                                                           "Dog.gene.name",
                                                           "Pig.gene.name"))

if(ADT){
  marker_list <- find_write_markers_orthologs(seurat_object = seurat_data,
                                    meta_col = "cluster_celltype",
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

# Tabula muris -----------------------------------------------------------------

## Celltype DE -----------------------------------------------------------------

# marker_list <- find_write_markers_orthologs(seurat_object = seurat_data,
#                                             meta_col = celltype_two,
#                                             pval = pval,
#                                             logfc = logfc,
#                                             assay = "RNA",
#                                             save_dir = save_dir,
#                                             mapping_file = mapping_file,
#                                             mapping_gene_col = "gene_id",
#                                             mapping_ortholog_col = c("Mouse.gene.name",
#                                                                      "Human.gene.name",
#                                                                      "Dog.gene.name",
#                                                                      "Pig.gene.name"))
# 
# if(ADT){
#   marker_list <- find_write_markers_orthologs(seurat_object = seurat_data,
#                                               meta_col = celltype_two,
#                                               pval = pval,
#                                               logfc = logfc,
#                                               assay = "ADT",
#                                               save_dir = save_dir)
# }

# RNA cluster DE ---------------------------------------------------------------

seurat_data$cluster_celltype_tm <- paste0(seurat_data[[clusters]][[1]], "_",
                                       seurat_data[[celltype_two]][[1]])

marker_list <- find_write_markers_orthologs(seurat_object = seurat_data,
                                            meta_col = "cluster_celltype_tm",
                                            pval = pval,
                                            logfc = logfc,
                                            assay = "RNA",
                                            save_dir = save_dir,
                                            mapping_file = mapping_file,
                                            mapping_gene_col = "gene_id",
                                            mapping_ortholog_col = c("Mouse.gene.name",
                                                                     "Human.gene.name",
                                                                     "Dog.gene.name",
                                                                     "Pig.gene.name"))

if(ADT){
  marker_list <- find_write_markers_orthologs(seurat_object = seurat_data,
                                              meta_col = "cluster_celltype_tm",
                                              pval = pval,
                                              logfc = logfc,
                                              assay = "ADT",
                                              save_dir = save_dir)
}

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
