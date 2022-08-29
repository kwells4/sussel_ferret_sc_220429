library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)
library(openxlsx)
library(gprofiler2)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

source("src/scripts/functions.R")

normalization_method <- "log" # can be SCT or log

sample <- "WT_D9"

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

plotDimRed(seurat_data, col_by = "RNA_cluster", plot_type = "rna.umap")

test_celltype_map <- c("0" = "ductal_1",
                       "1" = "ductal_2",
                       "2" = "ductal_3",
                       "3" = "ductal_1",
                       "4" = "ductal_4",
                       "5" = "ductal_1",
                       "6" = "transitional_to_acinar2",
                       "7" = "cnetroacinar",
                       "8" = "ductal_3",
                       "9" = "ductal_1",
                       "10" = "ductal_1",
                       "11" = "ductal_2",
                       "12" = "progenitor_like_cells",
                       "13" = "ductal_4",
                       "14" = "ductal_2")

seurat_data$test_celltype <- test_celltype_map[seurat_data$RNA_cluster]

plotDimRed(seurat_data, col_by = "test_celltype", plot_type = "rna.umap")

Idents(seurat_data) <- "test_celltype"

# Find markers
seurat_markers <- find_write_markers_orthologs(seurat_object = seurat_data,
                                               meta_col = "test_celltype",
                                               pval = 0.05,
                                               logfc = 0.5,
                                               assay = "RNA",
                                               save_dir = save_dir,
                                               mapping_file = mapping_file,
                                               mapping_gene_col = "gene_id",
                                               mapping_ortholog_col = c("Mouse.gene.name",
                                                                        "Human.gene.name",
                                                                        "Dog.gene.name",
                                                                        "Pig.gene.name"))

# Find overlaps between ductal populations

ductal_genes <- lapply(1:4, function(x){
  subset_val <- paste0("ductal_", x)
  de_genes <- seurat_markers %>%
    dplyr::filter(cluster == subset_val, p_val_adj < 0.05,
                  avg_log2FC > 0.5)
  return(unique(de_genes$gene))
})

names(ductal_genes) <- paste0("ductal_", 1:4)


# Number of intersects
intersects <- lapply(names(ductal_genes), function(x){
  gene_list <- ductal_genes[[x]]
  gene_list_length <- length(gene_list)
  all_intersects <- lapply(names(ductal_genes), function(y){
    second_gene_list <- ductal_genes[[y]]
    overlapping_length <- length(intersect(gene_list, second_gene_list))
    return_df <- data.frame(overlapping_length)
    colnames(return_df) <- y
    return(return_df)
  })
  all_intersects <- do.call(cbind, all_intersects)
  rownames(all_intersects) <- x
  return(all_intersects)
})

intersects <- do.call(rbind, intersects)


# Gene ontology
invisible(lapply(names(ductal_genes), function(x){
  gene_list <- ductal_genes[[x]]
  gene_ontology <- gost(query = gene_list,
                        organism = "mpfuro",
                        ordered_query = TRUE)$result
  
  ontology_res <- createWorkbook()
  invisible(lapply(unique(gene_ontology$source), function(y){
    save_df <- gene_ontology %>%
      dplyr::filter(p_value < 0.05,
                    source == y)
    sheet_name <- gsub(":", "_", y)
    openxlsx::addWorksheet(ontology_res, sheetName = sheet_name)
    openxlsx::writeData(wb = ontology_res, sheet = sheet_name, x = save_df)
  }))
  openxlsx::saveWorkbook(wb = ontology_res, file = file.path(save_dir,
                                                             "files",
                                                             "DE",
                                                             paste0(x,
                                                                    "_GO.xlsx")),
                         overwrite = TRUE)
}))


