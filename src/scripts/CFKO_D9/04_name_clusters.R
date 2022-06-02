library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(here)
library(scAnalysisR)
library(viridis)
library(clustifyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "CFKO_D9"

source(here("src/scripts/functions.R"))

all_ref_dir <-
  "/Users/wellskr/Documents/Analysis/references/single_cell_references"

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

#-------------------------------------------------------------------------------

########################
# Pancreas development #
########################

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pancreas/byrnes_2018_mouse")

ref_mat <- read.csv(file.path(ref_dir, "E14_epithelial_average.csv"),
                    header = TRUE, row.names = 1)


cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                   save_dir = save_dir,
                                   save_name = "celltype_byrnes",
                                   mapping_file = mapping_file,
                                   mapping_gene_col = "gene_id",
                                   mapping_ortholog_col = "Mouse.gene.name",
                                   assay = "RNA",
                                   nfeatures = 1500, clusters = "RNA_cluster",
                                   plot_type = "rna.umap")

seurat_data <- cluster_res$object

seurat_res_byrnes <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_celltype_byrnes",
           plot_type = "rna.umap")

write.csv(seurat_res_byrnes, file.path(save_dir,
                                "files/celltype_mapping_byrnes_2018.csv"))


#-------------------------------------------------------------------------------

##################
# Pancreas atlas #
##################

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pancreas/Baron_2016")

ref_mat <- read.csv(file.path(ref_dir, "clustifyr_mouse_reference.csv"),
                    header = TRUE, row.names = 1)

cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                   save_dir = save_dir,
                                   save_name = "baron_celltype",
                                   mapping_file = mapping_file,
                                   mapping_gene_col = "gene_id",
                                   mapping_ortholog_col = "Mouse.gene.name",
                                   assay = "RNA",
                                   nfeatures = 1500, clusters = "RNA_cluster",
                                   plot_type = "rna.umap")


seurat_data <- cluster_res$object

seurat_res_baron <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_baron_celltype", plot_type = "rna.umap")

write.csv(seurat_res_baron,
          file.path(save_dir, "files/celltype_mapping_baron_2016.csv"))

#-------------------------------------------------------------------------------

################
# Tabula muris #
################

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pancreas/tabula_muris")

ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference.csv"),
                    header = TRUE, row.names = 1)

cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                   save_dir = save_dir,
                                   save_name = "tabula_muris_celltype",
                                   mapping_file = mapping_file,
                                   mapping_gene_col = "gene_id",
                                   mapping_ortholog_col = "Mouse.gene.name",
                                   assay = "RNA",
                                   nfeatures = 2000, clusters = "RNA_cluster",
                                   plot_type = "rna.umap",
                                   cor_cutoff = 0.4)


seurat_data <- cluster_res$object

seurat_res_tabula_muris <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_tabula_muris_celltype",
           plot_type = "rna.umap")

write.csv(seurat_res_tabula_muris,
          file.path(save_dir, "files/celltype_mapping_tabula_muris.csv"))


#-------------------------------------------------------------------------------

############
# Combined #
############

full_res <- cbind(cbind(seurat_res_baron, seurat_res_byrnes),
                  seurat_res_tabula_muris)

seurat_cluster <- cor_to_call(full_res) %>% 
  mutate(type = ifelse(r < 0.4, "undetermined", type))

new_clusters <- seurat_cluster$type

names(new_clusters) <- seurat_cluster$cluster

seurat_data$RNA_combned_celltype <- new_clusters[seurat_data$RNA_cluster]

plotDimRed(seurat_data, col_by = "RNA_combned_celltype",
           plot_type = "rna.umap")

#-------------------------------------------------------------------------------

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))