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

sample <- "CFKO_D5"

source(here("src/scripts/functions.R"))

all_ref_dir <-
  "/Users/wellskr/Documents/Analysis/references/single_cell_references"

HTO <- FALSE
ADT <- FALSE
cor_cutoff <- 0.4
krentz_cor_cutoff <- 0.3

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

# Mouse ------------------------------------------------------------------------

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
                                   nfeatures = 2000, clusters = "RNA_cluster",
                                   plot_type = "rna.umap",
                                   cor_cutoff = cor_cutoff)

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
                                   nfeatures = 2000, clusters = "RNA_cluster",
                                   plot_type = "rna.umap",
                                   cor_cutoff = cor_cutoff)


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
                                   cor_cutoff = cor_cutoff)


seurat_data <- cluster_res$object

seurat_res_tabula_muris <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_tabula_muris_celltype",
           plot_type = "rna.umap")

write.csv(seurat_res_tabula_muris,
          file.path(save_dir, "files/celltype_mapping_tabula_muris.csv"))


# Human ------------------------------------------------------------------------


##################
# Pancreas atlas #
##################

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pancreas/Baron_2016")

ref_mat <- read.csv(file.path(ref_dir, "clustifyr_human_reference.csv"),
                    header = TRUE, row.names = 1)

cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                   save_dir = save_dir,
                                   save_name = "baron_human_celltype",
                                   mapping_file = mapping_file,
                                   mapping_gene_col = "gene_id",
                                   mapping_ortholog_col = "Human.gene.name",
                                   assay = "RNA",
                                   nfeatures = 2000, clusters = "RNA_cluster",
                                   plot_type = "rna.umap",
                                   cor_cutoff = cor_cutoff)


seurat_data <- cluster_res$object

seurat_res_baron_human <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_baron_human_celltype",
           plot_type = "rna.umap")

write.csv(seurat_res_baron_human,
          file.path(save_dir, "files/celltype_mapping_baron_human_2016.csv"))


#-------------------------------------------------------------------------------

####################
# Ductal reference #
####################

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pancreas/Qadir_2020_human")

ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference.csv"),
                    header = TRUE, row.names = 1)

cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                   save_dir = save_dir,
                                   save_name = "qadir_celltype",
                                   mapping_file = mapping_file,
                                   mapping_gene_col = "gene_id",
                                   mapping_ortholog_col = "Human.gene.name",
                                   assay = "RNA",
                                   nfeatures = 2000, clusters = "RNA_cluster",
                                   plot_type = "rna.umap",
                                   cor_cutoff = cor_cutoff)


seurat_data <- cluster_res$object

seurat_res_qadir <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_qadir_celltype",
           plot_type = "rna.umap")

write.csv(seurat_res_qadir,
          file.path(save_dir, "files/celltype_mapping_qadir_2020.csv"))



#-------------------------------------------------------------------------------

##################
# Pancreas atlas #
##################

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pancreas/Muraro_2016_human")

ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference.csv"),
                    header = TRUE, row.names = 1)

cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                   save_dir = save_dir,
                                   save_name = "muraro_celltype",
                                   mapping_file = mapping_file,
                                   mapping_gene_col = "gene_id",
                                   mapping_ortholog_col = "Human.gene.name",
                                   assay = "RNA",
                                   nfeatures = 2000, clusters = "RNA_cluster",
                                   plot_type = "rna.umap",
                                   cor_cutoff = cor_cutoff)


seurat_data <- cluster_res$object

seurat_res_muraro <- cluster_res$RNA

write.csv(seurat_res_muraro,
          file.path(save_dir, "files/celltype_mapping_muraro_2016.csv"))


#-------------------------------------------------------------------------------

########################
# Pancreas development #
########################

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pancreas/krentz_2018_human_mouse")

ref_mat <- read.csv(file.path(ref_dir,
                              "S6D1_GFP_clustifyr_reference_celltype.csv"),
                    header = TRUE, row.names = 1)

cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                   save_dir = save_dir,
                                   save_name = "krentz_celltype",
                                   mapping_file = mapping_file,
                                   mapping_gene_col = "gene_id",
                                   mapping_ortholog_col = "Human.gene.name",
                                   assay = "RNA",
                                   nfeatures = 2500, clusters = "RNA_cluster",
                                   plot_type = "rna.umap",
                                   cor_cutoff = krentz_cor_cutoff)


seurat_data <- cluster_res$object

seurat_res_krentz <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_krentz_celltype",
           plot_type = "rna.umap")

write.csv(seurat_res_krentz,
          file.path(save_dir, "files/celltype_mapping_krentz_2018_2020.csv"))



#-------------------------------------------------------------------------------

############
# Combined #
############

colnames(seurat_res_baron) <- paste0(colnames(seurat_res_baron), "_m")

colnames(seurat_res_baron_human) <- paste0(colnames(seurat_res_baron_human), "_h")

full_res <- cbind(cbind(seurat_res_baron, seurat_res_byrnes),
                  seurat_res_tabula_muris, seurat_res_baron_human,
                  seurat_res_qadir, seurat_res_muraro,
                  seurat_res_krentz)

seurat_cluster <- cor_to_call(full_res) %>% 
  mutate(type = ifelse(r < cor_cutoff, "undetermined", type))

new_clusters <- seurat_cluster$type

names(new_clusters) <- seurat_cluster$cluster

seurat_data$RNA_combined_celltype <- new_clusters[seurat_data$RNA_cluster]

plotDimRed(seurat_data, col_by = "RNA_combined_celltype",
           plot_type = "rna.umap")

seurat_data$RNA_combined_celltype <- gsub("_h$", "",
                                          seurat_data$RNA_combined_celltype)
seurat_data$RNA_combined_celltype <- gsub("_m$", "",
                                          seurat_data$RNA_combined_celltype)
seurat_data$RNA_combined_celltype <- gsub("Acinar", "acinar",
                                          seurat_data$RNA_combined_celltype)

seurat_data$RNA_combined_celltype <- gsub("Ductal", "ductal",
                                          seurat_data$RNA_combined_celltype)


plotDimRed(seurat_data, col_by = "RNA_combined_celltype",
           plot_type = "rna.umap")

#-------------------------------------------------------------------------------

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
