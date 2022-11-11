#I ran this to test the differences in cell types when I use clusters from the 
# combined (with a very high resoultion to simulate the cluster size from the
# independent experiments). For nearly all cases, the same cell types are
# identified when basing on the new clusters (and not rerunning clustifyr).
# We could use these plots to show that the cell types don't differ significantly
# when using all samples combined.


library(slingshot)
library(scAnalysisR)
library(Seurat)
library(here)
library(tidyverse)
library(tradeSeq)
library(pheatmap)
library(clustifyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

all_samples <- "All_combined"

all_sample_dir <- here("results", all_samples, "R_analysis")

save_dir <- all_sample_dir

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap")

plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "harmony.umap")


plotDimRed(merged_seurat, "RNA_combined_celltype_ref_name", plot_type = "harmony.umap")
plotDimRed(merged_seurat, "RNA_combined_celltype_ref_name", plot_type = "rna.umap")


# Corrected cluster vs uncorrected cluster
# RNA_combined cell type vs cluster and uncorrected cluster
matrix_1 <- confusionMatrix(merged_seurat$corrected_cluster,
                            merged_seurat$RNA_combined_celltype)

matrix_2 <- confusionMatrix(merged_seurat$uncorrected_cluster,
                            merged_seurat$RNA_combined_celltype)

matrix_1 <- matrix_1 / rowSums(matrix_1)

matrix_2 <- matrix_2 / rowSums(matrix_2)

pheatmap::pheatmap(matrix_1)

pheatmap::pheatmap(matrix_2)


plotDimRed(merged_seurat, "RNA_cluster", plot_type = "rna.umap")
plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap")

RNA_pcs <- 32
seurat_assay <- "RNA"
HTO <- FALSE

set.seed(0)
umap_data <- group_cells(merged_seurat, nPCs = RNA_pcs,
                         resolution = 10, assay = seurat_assay, HTO = HTO)

seurat_data <- umap_data$object

gene_plots <- umap_data$plots


# RNA_combined cell type vs cluster and uncorrected cluster
matrix_3 <- confusionMatrix(seurat_data$RNA_cluster,
                            seurat_data$RNA_combined_celltype)


matrix_3 <- matrix_3 / rowSums(matrix_3)

pheatmap::pheatmap(matrix_3)


new_clust <- colnames(matrix_3)[apply(matrix_3, 1 , which.max)]

clust_names <- new_clust

names(clust_names) <- rownames(matrix_3)

seurat_data$new_clusters <- clust_names[as.character(seurat_data$RNA_cluster)]

matrix_4 <- confusionMatrix(seurat_data$new_clusters,
                            seurat_data$RNA_combined_celltype)


matrix_4 <- matrix_4 / rowSums(matrix_4)

pheatmap::pheatmap(matrix_4)

plotDimRed(seurat_data, "new_clusters", plot_type = "rna.umap")
plotDimRed(seurat_data, "RNA_combined_celltype", plot_type = "rna.umap")

# Check clustifyr --------------------------------------------------------------

seurat_data$original_celltype <- seurat_data$RNA_combined_celltype

all_ref_dir <-
  "/Users/wellskr/Documents/Analysis/references/single_cell_references"

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

source(here("src/scripts/functions.R"))

cor_cutoff <- 0.4
krentz_cor_cutoff <- 0.3

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

plotDimRed(seurat_data, col_by = "RNA_muraro_celltype",
           plot_type = "rna.umap")

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
                                   nfeatures = 2500,
                                   clusters = "RNA_cluster",
                                   plot_type = "rna.umap",
                                   cor_cutoff = krentz_cor_cutoff)


seurat_data <- cluster_res$object

seurat_res_krentz <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_krentz_celltype",
           plot_type = "rna.umap")

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


# These don't agree super well, but it seems like the disagreements are
# most often in cases where the correlation is close.
plotDimRed(seurat_data, col_by = "RNA_combined_celltype",
           plot_type = "rna.umap")

plotDimRed(seurat_data, col_by = "original_celltype",
           plot_type = "rna.umap")

cm <- confusionMatrix(seurat_data$original_celltype,
      seurat_data$RNA_combined_celltype)


#-------------------------------------------------------------------------------
############################
# Combined with ref winner #
############################

colnames(seurat_res_baron) <- paste0(colnames(seurat_res_baron), "_baron")

colnames(seurat_res_baron_human) <- paste0(colnames(seurat_res_baron_human),
                                           "_baron")
colnames(seurat_res_byrnes) <- paste0(colnames(seurat_res_byrnes),
                                      "_m_byrnes")

colnames(seurat_res_tabula_muris) <- paste0(colnames(seurat_res_tabula_muris),
                                            "_m_tabula_muris")

colnames(seurat_res_qadir) <- paste0(colnames(seurat_res_qadir),
                                     "_h_qadir")

colnames(seurat_res_muraro) <- paste0(colnames(seurat_res_muraro),
                                      "_h_muraro")

colnames(seurat_res_krentz) <- paste0(colnames(seurat_res_krentz),
                                      "_h_krentz")

full_res <- cbind(cbind(seurat_res_baron, seurat_res_byrnes),
                  seurat_res_tabula_muris, seurat_res_baron_human,
                  seurat_res_qadir, seurat_res_muraro,
                  seurat_res_krentz)

seurat_cluster <- cor_to_call(full_res) %>% 
  mutate(type = ifelse(r < cor_cutoff, "undetermined", type))

new_clusters <- seurat_cluster$type

names(new_clusters) <- seurat_cluster$cluster

seurat_data$RNA_combined_celltype_ref_name <- new_clusters[seurat_data$RNA_cluster]

plotDimRed(seurat_data, col_by = "RNA_combined_celltype_ref_name",
           plot_type = "rna.umap")


write.csv(full_res, file = file.path(save_dir, "files", "all_clustifyr_cors.csv"))


plotDimRed(seurat_data, col_by = "RNA_cluster",
           plot_type = "rna.umap",
           highlight_group = TRUE,
           group = 10,
           meta_data_col = "RNA_cluster")

plotDimRed(seurat_data, col_by = "sample",
           plot_type = "rna.umap")


# Cluster 10 is now centroacinar, was ductal
new_res <- full_res["10",]

new_res[order(new_res)]

#-------------------------------------------------------------------------------

# Check some markers
ductal_markers <- c("CLDN1", "TFF1", "KRT19", "SLC3A1", "MUC1")

ductal_plots <- plotDimRed(seurat_data, col_by = ductal_markers,
           plot_type = "rna.umap")


ductal_genes <- c("ALDH1A3", "CFTR", "AQP1",
                  "DEFB1", "KRT19", "SPP1", "TSPAN8")

ductal_violins <- featDistPlot(seurat_data, geneset = ductal_genes,
                               combine = FALSE, sep_by = "RNA_combined_celltype")



ductal_violins <- featDistPlot(seurat_data, geneset = ductal_genes,
                               combine = FALSE, sep_by = "original_celltype")


acinar_genes <- c("ALB", "CELA3A", "PNLIP",
                  "CEL", "CTRB2", "CPA2")

acinar_violins <- featDistPlot(seurat_data, geneset = acinar_genes,
                               combine = FALSE, sep_by = "original_celltype")
