# Document information
# This document processes WT and CFKO separately to identify joint celltypes
# for all samples within one genotype. Because they don't cluster together
# I kept the genotypes separate so that genes driven by batch effect wouldn't
# contribute to the cell type identification. The two are then combined
# and some metrics are explored to ensure that the cell types are similar
# between this method and naming celltypes in each sample individually.
# It makes three important metadata columns
# new_celltype = This cell type is identified by taking the clustering from
#   each genotype and finding the cell type from the individual sample analysis
#   that the majority of cells in the cluster belong to. This is mostly used
#   to check for similarity between the analysis.
# final_ind_celltype = This is the cell type that is determined from running
#   clustifyr on all cells for one genotype together. The transitional to 
#   acinar and progenitor cells are then combined.
# ind_full_celltype = This is the same as the final_ind_celltype but the 
#   transitional to acinar and progenitor cells are left separate.

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

sample <- "All_combined_2_9"

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


all_colors2 <- c("acinar" = "#D4405B",
                 "ductal" = "#A5903E",
                 "Prolif_acinar" = "#55A470",
                 "Prolif_ductal" = "#767FC9",
                 "centroacinar_progenitor" = "#297878",
                 "centroacinar" = "#78295D")

sample_colors <- as.character(LaCroixColoR::lacroix_palette("Coconut", 10))
sample_colors[5] <- "#F4E3C7"
names(sample_colors) <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                          "CFKO_D14", "CFKO_D9", "CFKO_D7", "CFKO_D5", "CFKO_D2")
sample_colors <- sample_colors[c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                                 "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")]

sample_colors <- sample_colors[!(grepl("D14", names(sample_colors)))]

# Function -------------------------------------------------------------------

run_all_clustifyr <- function(seurat_data, save_dir,
                              all_var_features, mapping_file,
                              cor_cutoff){
  
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
                                     features = all_var_features,
                                     clusters = "RNA_cluster",
                                     plot_type = "rna.umap",
                                     cor_cutoff = cor_cutoff)
  
  seurat_data <- cluster_res$object
  
  seurat_res_byrnes <- cluster_res$RNA
  

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
                                     features = all_var_features,
                                     clusters = "RNA_cluster",
                                     plot_type = "rna.umap",
                                     cor_cutoff = cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_baron <- cluster_res$RNA
  

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
                                     features = all_var_features, 
                                     clusters = "RNA_cluster",
                                     plot_type = "rna.umap",
                                     cor_cutoff = cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_tabula_muris <- cluster_res$RNA
  

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
                                     features = all_var_features, 
                                     clusters = "RNA_cluster",
                                     plot_type = "rna.umap",
                                     cor_cutoff = cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_baron_human <- cluster_res$RNA
  

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
                                     features = all_var_features,
                                     clusters = "RNA_cluster",
                                     plot_type = "rna.umap",
                                     cor_cutoff = cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_qadir <- cluster_res$RNA

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
                                     features = all_var_features,
                                     clusters = "RNA_cluster",
                                     plot_type = "rna.umap",
                                     cor_cutoff = cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_muraro <- cluster_res$RNA

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
                                     features = all_var_features,
                                     clusters = "RNA_cluster",
                                     plot_type = "rna.umap",
                                     cor_cutoff = krentz_cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_krentz <- cluster_res$RNA
  

  
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

  seurat_data$RNA_combined_celltype <- gsub("_h$", "",
                                            seurat_data$RNA_combined_celltype)
  seurat_data$RNA_combined_celltype <- gsub("_m$", "",
                                            seurat_data$RNA_combined_celltype)
  seurat_data$RNA_combined_celltype <- gsub("Acinar", "acinar",
                                            seurat_data$RNA_combined_celltype)
  
  seurat_data$RNA_combined_celltype <- gsub("Ductal", "ductal",
                                            seurat_data$RNA_combined_celltype)
  

  return(seurat_data)

}


# Analysis ------------------------------------------------------------------

## WT alone ----------------------------------------------------------------
seurat_wt <- subset(seurat_data, subset = genotype == "WT")

RNA_pcs <- 32

seurat_wt <- seurat_wt %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

seurat_wt <- PCA_dimRed(seurat_wt, assay = seurat_assay)

umap_data <- group_cells(seurat_wt, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 3, assay = "RNA", HTO = FALSE)

seurat_wt <- umap_data$object

cM <- confusionMatrix(seurat_wt$RNA_cluster,
                      seurat_wt$final_celltype)

pheatmap::pheatmap(cM)
pheatmap::pheatmap(cM / rowSums(cM))


new_clust <- colnames(cM)[apply(cM, 1 , which.max)]

clust_names <- new_clust

names(clust_names) <- rownames(cM)

seurat_wt$new_celltype <- clust_names[as.character(seurat_wt$RNA_cluster)]


cm <- confusionMatrix(seurat_wt$RNA_cluster, 
                      seurat_wt$final_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

all_var_features <- lapply(unique(seurat_wt$sample), function(sample_name){
  subset_seurat <- subset(seurat_wt, subset = sample == sample_name)
  subset_seurat <- FindVariableFeatures(subset_seurat, nfeatures = 1500)
  return(VariableFeatures(subset_seurat))
})

all_var_features <- unique(unlist(all_var_features))

all_var_features <- VariableFeatures(seurat_wt)

seurat_wt <- run_all_clustifyr(seurat_data = seurat_wt, 
                               save_dir = save_dir,
                               all_var_features = all_var_features,
                               mapping_file = mapping_file,
                               cor_cutoff = cor_cutoff)

### Rename cell types -------------------------------------------------------
seurat_name_mapping <- unique(seurat_wt$RNA_combined_celltype)
names(seurat_name_mapping) <- seurat_name_mapping

seurat_name_mapping[seurat_name_mapping == "progenitor_like_cells"] = 
  "centroacinar_progenitor"

seurat_name_mapping[seurat_name_mapping == "transitional_to_acinar1"] = 
  "centroacinar_progenitor"

seurat_name_mapping[seurat_name_mapping == "transitional_to_acinar2"] = 
  "centroacinar_progenitor"

seurat_wt$final_wt_celltype <- seurat_name_mapping[as.character(seurat_wt$RNA_combined_celltype)]

seurat_wt$final_wt_celltype <- factor(seurat_wt$final_wt_celltype,
                                     levels = names(all_colors2))


### Check -------------------------------------------------------------------

cm <- confusionMatrix(seurat_wt$final_wt_celltype, 
                      seurat_wt$final_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)


cm <- confusionMatrix(seurat_wt$final_celltype, 
                      seurat_wt$final_wt_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)


cm <- confusionMatrix(seurat_wt$new_celltype, 
                      seurat_wt$final_wt_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

all_colors3 <- c("acinar" = "#D4405B",
                 "ductal" = "#A5903E",
                 "Prolif_acinar" = "#55A470",
                 "Prolif_ductal" = "#767FC9",
                 "centroacinar2" = "#297878",
                 "centroacinar1" = "#78295D")
plotDimRed(seurat_wt, "final_celltype", plot_type = "rna.umap",
           color = all_colors3)

plotDimRed(seurat_wt, "final_wt_celltype", plot_type = "rna.umap",
           color = all_colors2)

wt_meta <- seurat_wt[[]] %>%
  dplyr::select(new_celltype, final_wt_celltype, RNA_combined_celltype)

## CFKO alone ----------------------------------------------------------------
seurat_cfko <- subset(seurat_data, subset = genotype == "CFKO")

RNA_pcs <- 32

seurat_cfko <- seurat_cfko %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

seurat_cfko <- PCA_dimRed(seurat_cfko, assay = seurat_assay)

umap_data <- group_cells(seurat_cfko, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 3, assay = "RNA", HTO = FALSE)

seurat_cfko <- umap_data$object

cM <- confusionMatrix(seurat_cfko$RNA_cluster,
                      seurat_cfko$final_celltype)

pheatmap::pheatmap(cM)
pheatmap::pheatmap(cM / rowSums(cM))


new_clust <- colnames(cM)[apply(cM, 1 , which.max)]

clust_names <- new_clust

names(clust_names) <- rownames(cM)

seurat_cfko$new_celltype <- clust_names[as.character(seurat_cfko$RNA_cluster)]


cm <- confusionMatrix(seurat_cfko$RNA_cluster, 
                      seurat_cfko$new_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

cm <- confusionMatrix(seurat_cfko$RNA_cluster, 
                      seurat_cfko$orig_combined_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

all_var_features <- lapply(unique(seurat_cfko$sample), function(sample_name){
  subset_seurat <- subset(seurat_cfko, subset = sample == sample_name)
  subset_seurat <- FindVariableFeatures(subset_seurat, nfeatures = 1500)
  return(VariableFeatures(subset_seurat))
})

all_var_features <- unique(unlist(all_var_features))

all_var_features <- VariableFeatures(seurat_cfko)

seurat_cfko <- run_all_clustifyr(seurat_data = seurat_cfko, 
                               save_dir = save_dir,
                               all_var_features = all_var_features,
                               mapping_file = mapping_file,
                               cor_cutoff = cor_cutoff)


cm <- confusionMatrix(seurat_cfko$orig_combined_celltype, 
                      seurat_cfko$RNA_combined_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)
### Rename cell types -------------------------------------------------------
seurat_name_mapping <- unique(seurat_cfko$RNA_combined_celltype)
names(seurat_name_mapping) <- seurat_name_mapping

seurat_name_mapping[seurat_name_mapping == "progenitor_like_cells"] = 
  "centroacinar_progenitor"

seurat_name_mapping[seurat_name_mapping == "transitional_to_acinar1"] = 
  "centroacinar_progenitor"

seurat_name_mapping[seurat_name_mapping == "transitional_to_acinar2"] = 
  "centroacinar_progenitor"

seurat_name_mapping[seurat_name_mapping == "activated_stellate"] = 
  "centroacinar_progenitor"

seurat_cfko$final_cfko_celltype <- seurat_name_mapping[as.character(seurat_cfko$RNA_combined_celltype)]

seurat_cfko$final_cfko_celltype <- factor(seurat_cfko$final_cfko_celltype,
                                      levels = names(all_colors2))


### Check -------------------------------------------------------------------

cm <- confusionMatrix(seurat_cfko$final_cfko_celltype, 
                      seurat_cfko$final_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)


cm <- confusionMatrix(seurat_cfko$final_celltype, 
                      seurat_cfko$final_cfko_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)


cm <- confusionMatrix(seurat_cfko$new_celltype, 
                      seurat_cfko$final_cfko_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

all_colors3 <- c("acinar" = "#D4405B",
                 "ductal" = "#A5903E",
                 "Prolif_acinar" = "#55A470",
                 "Prolif_ductal" = "#767FC9",
                 "centroacinar1" = "#297878",
                 "centroacinar2" = "#78295D")
plotDimRed(seurat_cfko, "final_celltype", plot_type = "rna.umap",
           color = all_colors3)

plotDimRed(seurat_cfko, "final_cfko_celltype", plot_type = "rna.umap",
           color = all_colors2)

cfko_meta <- seurat_cfko[[]] %>%
  dplyr::select(new_celltype, final_cfko_celltype, RNA_combined_celltype)

# Combine -----------------------------------------------------------------
colnames(wt_meta) <- c("new_celltype", "final_ind_celltype", 
                       "ind_full_celltype")
colnames(cfko_meta) <- c("new_celltype", "final_ind_celltype", 
                         "ind_full_celltype")
full_meta <- rbind(wt_meta, cfko_meta)

seurat_data <- AddMetaData(seurat_data, metadata = full_meta)

## UMAP ------------------------------------------------------------------

plot1 <- plotDimRed(seurat_data, "final_celltype", plot_type = "rna.umap",
           color = all_colors3, ggrastr = TRUE)[[1]] +
  ggplot2::ggtitle("Previous labels")

plot2 <- plotDimRed(seurat_data, "new_celltype", plot_type = "rna.umap",
           color = all_colors3, ggrastr = TRUE)[[1]]

plot3 <- plotDimRed(seurat_data, "final_ind_celltype", plot_type = "rna.umap",
           color = all_colors2, ggrastr = TRUE)[[1]] +
  ggplot2::ggtitle("New labels")

cowplot::plot_grid(plot1, plot2, plot3)

## Barplots --------------------------------------------------------------

all_barplots1 <- scAnalysisR::stacked_barplots(seurat_data,
                                              meta_col = "final_celltype",
                                              split_by = "sample",
                                              return_values = TRUE,
                                              color = all_colors3)$barplot +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5, hjust=1)) +
  ggplot2::ggtitle("Previous labels")

all_barplots2 <- scAnalysisR::stacked_barplots(seurat_data,
                                               meta_col = "new_celltype",
                                               split_by = "sample",
                                               return_values = TRUE,
                                               color = all_colors3)$barplot +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5, hjust=1))

all_barplots3 <- scAnalysisR::stacked_barplots(seurat_data,
                                               meta_col = "final_ind_celltype",
                                               split_by = "sample",
                                               return_values = TRUE,
                                               color = all_colors2)$barplot +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5, hjust=1)) +
  ggplot2::ggtitle("New labels")


cowplot::plot_grid(all_barplots1, all_barplots2, all_barplots3)


saveRDS(seurat_data, file.path(save_dir, "rda_obj/seurat_processed.rds"))
saveRDS(seurat_wt, file.path(save_dir, "rda_obj/seurat_wt.rds"))
saveRDS(seurat_cfko, file.path(save_dir, "rda_obj/seurat_cfko.rds"))

## Plots -----------------------------------------------------------------------

krt7_before <- featDistPlot(seurat_data, "KRT7", combine = FALSE, 
                     sep_by = "final_celltype",
                     color = all_colors3)[[1]] +
  ggplot2::ggtitle("Previous labels")

krt7_after <- featDistPlot(seurat_data, "KRT7", combine = FALSE, 
                            sep_by = "final_ind_celltype",
                            color = all_colors2)[[1]] +
  ggplot2::ggtitle("New labels")

seurat_sub <- subset(seurat_data, 
                     subset = final_ind_celltype %in%
                       c("centroacinar_progenitor",
                         "centroacinar"))

seurat_sub2 <- subset(seurat_data, 
                      subset = final_celltype %in%
                        c("centroacinar1",
                          "centroacinar2"))


aldh1a1_before <- featDistPlot(seurat_sub2, "ALDH1A1", col_by = "sample",
                        sep_by = "final_celltype", combine = FALSE,
                        color = sample_colors)[[1]] +
  ggplot2::ggtitle("Previous labels")

aldh1a1_after <- featDistPlot(seurat_sub, "ALDH1A1", col_by = "sample",
                               sep_by = "final_ind_celltype", combine = FALSE,
                               color = sample_colors)[[1]] +
  ggplot2::ggtitle("New labels")



cm <- confusionMatrix(seurat_data$final_celltype, 
                      seurat_data$final_ind_celltype)

cm <- cm / rowSums(cm)

save_heatmap <- pheatmap::pheatmap(cm, main = "original vs new celltypes")


cm <- confusionMatrix(seurat_data$new_celltype, 
                      seurat_data$final_ind_celltype)

cm <- cm / rowSums(cm)

#pheatmap::pheatmap(cm)


cm <- confusionMatrix(seurat_data$orig_combined_celltype, 
                      seurat_data$ind_full_celltype)

cm <- cm / rowSums(cm)

save_heatmap2 <- pheatmap::pheatmap(cm, main = "original vs new celltypes")


# Save plots
save_dir_images <- file.path(save_dir, "images",
                             "combined_celltypes")

ifelse(!dir.exists(save_dir_images), dir.create(save_dir_images), FALSE)

pdf(file.path(save_dir_images, "umap_comparison.pdf"),
    width = 15, height = 8)
print(cowplot::plot_grid(plot1, plot3))

dev.off()


pdf(file.path(save_dir_images, "barplot_comparison.pdf"),
    width = 12, height = 6)
print(cowplot::plot_grid(all_barplots1, all_barplots3)
)

dev.off()


pdf(file.path(save_dir_images, "cell_type_mapping.pdf"))

print(save_heatmap)
grid::grid.newpage()
print(save_heatmap2)

dev.off()


pdf(file.path(save_dir_images, "aldh1a1_expression.pdf"))

print(cowplot::plot_grid(aldh1a1_before, aldh1a1_after,
                         nrow = 2, ncol = 1))

dev.off()

pdf(file.path(save_dir_images, "krt7_expression.pdf"))
print(cowplot::plot_grid(krt7_before, krt7_after,
                         nrow = 2, ncol = 1))
dev.off()

seurat_data$celltype_sample <- paste(seurat_data$ind_full_celltype,
                                     seurat_data$sample, sep = "_")
pseudobulk <- AverageExpression(seurat_data,
                                features = VariableFeatures(seurat_data),
                                group.by = "celltype_sample")

pseudobulk <- pseudobulk$RNA


all_corr <- cor(pseudobulk)
all_corr <- round(all_corr, digits = 2)

pdf(file.path(save_dir_images, "celltype_correlations_updated.pdf"),
    width = 20, height = 20)
print(pheatmap::pheatmap(all_corr, display_numbers = all_corr, cluster_rows = FALSE,
                   cluster_cols = FALSE))

dev.off()

