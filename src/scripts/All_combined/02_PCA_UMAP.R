library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(clustree)
library(harmony)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "All_combined"

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

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj/seurat_start.rds"))

# PCA --------------------------------------------------------------------------


# PCA of gene expression
seurat_data <- PCA_dimRed(seurat_data, assay = seurat_assay)

RNA_plots <- plot_PCA(HTO = HTO, assay = seurat_assay,
                      sample_object = seurat_data,
                      jackstraw = FALSE)

pdf(file.path(save_dir, "images/RNA_pca.pdf"))
RNA_plots
dev.off()

if(ADT){
  # PCA of surface protein
  seurat_data <- PCA_dimRed(seurat_data, assay = "ADT")
  
  ADT_plots <- plot_PCA(HTO = HTO, assay = "ADT", sample_object = seurat_data)
  
  pdf(paste0(save_dir, "images/RNA_pca.pdf"))
  plot(ADT_plots)
  dev.off()
}

# UMAP -------------------------------------------------------------------------

RNA_pcs <- 32
ADT_pcs <- 8

seurat_data$RNA_cluster <- NULL

# Remove previous clustering
remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                               colnames(seurat_data[[]]))]

for (i in remove_cols){
  seurat_data[[i]] <- NULL
}

set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = seurat_assay, HTO = HTO)

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

seurat_data <- FindClusters(seurat_data, resolution = c(0.2, 0.5, 0.8, 1,
                                                        1.2, 1.4))
clustree(seurat_data)

# UMAP of gene expression
set.seed(0)
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 1.2, assay = seurat_assay, HTO = HTO)

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

seurat_data$day <- gsub(".*_", "", seurat_data$sample)
seurat_data$genotype <- gsub("_D[0-9]+", "", seurat_data$sample)

plot_values <- c("sample", "day", "genotype", "RNA_celltype_byrnes",
                 "RNA_baron_celltype", "RNA_tabula_muris_celltype",
                 "RNA_combned_celltype", "RNA_baron_human_celltype",
                 "RNA_qadir_celltype", "RNA_muraro_celltype")

all_plots <- plotDimRed(seurat_data, col_by = plot_values,
                        plot_type = "rna.umap")


individual_plots <- lapply(unique(seurat_data$sample), function(x){
  plotDimRed(seurat_data, col_by = "sample", plot_type = "rna.umap",
             highlight_group = TRUE, group = x, meta_data_col = "sample")
})

seurat_data <- BuildClusterTree(seurat_data, dims = 1:RNA_pcs)
PlotClusterTree(seurat_data)

seurat_data$uncorrected_cluster <- seurat_data$RNA_cluster

seurat_data$RNA_cluster <- NULL

cluster_data <- data.frame(table(paste0(seurat_data$uncorrected_cluster,
                                    "_", seurat_data$RNA_combined_celltype))) %>%
  mutate(cluster = gsub("_.*", "", Var1)) %>%
  group_by(cluster) %>%
  add_count(wt = Freq) %>%
  mutate(percent = Freq/n) %>%
  slice_max(percent) %>%
  data.frame()

# Batch correction -------------------------------------------------------------

seurat_data <- RunHarmony(seurat_data, c("genotype"),
                          plot_convergence = TRUE,
                          theta = 3)

harmony_plots <- plotDimRed(seurat_data, col_by = "sample",
                           plot_type = "harmony")

# UMAP -------------------------------------------------------------------------

RNA_pcs <- 30
ADT_pcs <- 8

# Remove previous clustering
remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                               colnames(seurat_data[[]]))]

for (i in remove_cols){
  seurat_data[[i]] <- NULL
}


set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = seurat_assay, HTO = HTO,
                         reduction = "harmony")

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

seurat_data <- FindClusters(seurat_data, resolution = c(0.5, 0.8, 1, 1.2))
clustree(seurat_data)

# UMAP of gene expression
set.seed(0)
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 1.2, assay = seurat_assay, HTO = HTO,
                         reduction = "harmony")


seurat_data <- umap_data$object

gene_plots <- umap_data$plots

all_data <- data.frame(table(paste0(seurat_data$RNA_cluster,
                         "_", seurat_data$RNA_combined_celltype))) %>%
  mutate(cluster = gsub("_.*", "", Var1)) %>%
  group_by(cluster) %>%
  add_count(wt = Freq) %>%
  mutate(percent = Freq/n) %>%
  slice_max(percent) %>%
  data.frame()

all_plots_h <- plotDimRed(seurat_data, col_by = plot_values,
                        plot_type = "harmony.umap")

individual_plots_h <- lapply(unique(seurat_data$sample), function(x){
  plotDimRed(seurat_data, col_by = "sample", plot_type = "harmony.umap",
             highlight_group = TRUE, group = x, meta_data_col = "sample")
})


seurat_data <- BuildClusterTree(seurat_data, dims = 1:RNA_pcs)
PlotClusterTree(seurat_data)


if(ADT){
  # UMAP of surface protein
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = ADT_pcs,
                           resolution = 0.6, assay = "ADT", HTO = TRUE)
  
  seurat_data <- umap_data$object
  
  adt_plots <- umap_data$plots
  
  
  # UMAP of combined
  # Identify multimodal neighbors. These will be stored in the neighbors slot, 
  # and can be accessed using bm[['weighted.nn']]
  # The WNN graph can be accessed at bm[["wknn"]], 
  # and the SNN graph used for clustering at bm[["wsnn"]]
  # Cell-specific modality weights can be accessed at bm$RNA.weight
  if(SCT){
    pca_slot <- "sctpca"
    weight_name <- "SCT.weight"
  } else{
    pca_slot <- "pca"
    weight_name <- "RNA.weight"
  }
  seurat_data <- FindMultiModalNeighbors(
    seurat_data, reduction.list = list(pca_slot, "apca"), 
    dims.list = list(1:RNA_pcs, 1:ADT_pcs),
    modality.weight.name = c(weight_name, "ADT.weight")
  )
  
  
  seurat_data <- RunUMAP(seurat_data, nn.name = "weighted.nn",
                         reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  seurat_data <- FindClusters(seurat_data, graph.name = "wsnn",
                              algorithm = 3, resolution = 2, verbose = FALSE)
  
  seurat_data[["combined_cluster"]] <- Idents(seurat_data)
  col_by_list <- c("combined_cluster", "orig.ident")
  if(HTO){
    col_by_list <- c(col_by_list, "HTO_classification")
  }
  save_plot <- file.path(save_dir, "images/combined_umap.pdf")
  plot_list <- plotDimRed(sample_object = seurat_data,
                          save_plot = NULL,
                          col_by = col_by_list, return_plot = TRUE,
                          plot_type = "wnn.umap")
}

saveRDS(seurat_data, file.path(save_dir, "rda_obj/seurat_processed.rds"))
