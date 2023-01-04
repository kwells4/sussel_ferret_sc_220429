library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
library(scran)

source("src/scripts/functions.R")

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "All_combined"

pval <- 0.05
logfc <- 0.5
cell_cutoff <- 20

all_colors <- c("acinar" = "#D4405B",
                "ductal" = "#A5903E",
                "Prolif_acinar" = "#55A470",
                "Prolif_ductal" = "#767FC9",
                "progenitor_like_cells" = "#297878",
                "transitional_to_acinar" = "#874652",
                "centroacinar" = "#78295D")


mapping_file <- read.csv(here("files/species_mapping_file.csv"))

sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

# Updated cluster calls in "best_celltype_markers.R" but get rid of activated
# ductal

merged_seurat$final_celltype <- merged_seurat$new_celltype

merged_seurat$final_celltype[merged_seurat$final_celltype == 
                               "activated_ductal"] <- "ductal"

# Scran find markers with both, make a heatmap
merged_sce <- as.SingleCellExperiment(merged_seurat)

scran_markers <- scran::scoreMarkers(merged_sce, merged_sce$RNA_combined_celltype,
                                   block = merged_sce$orig.ident)

scran_markers$acinar %>%
  data.frame() %>%
  top_n(n = -5, wt = rank.AUC)

# Get top vals 

top_scran_markers <- lapply(scran_markers, function(marker_list){
  return_markers <- marker_list %>%
    data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(mean.logFC.cohen > 0) %>%
    dplyr::top_n(n = -50, wt = rank.AUC)
})

top_scran_markers <- do.call(rbind, top_scran_markers)

scAnalysisR::plot_heatmap(merged_seurat, 
                          gene_list = unique(top_scran_markers$gene),
                          meta_col = "RNA_combined_celltype",
                          colors = all_colors,
                          average_expression = TRUE,
                          plot_rownames = FALSE)


Idents(merged_seurat) <- "RNA_combined_celltype"

# conserved_markers <- lapply(unique(Idents(merged_seurat)), function(cell_type){
#   all_markers <- FindConservedMarkers(merged_seurat, cell_type,
#                                       grouping.var = "orig.ident")
#   
#   all_markers <- all_markers %>%
# })
# 
# conserved_markers <- conserved_markers %>%
#   dplyr::select(contains("log2FC"), contains("p_val_adj"))
# 
# all_pvals <- conserved_markers %>%
#   dplyr::select(contains("p_val_adj")) %>%
#   dplyr::mutate(avg = )
# 
# conserved_markers <- FindConservedMarkers(merged_seurat, "acinar",
#                                     grouping.var = "orig.ident")


all_markers <- FindAllMarkers(merged_seurat)

all_markers_new <- all_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 50, wt = avg_log2FC)

scAnalysisR::plot_heatmap(merged_seurat, 
                          gene_list = unique(all_markers_new$gene),
                          meta_col = "RNA_combined_celltype",
                          colors = all_colors,
                          average_expression = TRUE,
                          plot_rownames = FALSE)


# Overlap with panglao db
ref_path <- file.path("/Users/wellskr/Documents/Analysis",
                      "references/single_cell_references/gene_lists",
                      "PanglaoDB_markers_27_Mar_2020.tsv.gz")

ref_data <- read.table(ref_path, sep = "\t", header = TRUE) %>%
  dplyr::select(species, official.gene.symbol, cell.type) %>%
  dplyr::rename(gene = official.gene.symbol)

ref_data <- merge(ref_data, all_markers_new, all.x = FALSE, all.y = TRUE)

# Scran find markers with both, make a heatmap
merged_sce <- as.SingleCellExperiment(merged_seurat)

all_markers <- scran::scoreMarkers(merged_sce, merged_sce$final_celltype,
                                   block = merged_sce$orig.ident)

all_markers$acinar %>%
  data.frame() %>%
  arrange(rank.AUC) %>%
  head()

# Get top vals 

featDistPlot(merged_seurat, geneset = "IGFBP7", combine = FALSE,
             sep_by = "final_celltype", col_by = "final_celltype")

featDistPlot(merged_seurat, geneset = "IGFBP7", combine = FALSE,
             sep_by = "orig.ident", col_by = "final_celltype")


Idents(merged_seurat) <- "final_celltype"
all_markers <- FindConservedMarkers(merged_seurat, "acinar",
                                    grouping.var = "orig.ident")


all_markers <- FindAllMarkers(merged_seurat)
