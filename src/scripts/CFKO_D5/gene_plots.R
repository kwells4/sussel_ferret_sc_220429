library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)

source("src/scripts/functions.R")

# Setup ------------------------------------------------------------------------

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "CFKO_D5"

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

gene_path <- here("files/GSEA_signaling_pathways_with_orthologs.xlsx")

all_sheets <- openxlsx::getSheetNames(gene_path)

gene_lists <- lapply(all_sheets, function(x){
  gene_df <- openxlsx::readWorkbook(gene_path, sheet = x)
  return(unique(gene_df$gene_id))
})

names(gene_lists) <- all_sheets

gene_dfs <- lapply(all_sheets, function(x){
  gene_df <- openxlsx::readWorkbook(gene_path, sheet = x) %>%
    group_by(ortholog_id) %>%
    mutate(id = row_number()) %>%
    mutate(save_name = paste(ortholog_id, id, sep = "_"))
  return(gene_df)
})

names(gene_dfs) <- all_sheets


# Make directories
ifelse(!dir.exists(file.path(save_dir, "images", "gene_plots")),
       dir.create(file.path(save_dir, "images", "gene_plots")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images", "gene_plots", "modules")),
       dir.create(file.path(save_dir, "images", "gene_plots", "modules")),
       FALSE)

# Colors -----------------------------------------------------------------------
nColors <- length(unique(seurat_data$cluster_celltype))

cluster_colors <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(9, "Set1"))(nColors)
names(cluster_colors) <- unique(seurat_data$cluster_celltype)

byrnes_colors <- c("Acinar" = "#D4405B",
                   "Fev_hi" = "#75B140",
                   "Alpha" = "#885FCF",
                   "Ductal" = "#A5903E",
                   "Ngn3_pos" = "#C65CAC",
                   "Prolif_acinar" = "#55A470",
                   "Prolif_ductal" = "#767FC9",
                   "undetermined" = "#D3D3D3")

baron_colors <- c("ductal" = "#A5903E",
                  "acinar" = "#D4405B")

tabula_muris_colors <- c("pancreatic_acinar_cell" = "#D4405B",
                         "pancreatic_stellate_cell" = "#CC6D3A",
                         "pancreatic_ductal_cell" = "#A5903E",
                         "endothelial_cell" = "#45B0CF",
                         "undetermined" = "#D3D3D3")

qadir_colors <- c("immune_cells" = "#D3D3D3",
                  "progenitor_like_cells" = "#297878",
                  "transitional_to_acinar1" = "#874652",
                  "transitional_to_acinar2" = "#CC1B3B",
                  "centroacinar" = "#78295D",
                  "Small_ducts" = "#B37422")

muraro_colors <- c("ductal" = "#A5903E",
                   "acinar" = "#D4405B",
                   "undetermined" = "#D3D3D3")

krentz_colors <- c("AFP" = "#5C6AE0",
                   "Endo1" = "#45B0CF",
                   "undetermined" = "#D3D3D3")

all_colors <- c("Acinar" = "#D4405B",
                "Ductal" = "#A5903E",
                "ductal" = "#A5903E",
                "Prolif_acinar" = "#55A470",
                "Prolif_ductal" = "#767FC9",
                "progenitor_like_cells" = "#297878",
                "transitional_to_acinar1" = "#874652",
                "transitional_to_acinar2" = "#CC1B3B",
                "centroacinar" = "#78295D")

# Sample colors
sample_colors <- as.character(LaCroixColoR::lacroix_palette("Coconut", 10))
sample_colors[5] <- "#F4E3C7"
names(sample_colors) <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                          "CFKO_D14", "CFKO_D9", "CFKO_D7", "CFKO_D5", "CFKO_D2")
sample_colors <- sample_colors[c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                                 "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")]


# Add module score -------------------------------------------------------------

# Module score of all
for(i in names(gene_lists)){
  # print(i)
  # print(length(intersect(gene_lists[[i]], rownames(seurat_data))))
  # print(length(gene_lists[[i]]))
  seurat_data <- AddModuleScore(seurat_data, features = gene_lists[i],
                                name = paste0("module_", i))
}

all_plots <- lapply(names(gene_lists), function(x){
  module_name <- paste0("module_", make.names(x), 1)
  save_path <- file.path(save_dir, "images", "gene_plots", "modules",
                         paste0(module_name, ".pdf"))
  return_plot <- make_plots(seurat_object = seurat_data,
                            cluster_name = "cluster_celltype",
                            celltype_name = "RNA_combined_celltype",
                            cluster_colors = cluster_colors,
                            celltype_colors = all_colors,
                            gene = module_name,
                            assay = NULL, plot_type = "rna.umap",
                            save_name = save_path)
})

pdf(file.path(save_dir, "images", "gene_plots", "modules",
              "module_heatmaps.pdf"),
    height = 8, width = 7)

heatmaps <- lapply(names(gene_lists), function(x){
  print(plot_heatmap(seurat_data, gene_lists[[x]], "cluster_celltype",
                     colors = cluster_colors, main = x,
                     cluster_rows = TRUE))
  grid::grid.newpage()
})

dev.off()

# Gene umaps -------------------------------------------------------------------

invisible(lapply(names(gene_lists), function(x){
  print(x)
  ifelse(!dir.exists(file.path(save_dir, "images", "gene_plots", x)),
         dir.create(file.path(save_dir, "images", "gene_plots", x)),
         FALSE)
  
  invisible(lapply(gene_lists[[x]], function(y){
    if(y %in% rownames(seurat_data)){
      mapping_df <- gene_dfs[[x]]
      save_name <- mapping_df %>%
        filter(gene_id == y)
      save_name <- save_name$save_name
      
      save_path <- file.path(save_dir, "images", "gene_plots", x,
                             paste0(save_name, ".pdf"))
      
      return_plot <- make_plots(seurat_object = seurat_data,
                                cluster_name = "cluster_celltype",
                                celltype_name = "RNA_combined_celltype",
                                cluster_colors = cluster_colors,
                                celltype_colors = all_colors,
                                gene = y,
                                assay = NULL, plot_type = "rna.umap",
                                save_name = save_path)
      
      return("done")
      }
    }))
}))

