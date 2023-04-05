library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
library(ggridges)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "All_combined"

wt_samples <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14")

cfko_samples <- c("CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")

sample_colors <- as.character(LaCroixColoR::lacroix_palette("Coconut", 10))
sample_colors[5] <- "#F4E3C7"


names(sample_colors) <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                          "CFKO_D14", "CFKO_D9", "CFKO_D7", "CFKO_D5", "CFKO_D2")
sample_colors <- sample_colors[c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                                 "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")]



all_colors <- c("acinar" = "#D4405B",
                "ductal" = "#A5903E",
                "Prolif_acinar" = "#55A470",
                "Prolif_ductal" = "#767FC9",
                "progenitor_like_cells" = "#297878",
                "transitional_to_acinar" = "#874652",
                "centroacinar" = "#78295D")

phase_colors <- LaCroixColoR::lacroix_palette("CranRaspberry", n = 3)
names(phase_colors) <- c("S", "G1", "G2M")

sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))


merged_seurat$sample <- merged_seurat$orig.ident

# Gene sets
gene_path <- here("files/GSEA_signaling_pathways_with_orthologs.xlsx")

all_sheets <- openxlsx::getSheetNames(gene_path)

gene_lists <- lapply(all_sheets, function(x){
  gene_df <- openxlsx::readWorkbook(gene_path, sheet = x)
  #all_genes <- unique(gene_df$gene_id)
  return(unique(gene_df$gene_id))
})

names(gene_lists) <- all_sheets
all_samples <- c(wt_samples, cfko_samples)

# Umap of cell types -----------------------------------------------------------
all_plots <- lapply(all_samples, function(x){
  plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap",
             highlight_group = TRUE, group = x, meta_data_col = "sample",
             color = all_colors)[[1]]
})

names(all_plots) <- all_samples


all_cell_types <- cowplot::plot_grid(all_plots$WT_D2, all_plots$WT_D5,
                                     all_plots$WT_D7, all_plots$WT_D9,
                                     all_plots$WT_D14, NULL, all_plots$CFKO_D2,
                                     all_plots$CFKO_D5, all_plots$CFKO_D7,
                                     all_plots$CFKO_D9, all_plots$CFKO_D14,
                                     NULL, nrow = 4, ncol = 3)


pdf(file.path(sample_dir, "images", "combined_cell_type_umap.pdf"),
    width = 15, height = 10)
print(all_cell_types)

dev.off()

combined_plot <- plotDimRed(merged_seurat, "RNA_combined_celltype", 
                            plot_type = "rna.umap",
                            color = all_colors)[[1]]

pdf(file.path(sample_dir, "images", "combined_cell_type_umap_all.pdf"),
    width = 15, height = 10)
print(combined_plot)

dev.off()


wt_plot <- plotDimRed(merged_seurat, "RNA_combined_celltype", 
                      plot_type = "rna.umap",
                      color = all_colors,
                      highlight_group = TRUE, group = "WT",
                      meta_data_col = "genotype")[[1]]

pdf(file.path(sample_dir, "images", "combined_cell_type_umap_wt.pdf"),
    width = 15, height = 10)
print(wt_plot)

dev.off()

cfko_plot <- plotDimRed(merged_seurat, "RNA_combined_celltype", 
                        plot_type = "rna.umap",
                        color = all_colors,
                        highlight_group = TRUE, group = "CFKO",
                        meta_data_col = "genotype")[[1]]

pdf(file.path(sample_dir, "images", "combined_cell_type_umap_cfko.pdf"),
    width = 15, height = 10)
print(cfko_plot)

dev.off()
