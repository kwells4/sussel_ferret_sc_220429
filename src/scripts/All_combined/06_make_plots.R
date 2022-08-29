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

sample <- "All_combined"

wt_samples <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14")

cfko_samples <- c("CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")

sample_colors <- as.character(LaCroixColoR::lacroix_palette("Coconut", 10))
sample_colors[5] <- "#F4E3C7"
names(sample_colors) <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                          "CFKO_D14", "CFKO_D9", "CFKO_D7", "CFKO_D5", "CFKO_D2")
sample_colors <- sample_colors[c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                                 "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")]



all_colors <- c("Acinar" = "#D4405B",
                "acinar" = "#D4405B",
                "Ductal" = "#A5903E",
                "ductal" = "#A5903E",
                "Prolif_acinar" = "#55A470",
                "Prolif_ductal" = "#767FC9",
                "progenitor_like_cells" = "#297878",
                "transitional_to_acinar1" = "#874652",
                "transitional_to_acinar2" = "#CC1B3B",
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


# Add psuedotime to object
psuper_obj_wt <- readRDS(file.path(sample_dir, "rda_obj", "wt_psupertime.rds"))

psuper_obj_cfko <- readRDS(file.path(sample_dir, "rda_obj",
                                   "cfko_psupertime.rds"))

project_info_wt <- psuper_obj_wt$proj_dt %>%
  tibble::column_to_rownames("cell_id")

project_info_cfko <- psuper_obj_cfko$proj_dt %>%
  tibble::column_to_rownames("cell_id")

project_info <- rbind(project_info_cfko, project_info_wt)

merged_seurat <- AddMetaData(merged_seurat, metadata = project_info)

all_samples <- c(wt_samples, cfko_samples)

cluster_plots <- lapply(all_samples, function(x){
  plotDimRed(merged_seurat, col_by = "RNA_cluster", plot_type = "rna.umap",
             highlight_group = TRUE, group = x, meta_data_col = "sample")[[1]]
})


names(cluster_plots) <- all_samples

pdf(file.path(sample_dir, "images", "cluster_plots_separate.pdf"),
    width = 8, height = 8)
print(cluster_plots)
dev.off()

pdf(file.path(sample_dir, "images", "cluster_plots_combined.pdf"),
    width = 15, height = 12)

print(cowplot::plot_grid(cluster_plots$WT_D2, cluster_plots$WT_D5,
                   cluster_plots$WT_D7, cluster_plots$CFKO_D2,
                   cluster_plots$CFKO_D5, cluster_plots$CFKO_D7,
                   cluster_plots$WT_D9, cluster_plots$WT_D14,
                   NULL, cluster_plots$CFKO_D9, cluster_plots$CFKO_D14,
                   NULL, nrow = 4, ncol = 3))

dev.off()


pdf(file.path(sample_dir, "images", "pseudotime_plots.pdf"),
    width = 8, height = 8)

psupertime_plot_wt <- plotDimRed(merged_seurat, col_by = "psuper",
                                 plot_type = "rna.umap",
                                 highlight_group = TRUE, group = "WT",
                                 meta_data_col = "genotype")


psupertime_plot_cfko <- plotDimRed(merged_seurat, col_by = "psuper",
                                 plot_type = "rna.umap",
                                 highlight_group = TRUE, group = "CFKO",
                                 meta_data_col = "genotype")

print(psupertime_plot_wt)
print(psupertime_plot_cfko)

dev.off()

# DE heatmaps ------------------------------------------------------------------

# First some building

merged_seurat$sample <- factor(merged_seurat$sample,
                               levels = names(sample_colors))

# Read in DE genes
de_directory <- file.path(sample_dir, "files", "DE")

## Sample DE -------------------------------------------------------------------

de_tests <- c("D2_all", "D5_all", "D7_all", "D9_all", "D14_all")

# Make plots for all tests 
invisible(lapply(de_tests, function(de_test){
  print(de_test)
  
  excel_file <- file.path(de_directory, paste0(de_test, ".xlsx"))
  
  excel_sheets <- openxlsx::getSheetNames(excel_file)
  
  de_genes <- lapply(excel_sheets, function(x){
    de_df <- openxlsx::readWorkbook(excel_file, sheet = x)
    de_df$up_in <- x
    return(de_df)
  })
  
  names(de_genes) <- excel_sheets
  
  de_genes <- do.call(rbind, de_genes)
  
  
  fig_height <- round(nrow(de_genes) / 10)
  
  # Plot heatmap across all samples
  pdf(file.path(sample_dir, "images", "sample_heatmaps",
                paste0(de_test, "_heatmap_average.pdf")))
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample", average_expression = TRUE,
                     colors = sample_colors, plot_rownames = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "sample_heatmaps",
                paste0(de_test, "_heatmap_average_labeled.pdf")),
      width = 8, height = fig_height)
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample", average_expression = TRUE,
                     colors = sample_colors, plot_rownames = TRUE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "sample_heatmaps",
                paste0(de_test, "_heatmap_cells.pdf")))
  
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample", average_expression = FALSE,
                     colors = sample_colors, plot_rownames = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "sample_heatmaps",
                paste0(de_test, "_heatmap_cells_labeled.pdf")),
      width = 8, height = fig_height)
  
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample", average_expression = FALSE,
                     colors = sample_colors, plot_rownames = TRUE))
  
  dev.off()
  
}))

## Celltype DE -----------------------------------------------------------------
de_tests <- c("D2_combined_celltype", "D5_combined_celltype",
              "D7_combined_celltype", "D9_combined_celltype",
              "D14_combined_celltype")


merged_seurat$sample_celltype <- paste(merged_seurat$sample,
                                       merged_seurat$RNA_combined_celltype,
                                       sep = "_")

# Make plots for all tests 
invisible(lapply(de_tests, function(de_test){
  print(de_test)
  
  excel_file <- file.path(de_directory, paste0(de_test, ".xlsx"))
  
  excel_sheets <- openxlsx::getSheetNames(excel_file)
  
  de_genes <- lapply(excel_sheets, function(x){
    de_df <- openxlsx::readWorkbook(excel_file, sheet = x)
    de_df$up_in <- x
    return(de_df)
  })
  
  names(de_genes) <- excel_sheets
  
  de_genes <- do.call(rbind, de_genes)
  
  
  fig_height <- round(nrow(de_genes) / 10)
  
  # Keep only the one time point
  subset_seurat <- subset(merged_seurat,
                          subset = day == gsub("_.*", "", de_test))
  
  # Make meta data and colors
  meta_df <- subset_seurat[[c("sample", "RNA_combined_celltype",
                              "sample_celltype")]]
  meta_ave <- data.frame(table(paste0(subset_seurat$sample, ";",
                                      subset_seurat$RNA_combined_celltype)))
  
  meta_ave$RNA_combined_celltype <- gsub(".*_D[0-9]+;", "", meta_ave$Var1)
  meta_ave$sample <- gsub(";.*$", "", meta_ave$Var1)
  meta_ave$Var1 <- NULL
  meta_ave$Freq <- NULL
  meta_ave$sample_celltype <- paste(meta_ave$sample,
                                    meta_ave$RNA_combined_celltype,
                                    sep = "_")
  
  meta_ave$sample_celltype <- factor(meta_ave$sample_celltype)
  
  rownames(meta_ave) <- meta_ave$sample_celltype
  color_list <- list("sample" = sample_colors,
                     "RNA_combined_celltype" = all_colors)

  
  # Plot heatmap across all cell types
  pdf(file.path(sample_dir, "images", "cell_type_heatmaps",
                paste0(de_test, "_heatmap_average.pdf")))
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample_celltype", average_expression = TRUE,
                     colors = sample_colors, plot_rownames = FALSE,
                     meta_df = meta_ave, color_list = color_list,
                     plot_meta_col = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "cell_type_heatmaps",
                paste0(de_test, "_heatmap_average_labeled.pdf")),
      width = 8, height = fig_height)
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample_celltype", average_expression = TRUE,
                     colors = sample_colors, plot_rownames = TRUE,
                     meta_df = meta_ave, color_list = color_list,
                     plot_meta_col = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "cell_type_heatmaps",
                paste0(de_test, "_heatmap_cells.pdf")))
  
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample_celltype", average_expression = FALSE,
                     colors = sample_colors, plot_rownames = FALSE,
                     meta_df = meta_df, color_list = color_list,
                     plot_meta_col = FALSE))
  
  dev.off()
  
  pdf(file.path(sample_dir, "images", "cell_type_heatmaps",
                paste0(de_test, "_heatmap_cells_labeled.pdf")),
      width = 8, height = fig_height)
  
  print(plot_heatmap(merged_seurat, gene_list = unique(de_genes$gene),
                     meta_col = "sample_celltype", average_expression = FALSE,
                     colors = sample_colors, plot_rownames = TRUE,
                     meta_df = meta_df, color_list = color_list,
                     plot_meta_col = FALSE))
  
  dev.off()
  
  
  # Plot heatmaps across one cell type all samples
  de_genes$celltype <- gsub(".*_D[0-9]+_", "", de_genes$up_in)
  
  invisible(lapply(unique(de_genes$celltype), function(y){
    celltype_seurat <- subset(merged_seurat,
                              subset = RNA_combined_celltype == y)
    
    celltype_seurat$sample <- droplevels(celltype_seurat$sample)
    
    new_de_genes <- de_genes %>%
      dplyr::filter(celltype == y)
    
    fig_height <- round(nrow(new_de_genes) / 10)
    
    # Plot heatmap across all samples
    pdf(file.path(sample_dir, "images", "cell_type_heatmaps",
                  paste0(de_test, "_", y,
                                               "_heatmap_average.pdf")))
    print(plot_heatmap(celltype_seurat, gene_list = unique(new_de_genes$gene),
                       meta_col = "sample", average_expression = TRUE,
                       colors = sample_colors, plot_rownames = FALSE))
    
    dev.off()
    
    pdf(file.path(sample_dir, "images", "cell_type_heatmaps",
                  paste0(de_test, "_", y, "_heatmap_average_labeled.pdf")),
        width = 8, height = fig_height)
    print(plot_heatmap(celltype_seurat, gene_list = unique(new_de_genes$gene),
                       meta_col = "sample", average_expression = TRUE,
                       colors = sample_colors, plot_rownames = TRUE))
    
    dev.off()
    
    pdf(file.path(sample_dir, "images", "cell_type_heatmaps",
                  paste0(de_test, "_", y,
                                               "_heatmap_cells.pdf")))
    
    print(plot_heatmap(celltype_seurat, gene_list = unique(new_de_genes$gene),
                       meta_col = "sample", average_expression = FALSE,
                       colors = sample_colors, plot_rownames = FALSE))
    
    dev.off()
    
    pdf(file.path(sample_dir, "images", "cell_type_heatmaps",
                  paste0(de_test, "_", y, "_heatmap_cells_labeled.pdf")),
        width = 8, height = fig_height)
    
    print(plot_heatmap(celltype_seurat, gene_list = unique(new_de_genes$gene),
                       meta_col = "sample", average_expression = FALSE,
                       colors = sample_colors, plot_rownames = TRUE))
    
    dev.off()
    
  }))
  
  
}))



# GOI genesets -----------------------------------------------------------------

## UMAP ------------------------------------------------------------------------
# Color palette - 2 versions. One my default. One that is grey/red

names(gene_lists) <- make.names(names(gene_lists))

# add module scores
gene_scores <- lapply(names(gene_lists), function(x){
  print(x)
  genes <- gene_lists[[x]]
  genes <- list(genes)
  new_seurat <- AddModuleScore(merged_seurat,
                               features = genes,
                               name = paste0(x, "_"))
  return_df <- new_seurat[[paste0(x, "_1")]]
  return(return_df)  
})

gene_scores <- do.call(cbind, gene_scores)

merged_seurat <- AddMetaData(merged_seurat, gene_scores)

all_umaps <- lapply(names(gene_lists), function(x){
  plot_one <- plotDimRed(merged_seurat, col_by = paste0(x, "_1"),
                         plot_type = "rna.umap")[[1]]
  plot_two <- plotDimRed(merged_seurat, col_by = paste0(x, "_1"),
                         plot_type = "rna.umap", color = c("grey", "red"))[[1]]
  
  return(list("default_color" = plot_one, "grey_red" = plot_two))
})

names(all_umaps) <- names(gene_lists)

combined_umaps <- cowplot::plot_grid(all_umaps$EP.genes$default_color,
                                     all_umaps$TGFB$default_color,
                                     all_umaps$WNT$default_color,
                                     all_umaps$PI3K.AKT.mTOR$default_color,
                                     all_umaps$PTEN$default_color,
                                     all_umaps$NOTCH$default_color,
                                     nrow = 3, ncol = 2)

combined_umaps_two <- cowplot::plot_grid(all_umaps$EP.genes$grey_red,
                                         all_umaps$TGFB$grey_red,
                                         all_umaps$WNT$grey_red,
                                         all_umaps$PI3K.AKT.mTOR$grey_red,
                                         all_umaps$PTEN$grey_red,
                                         all_umaps$NOTCH$grey_red,
                                         nrow = 3, ncol = 2)


pdf(file.path(sample_dir, "images", "gene_set_umap.pdf"),
    width = 8, height = 8)
print(combined_umaps)
print(combined_umaps_two)

dev.off()



## Pseudotime ------------------------------------------------------------------
make_pseudotime_plots <- function(seurat_object, genotype_select, geneset,
                                  pseudotime_column = "psuper"){
  meta_data <- seurat_object[[]] %>%
    dplyr::filter(genotype == genotype_select) %>%
    dplyr::select(dplyr::all_of(c(pseudotime_column, geneset)))
    
  plot <- ggplot2::ggplot(meta_data, ggplot2::aes_string(x = pseudotime_column,
                                                 y = geneset,
                                                 color = pseudotime_column)) +
    ggplot2::geom_point() +
    viridis::scale_color_viridis(option = "magma")
  
  return(plot)
}

all_plots <- lapply(names(gene_lists), function(x){
  geneset <- paste0(x, "_1")
  wt_plot <- make_pseudotime_plots(merged_seurat, genotype_select = "WT",
                                   geneset = geneset)
  
  cfko_plot <- make_pseudotime_plots(merged_seurat, genotype_select = "CFKO",
                                     geneset = geneset)
  
  return(list("wt_plot" = wt_plot, "cfko_plot" = cfko_plot))
  
})

names(all_plots) <- names(gene_lists)

wt_pseudotime <- cowplot::plot_grid(all_plots$EP.genes$wt_plot,
                                    all_plots$TGFB$wt_plot,
                                    all_plots$WNT$wt_plot,
                                    all_plots$PI3K.AKT.mTOR$wt_plot,
                                    all_plots$PTEN$wt_plot,
                                    all_plots$NOTCH$wt_plot,
                                    nrow = 3, ncol = 2)

cfko_pseudotime <- cowplot::plot_grid(all_plots$EP.genes$cfko_plot,
                                    all_plots$TGFB$cfko_plot,
                                    all_plots$WNT$cfko_plot,
                                    all_plots$PI3K.AKT.mTOR$cfko_plot,
                                    all_plots$PTEN$cfko_plot,
                                    all_plots$NOTCH$cfko_plot,
                                    nrow = 3, ncol = 2)

pdf(file.path(sample_dir, "images", "gene_set_pseudotime.pdf"),
    width = 8, height = 8)
print(wt_pseudotime)
print(cfko_pseudotime)

dev.off()

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


# Umap of cell types -----------------------------------------------------------
color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
                                                                      "Set1"))

all_celltypes <- unique(merged_seurat$RNA_combined_celltype_ref_name)

celltype_cols <- color_palette(length(all_celltypes))

names(celltype_cols) <- all_celltypes

all_plots <- lapply(all_samples, function(x){
  plotDimRed(merged_seurat, "RNA_combined_celltype_ref_name",
             plot_type = "rna.umap",
             highlight_group = TRUE, group = x, meta_data_col = "sample",
             color = celltype_cols)[[1]]
})

names(all_plots) <- all_samples


all_cell_types <- cowplot::plot_grid(all_plots$WT_D2, all_plots$WT_D5,
                                     all_plots$WT_D7, all_plots$WT_D9,
                                     all_plots$WT_D14, NULL, all_plots$CFKO_D2,
                                     all_plots$CFKO_D5, all_plots$CFKO_D7,
                                     all_plots$CFKO_D9, all_plots$CFKO_D14,
                                     NULL, nrow = 4, ncol = 3)


pdf(file.path(sample_dir, "images", "combined_cell_type_reference_umap.pdf"),
    width = 15, height = 15)
print(all_cell_types)

dev.off()

## Violin plots ----------------------------------------------------------------
pdf(file.path(sample_dir, "images", "DE_violin_plots.pdf"),
    height = 11, width = 8)

# D2
plot(featDistPlot(merged_seurat, geneset = c("ATOX1", "CRABP1",
                                             "GATA3", "PAX6"),
                  sep_by = "sample", color = sample_colors))

# D5
plot(featDistPlot(merged_seurat, geneset = c("FN1", "PALLD",
                                             "TMSB10", "HES1"),
                  sep_by = "sample", color = sample_colors))

# D7
plot(featDistPlot(merged_seurat, geneset = c("LAMB1", "PDIA3",
                                             "GAS6", "ITGB1"),
                  sep_by = "sample", color = sample_colors))

# D9
plot(featDistPlot(merged_seurat, geneset = c("DKK1", "IGFBP7",
                                             "TAGLN", "ANKRD1"),
                  sep_by = "sample", color = sample_colors))

# D14
plot(featDistPlot(merged_seurat, geneset = c("SLC25A13", "SPP1",
                                             "CTNNB1", "LASP1"),
                  sep_by = "sample", color = sample_colors))

dev.off()


# Umap of cell cycle phase -----------------------------------------------------

all_plots <- lapply(all_samples, function(x){
  plotDimRed(merged_seurat, "Phase", plot_type = "rna.umap",
             highlight_group = TRUE, group = x, meta_data_col = "sample",
             color = phase_colors)[[1]]
})

names(all_plots) <- all_samples


all_cell_types <- cowplot::plot_grid(all_plots$WT_D2, all_plots$WT_D5,
                                     all_plots$WT_D7, all_plots$WT_D9,
                                     all_plots$WT_D14, NULL, all_plots$CFKO_D2,
                                     all_plots$CFKO_D5, all_plots$CFKO_D7,
                                     all_plots$CFKO_D9, all_plots$CFKO_D14,
                                     NULL, nrow = 4, ncol = 3)


pdf(file.path(sample_dir, "images", "cell_cycle_umap.pdf"),
    width = 10, height = 10)
print(all_cell_types)

dev.off()

pdf(file.path(sample_dir, "images", "cell_cycle_umap_individual.pdf"),
    width = 6, height = 6)

print(all_plots)

dev.off()



saveRDS(merged_seurat, file.path(sample_dir, "rda_obj/seurat_processed.rds"))
