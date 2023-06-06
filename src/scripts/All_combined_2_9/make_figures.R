library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(viridis)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "All_combined_2_9"

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

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

fig_dir <- file.path(save_dir, "images", "final_figures")

ifelse(!dir.exists(fig_dir), dir.create(fig_dir), FALSE)

seurat_data$celltype <- seurat_data$final_ind_celltype

seurat_data$sample <- factor(seurat_data$sample,
                             levels = names(sample_colors))

# Main figures -----------------------------------------------------------------

## UMAPs -----------------------------------------------------------------------
# Separated by genotype and differention day, colored by celltype
all_plots <- lapply(names(sample_colors), function(x){
  plotDimRed(seurat_data, col_by = "celltype", color = all_colors2,
             highlight_group = TRUE, group = x, meta_data_col = "orig.ident",
             plot_type = "rna.umap", ggrastr = TRUE)[[1]]
})

names(all_plots) <- names(sample_colors)

all_cell_types <- cowplot::plot_grid(all_plots$WT_D2, all_plots$WT_D5,
                                     all_plots$WT_D7, all_plots$WT_D9,
                                     all_plots$CFKO_D2,
                                     all_plots$CFKO_D5, all_plots$CFKO_D7,
                                     all_plots$CFKO_D9, nrow = 2, ncol = 4)

pdf(file.path(fig_dir, "combined_cell_type_umap.pdf"),
    width = 20, height = 10)
print(all_cell_types)

dev.off()

## Barplots --------------------------------------------------------------------
# Barplots of cell types by day

all_barplots <- scAnalysisR::stacked_barplots(seurat_data,
                                              meta_col = "celltype",
                                              split_by = "sample",
                                              return_values = TRUE,
                                              color = all_colors2)$barplot +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust=1))

pdf(file.path(fig_dir, "cell_type_barplot.pdf"),
    width = 4, height = 3)
print(all_barplots)

dev.off()

## DE genes --------------------------------------------------------------------
# DE genes for each day, plot all days, group by genotype
de_tests <- c("D2_all", "D5_all", "D7_all", "D9_all")

de_directory <- file.path(save_dir, "files", "DE")

# Make plots for all tests 
invisible(lapply(de_tests, function(de_test){
  print(de_test)
  
  excel_file <- file.path(de_directory, paste0(de_test, ".xlsx"))
  
  excel_sheets <- openxlsx::getSheetNames(excel_file)
  
  excel_sheets <- excel_sheets[!grepl("gse", excel_sheets)]
  
  de_genes <- lapply(excel_sheets, function(x){
    de_df <- openxlsx::readWorkbook(excel_file, sheet = x)
    de_df$up_in <- x
    return(de_df)
  })
  
  names(de_genes) <- excel_sheets
  
  de_genes <- do.call(rbind, de_genes)
  
  
  fig_height <- round(nrow(de_genes) / 10)
  
  # Plot heatmap across all samples
  pdf(file.path(fig_dir,
                paste0(de_test, "_heatmap_average.pdf")))
  
  print(plot_heatmap(seurat_data, gene_list = unique(de_genes$gene),
                     meta_col = "sample", average_expression = TRUE,
                     colors = sample_colors, plot_rownames = FALSE))
  
  dev.off()
  
}))

### Cell type DE ---------------------------------------------------------------
seurat_data$sample_celltype <- paste(seurat_data$sample,
                                     seurat_data$celltype,
                                     sep = "_")

# Grab out results for that test
de_test <- "D9_combined_celltype"

excel_file <- file.path(de_directory, paste0(de_test, ".xlsx"))

excel_sheets <- openxlsx::getSheetNames(excel_file)

excel_sheets <- excel_sheets[!grepl("gse", excel_sheets)]


de_genes <- lapply(excel_sheets, function(x){
  de_df <- openxlsx::readWorkbook(excel_file, sheet = x)
  de_df$up_in <- x
  return(de_df)
})

names(de_genes) <- excel_sheets

de_genes <- do.call(rbind, de_genes)


#### Progenitor centroacinar DE at day 9 ---------------------------------------
keep_celltype <- "centroacinar_progenitor"

test_de <- de_genes %>%
  dplyr::filter(grepl(paste0(keep_celltype, "$"), up_in))

fig_height <- round(nrow(test_de) / 10)

# Plot heatmaps across one cell type all samples
test_de$celltype <- gsub(".*_D[0-9]+_", "", 
                         test_de$up_in)

celltype_seurat <- subset(seurat_data,
                          subset = celltype == keep_celltype)
  
celltype_seurat$sample <- droplevels(celltype_seurat$sample)

# Plot heatmap across all samples
pdf(file.path(fig_dir,
              paste0(de_test, "_", keep_celltype,
                     "_heatmap_average.pdf")))
print(plot_heatmap(celltype_seurat, gene_list = unique(test_de$gene),
                   meta_col = "sample", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = FALSE))

dev.off()


# Centroacinar DE at day 9
keep_celltype <- "centroacinar"

test_de <- de_genes %>%
  dplyr::filter(grepl(paste0(keep_celltype, "$"), up_in))

fig_height <- round(nrow(test_de) / 10)

# Plot heatmaps across one cell type all samples
test_de$celltype <- gsub(".*_D[0-9]+_", "", 
                         test_de$up_in)

celltype_seurat <- subset(seurat_data,
                          subset = celltype == keep_celltype)

celltype_seurat$sample <- droplevels(celltype_seurat$sample)

# Plot heatmap across all samples
pdf(file.path(fig_dir,
              paste0(de_test, "_", keep_celltype,
                     "_heatmap_average.pdf")))
print(plot_heatmap(celltype_seurat, gene_list = unique(test_de$gene),
                   meta_col = "sample", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = FALSE))

dev.off()

## Violin plots ----------------------------------------------------------------
violin_genes <- c("ATOX1", "GATA3", "LAMB1", "PAX6",
                  "PDIA3", "TMSB10", "DKK1", "HES1",
                  "IGFBP7", "GAS6", "SLC25A13",
                  "ITGB1", "SPP1", "LASP1")

all_violins <- featDistPlot(seurat_data, violin_genes,
                            sep_by = "sample", 
                            color = sample_colors,
                            combine = FALSE)

invisible(lapply(names(all_violins), function(x){
  violin_plot <- all_violins[[x]] +
    ggplot2::ylab("") +
    ggplot2::ggtitle(x)
  
  pdf(file.path(fig_dir, paste0(x, "_violin.pdf")))
  print(violin)
}))

## Pseudotime plots ------------------------------------------------------------


# Supplement -------------------------------------------------------------------

## UMAPs -----------------------------------------------------------------------
# Separated by genotype, colored by celltype

all_plots <- lapply(unique(seurat_data$genotype), function(x){
  plotDimRed(seurat_data, col_by = "celltype", color = all_colors2,
             highlight_group = TRUE, group = x, 
             meta_data_col = "genotype",
             plot_type = "rna.umap", ggrastr = TRUE)[[1]]
})

names(all_plots) <- unique(seurat_data$genotype)

all_cell_types <- cowplot::plot_grid(all_plots$WT,
                                     all_plots$CFKO,
                                     nrow = 1, ncol = 2)

pdf(file.path(fig_dir, "genotype_combined_cell_type_umap.pdf"),
    width = 10, height = 5)
print(all_cell_types)

dev.off()
## Heatmaps of gene lists ------------------------------------------------------
# Plot heatmap across all cell types in only day 9
# Keep working here:
# 1. Remove celltypes that are only in 1 sample
# 2. Order based on the order of colors2
# Make for all gene lists

gene_path <- here("files/GSEA_signaling_pathways_with_orthologs.xlsx")

all_sheets <- openxlsx::getSheetNames(gene_path)

gene_lists <- lapply(all_sheets, function(x){
  gene_df <- openxlsx::readWorkbook(gene_path, sheet = x)
  return(unique(gene_df$gene_id))
})

names(gene_lists) <- all_sheets

# Subset seurat data
seurat_sub <- subset(seurat_data, subset = day == "D9")
seurat_sub$celltype_sample <- paste(seurat_sub$celltype,
                                    seurat_sub$sample, 
                                    sep = "_")

# Figure out what cell types to plot
cell_types <- seurat_sub[[]] %>%
  dplyr::select(sample, celltype) %>%
  dplyr::distinct() %>%
  dplyr::group_by(celltype) %>%
  dplyr::add_count(name = "sample_count") %>%
  dplyr::filter(sample_count > 1)

# Subset seurat object to just those cell types
seurat_sub <- subset(seurat_sub, subset = celltype %in%
                       unique(cell_types$celltype))

# Make meta df
meta_df <- seurat_sub[[c("sample", "celltype",
                          "celltype_sample")]]
meta_ave <- data.frame(table(paste0(seurat_sub$sample, ";",
                                    seurat_sub$celltype)))

meta_ave$celltype <- gsub(".*_D[0-9]+;", "", 
                                       meta_ave$Var1)
meta_ave$sample <- gsub(";.*$", "", meta_ave$Var1)
meta_ave$Var1 <- NULL
meta_ave$Freq <- NULL
meta_ave$celltype_sample <- paste(meta_ave$celltype,
                                  meta_ave$sample,
                                  sep = "_")

# Get cell type order
meta_ave$sample <- factor(meta_ave$sample,
                          levels = names(sample_colors))

meta_ave$celltype <- factor(meta_ave$celltype,
                            levels = names(all_colors2))

meta_ave <- meta_ave %>%
  dplyr::arrange(celltype, sample)

meta_ave$celltype_sample <- factor(meta_ave$celltype_sample,
                                   levels = meta_ave$celltype_sample)

rownames(meta_ave) <- meta_ave$celltype_sample
keep_colors_sample <- sample_colors[names(sample_colors) %in%
                                      unique(meta_ave$sample)]

keep_colors_celltype <- all_colors2[names(all_colors2) %in%
                                      unique(meta_ave$celltype)]

color_list <- list("sample" = keep_colors_sample,
                   "celltype" = keep_colors_celltype)

invisible(lapply(names(gene_lists), function(x){
  print(x)
  pdf(file.path(fig_dir,
                paste0(x, "_heatmap_average.pdf")))
  print(plot_heatmap(seurat_sub, gene_list = gene_lists[[x]],
                     meta_col = "celltype_sample", 
                     average_expression = TRUE,
                     colors = sample_colors,
                     plot_rownames = FALSE,
                     meta_df = meta_ave, 
                     color_list = color_list,
                     plot_meta_col = FALSE,
                     cluster_rows = TRUE))
  
  dev.off()
  
  fig_height <- round(length(gene_lists[[x]]) / 10)

  pdf(file.path(fig_dir,
                paste0(x, "_heatmap_average_labeled.pdf")),
      width = 10, height = fig_height)
  print(plot_heatmap(seurat_sub, gene_list = gene_lists[[x]],
                     meta_col = "celltype_sample", 
                     average_expression = TRUE,
                     colors = sample_colors,
                     plot_rownames = TRUE,
                     meta_df = meta_ave, 
                     color_list = color_list,
                     plot_meta_col = FALSE,
                     cluster_rows = TRUE))
  
  dev.off()
}))

