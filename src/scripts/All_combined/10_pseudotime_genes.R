library(slingshot)
library(scAnalysisR)
library(Seurat)
library(here)
library(tidyverse)
library(tradeSeq)
library(pheatmap)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

all_samples <- "All_combined"

all_sample_dir <- here("results", all_samples, "R_analysis")

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

sample_colors <- as.character(LaCroixColoR::lacroix_palette("Coconut", 10))
sample_colors[5] <- "#F4E3C7"
names(sample_colors) <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                          "CFKO_D14", "CFKO_D9", "CFKO_D7", "CFKO_D5", "CFKO_D2")
sample_colors <- sample_colors[c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                                 "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")]


# Functions --------------------------------------------------------------------
# Seurat object = seurat object
# lineages = list of lineages that are in the metadata of the seurat object
# Gene = gene to plot across pseudotime
# col_by = color for the plot (can be "lineage" to color by lineage,
# "pseudotime" to color by the pseudotime vals or
# a value from the metadata)
# color = custom colors to pass
plotPseudotime <- function(seurat_object, lineages, gene_list,
                           col_by = "lineage", color = NULL,
                           save_plot = NULL, ...){
  
  plot_list <- lapply(gene_list, function(x) {
    print(x)
    plotPseudotimeSingle(seurat_object = seurat_object,
                         lineages = lineages, gene = x,
                         col_by = col_by, color = color, ...)
  })
  if (!is.null(save_plot)){
    pdf(save_plot)
    print(plot_list)
    dev.off()
  }
  return(plot_list)
}
 
plotPseudotimeSingle <- function(seurat_object, lineages, gene,
                                 col_by = "lineage", color = NULL,
                                 max_pseudotime = NULL, ...){
  
  # Get data
  # If col_by is lineage or pseudotime, do this,
  # otherwise add col_by to the call
  if (col_by %in% c("pseudotime", "lineage")){
    data_info <- Seurat::FetchData(seurat_object, vars = c(lineages, gene))
  } else {
    data_info <- Seurat::FetchData(seurat_object,
                                   vars = c(lineages, gene, col_by))
  }
  
  # Here I'll need to pivot_longer by the lineages
  data_info <- data_info %>%
    tidyr::pivot_longer(cols = all_of(lineages), names_to = "lineage",
                 values_to = "pseudotime") %>%
    dplyr::rename(gene = dplyr::all_of(gene)) %>%
    dplyr::filter(!is.na(pseudotime))

  data_info$col_by <- data_info[[col_by]]
  
  if(!is.null(max_pseudotime)){
    data_info <- data_info %>%
      dplyr::filter(pseudotime < max_pseudotime)
  }

  # Make plot
  base_plot <- ggplot2::ggplot(data = data_info, ggplot2::aes(x = pseudotime,
                                                              y = gene,
                                                              color = col_by,
                                                              group = lineage)) +
    ggplot2::geom_point(alpha = 0.75, size = 0.25) +
    ggplot2::geom_smooth(method = "gam", color = "black") +
    ggplot2::ylab(gene)
  
  # Next figure out colors, if col_by is continuous, default is 
  # viridis. If it's not continuous do RColorBrewer as I like
  if (is.numeric(data_info$col_by)){
    if (is.null(color)) {
      base_plot <- base_plot +
        ggplot2::scale_color_viridis_c(option = "magma") +
        labs(color = col_by)
    } else {
      low <- color[1]
      high <- color[2]
      base_plot <- base_plot + 
        ggplot2::scale_color_gradient(low = low, high = high, name = col_by)
    }
  } else {
    # Make colors if they don't already exist
    if (is.null(color)){
      # Factor in case it isn't already
      if(!is.factor(data_info$col_by)){
        data_info$col_by <- factor(data_info$col_by)
      }
      nColors <- length(levels(data_info$col_by))
      color <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors)
      names(color) <- levels(data_info$col_by)
    }
    base_plot <- base_plot + scale_color_manual(values = color,
                                                name = col_by)
  }
  
  return(base_plot)

}


make_lineage_heatmaps <- function(seurat_object, lineage_name,
                                  lineage_number, lineage_asso_res,
                                  image_dir, sample_colors,
                                  waldcutoff = 200, pval_cutoff = 0.05){
  
  keep_cols <- c("waldStat", "pvalue", "df")
  keep_cols <- paste0(keep_cols, "_", lineage_number)
  
  # Filter by pvalue
  filtered_res <- lineage_asso_res %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    setNames(gsub("_[0-9]*", "", names(.))) %>%
    dplyr::filter(pvalue < pval_cutoff)
  
  # Filter by wald cutoff (fewer genes)
  filtered_res_high <- filtered_res %>%
    dplyr::filter(waldStat > waldcutoff)
  
  # Subset to just the lineage
  cells <- seurat_object[[]] %>%
    dplyr::select(dplyr::all_of(lineage_name)) %>%
    tidyr::drop_na() %>%
    rownames()
  
  # Keep only cells in the lineage
  seurat_sub <- subset(seurat_object, cells = cells)
  
  # Order the cells by pseudotime
  order_cells <- seurat_sub[[]] %>%
    dplyr::select(all_of(lineage_name)) %>%
    dplyr::arrange_all()
  
  # Make plot of all sig genes
  pdf(file.path(image_dir, paste0(lineage_name, "_heatmap.pdf")))
  print(plot_heatmap(seurat_sub, gene_list = rownames(filtered_res),
                     meta_col = "sample", colors = sample_colors,
                     cell_order = rownames(order_cells),
                     cluster_rows = TRUE, plot_rownames = FALSE))
  dev.off()
  
  
  
  # Make plot of all genes with high wald
  pdf(file.path(image_dir,
                paste0(lineage_name, "_heatmap_wald", waldcutoff, ".pdf")))
  
  print(plot_heatmap(seurat_sub, gene_list = rownames(filtered_res_high),
                     meta_col = "sample", colors = sample_colors,
                     cell_order = rownames(order_cells),
                     cluster_rows = TRUE, plot_rownames = FALSE))
  dev.off()
  
  
  # Make plot of all genes with high wald as png
  png(file.path(image_dir,
                paste0(lineage_name, "_heatmap_wald", waldcutoff, ".png")))
  
  print(plot_heatmap(seurat_sub, gene_list = rownames(filtered_res_high),
                     meta_col = "sample", colors = sample_colors,
                     cell_order = rownames(order_cells),
                     cluster_rows = TRUE, plot_rownames = FALSE))
  dev.off()
  
  
  # Make a heatmap with rownames, return so that the size can be altered
  return_heatmap <- plot_heatmap(seurat_sub, 
                                 gene_list = rownames(filtered_res_high),
                                 meta_col = "sample", colors = sample_colors,
                                 cell_order = rownames(order_cells),
                                 cluster_rows = TRUE, plot_rownames = TRUE)
  
  return(list(sig_genes = filtered_res,
              sig_high_genes = filtered_res_high,
              heatmap = return_heatmap))
  
}

# Analysis ---------------------------------------------------------------------

# cfko slingshot
cfko_slingshot <- readRDS(file.path(all_sample_dir, "files",
                                    "slingshot", "CFKO_res.rds"))


cfko_pseudotime <- slingPseudotime(cfko_slingshot)

cfko_fitgam <- readRDS(file.path(all_sample_dir, "files", "slingshot", 
                                 "CFKO_fitgam_sce.rda"))

# Find converged genes
table(rowData(cfko_fitgam)$tradeSeq$converged)

# Testing weather the average gene expression is significantly changing
# across pseudotime
assoRes <- associationTest(cfko_fitgam)

lineage_assoRes <- associationTest(cfko_fitgam, lineages = TRUE)

write.csv(lineage_assoRes, file.path(all_sample_dir, "files",
                                             "slingshot",
                                             "cfko_association_test.csv"))

# lineage_2_genes <- lineage_assoRes %>%
#   dplyr::filter(!is.na(pvalue_2)) %>%
#   dplyr::filter(pvalue_2 < 0.05)
# 
# 
# lineage_7_genes <- lineage_assoRes %>%
#   dplyr::filter(!is.na(pvalue_7)) %>%
#   dplyr::filter(pvalue_7 < 0.05)

image_dir <- file.path(all_sample_dir, "images", "slingshot")

# Lineage 2 heatmaps
lineage_2_info <- make_lineage_heatmaps(seurat_object = merged_seurat,
                                        lineage_name = "cfko_Lineage2",
                                        lineage_number = "2",
                                        lineage_asso_res = lineage_assoRes,
                                        image_dir = image_dir,
                                        sample_colors = sample_colors,
                                        waldcutoff = 200, pval_cutoff = 0.05)

png(file.path(all_sample_dir, "images", "slingshot",
              "cfko_lineage2_heatmap_wald200_names.png"),
    width = 900, height = 1300)

print(lineage_2_info$heatmap)

dev.off()


# Lineage 7 heatmaps
lineage_7_info <- make_lineage_heatmaps(seurat_object = merged_seurat,
                                        lineage_name = "cfko_Lineage7",
                                        lineage_number = "7",
                                        lineage_asso_res = lineage_assoRes,
                                        image_dir = image_dir,
                                        sample_colors = sample_colors,
                                        waldcutoff = 200, pval_cutoff = 0.05)

png(file.path(all_sample_dir, "images", "slingshot",
              "cfko_lineage7_heatmap_wald200_names.png"),
    width = 900, height = 900)

print(lineage_7_info$heatmap)

dev.off()

# Lineage 8 heatmaps
lineage_8_info <- make_lineage_heatmaps(seurat_object = merged_seurat,
                                        lineage_name = "cfko_Lineage8",
                                        lineage_number = "8",
                                        lineage_asso_res = lineage_assoRes,
                                        image_dir = image_dir,
                                        sample_colors = sample_colors,
                                        waldcutoff = 200, pval_cutoff = 0.05)

png(file.path(all_sample_dir, "images", "slingshot",
              "cfko_lineage8_heatmap_wald200_names.png"),
    width = 900, height = 400)

print(lineage_8_info$heatmap)

dev.off()

# Lineage 5 heatmaps
lineage_5_info <- make_lineage_heatmaps(seurat_object = merged_seurat,
                                        lineage_name = "cfko_Lineage5",
                                        lineage_number = "5",
                                        lineage_asso_res = lineage_assoRes,
                                        image_dir = image_dir,
                                        sample_colors = sample_colors,
                                        waldcutoff = 200, pval_cutoff = 0.05)

png(file.path(all_sample_dir, "images", "slingshot",
              "cfko_lineage5_heatmap_wald200_names.png"),
    width = 900, height = 1300)

print(lineage_7_info$heatmap)

dev.off()

colors <- LaCroixColoR::lacroix_palette(name = "Pamplemousse", n = 2)

names(colors) <- c("cfko_Lineage2", "cfko_Lineage5")

pdf(file.path(image_dir, "example_pseudotime_plots.pdf"))

pseudotime_plots1 <- plotPseudotime(seurat_object,
                                   lineages = c("cfko_Lineage5",
                                                "cfko_Lineage2"),
                                   gene_list = c("MMP13", "ECRG4", "CLU"),
                                   col_by = "lineage", color = colors,
                                   max_pseudotime = 100)

pseudotime_plots2 <- plotPseudotime(seurat_object,
                                   lineages = c("cfko_Lineage5"),
                                   gene_list = c("MMP13", "ECRG4", "CLU"),
                                   col_by = "lineage", color = colors,
                                   max_pseudotime = 100)

pseudotime_plots3 <- plotPseudotime(seurat_object,
                                   lineages = c("cfko_Lineage2"),
                                   gene_list = c("FABP5", "ETS2", "CD46"),
                                   col_by = "lineage", color = colors,
                                   max_pseudotime = 100)

print(pseudotime_plots1)
print(pseudotime_plots2)
print(pseudotime_plots3)

dev.off()

# Testing difference in expression between two end points
diffEnd <- diffEndTest(cfko_fitgam)

pseudotime_plots <- plotPseudotime(seurat_object,
                                   lineages = c("cfko_Lineage2",
                                                "cfko_Lineage7"),
                                   gene_list = c("UPK1A", "LY6D"),
                                   col_by = "lineage", color = NULL)


pseudotime_plot2 <- plotPseudotime(seurat_object,
                                   lineages = c("cfko_Lineage2"),
                                   gene_list = c("UPK1A", "LY6D"),
                                   col_by = "pseudotime", color = NULL)


pseudotime_plot3 <- plotPseudotime(seurat_object,
                                   lineages = c("cfko_Lineage2"),
                                   gene_list = c("UPK1A", "LY6D"),
                                   col_by = "lineage", color = NULL)


pseudotime_plot4 <- plotPseudotime(seurat_object,
                                   lineages = c("cfko_Lineage2"),
                                   gene_list = c("UPK1A", "LY6D"),
                                   col_by = "sample", color = sample_colors)
