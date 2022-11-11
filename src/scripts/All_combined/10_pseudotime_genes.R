library(slingshot)
library(scAnalysisR)
library(Seurat)
library(here)
library(tidyverse)
library(tradeSeq)
library(pheatmap)
library(MetBrewer)

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
# Max pseudotime = for aesthetics - should a few cells that are outliers 
# be removed? This should be improved
plotPseudotime <- function(seurat_object, lineages, gene_list,
                           col_by = "lineage", color = NULL,
                           save_plot = NULL, max_pseudotime = NULL, ...){
  
  plot_list <- lapply(gene_list, function(x) {
    print(x)
    plotPseudotimeSingle(seurat_object = seurat_object,
                         lineages = lineages, gene = x,
                         col_by = col_by, color = color,
                         max_pseudotime = max_pseudotime, ...)
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
                                 max_pseudotime = NULL,
                                 line_color = "black", alpha = 0.75,
                                 ...){
  
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
    ggplot2::geom_point(alpha = alpha, size = 0.25) +
    ggplot2::ylab(gene)
  
  if(line_color == "lineage"){
    base_plot <- base_plot +
      ggplot2::geom_smooth(method = "gam")
  } else {
    base_plot <- base_plot + 
      ggplot2::geom_smooth(method = "gam", color = line_color)
  }
  
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

plot_overlapping_lineages <- function(seurat_object,
                                      cfko_lineage_genes,
                                      wt_lineage_genes,
                                      cfko_lineage_name,
                                      wt_lineage_name,
                                      save_dir = NULL,
                                      colors = NULL){
  
  if(is.null(colors)){
    lineage_colors <- met.brewer("Egypt", 2)
    names(lineage_colors) <- c(wt_lineage_name, cfko_lineage_name)
    
  }
  cfko_lineage_df <- data.frame("cfko_cluster" = cfko_lineage_genes,
                                 "gene" = names(cfko_lineage_genes))
  
  wt_lineage_df <- data.frame("wt_cluster" = wt_lineage_genes,
                               "gene" = names(wt_lineage_genes))
  
  combined_lineage <- merge(cfko_lineage_df, wt_lineage_df,
                            by = "gene", all.x = FALSE, all.y = FALSE)
  
  
  if(!is.null(save_dir)){
    all_plots <- plotPseudotime(seurat_object,
                                lineages = c(cfko_lineage_name,
                                             wt_lineage_name),
                                gene_list = combined_lineage$gene,
                                col_by = "lineage", color = lineage_colors,
                                line_color = "lineage",
                                alpha = 0.25)
    
    
    pdf(file.path(save_dir,
                  paste0(wt_lineage_name, "_and_", cfko_lineage_name, ".pdf")),
        width = 4, height = 4)
    
    print(all_plots)
    
    dev.off()    
  }

  
  return(combined_lineage)
  
}

# Analysis ---------------------------------------------------------------------


## CFKO ------------------------------------------------------------------------
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

# Get genes in each cluster
gene_info_cfko_lineage_2 <- sort(cutree(lineage_2_info$heatmap$tree_row, k=4))

# Cluster mapping
# 1 = high late
# 2 = high early
# 3 = high for 50% low
# 4 = Low high low

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

gene_info_cfko_lineage_7 <- sort(cutree(lineage_7_info$heatmap$tree_row, k=3))

# Cluster mapping
# 1 = high late
# 2 = high early
# 3 = random

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

gene_info_cfko_lineage_8 <- sort(cutree(lineage_8_info$heatmap$tree_row, k=2))

# Cluster mapping
# 1 = high late end
# 2 = high late longer

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

print(lineage_5_info$heatmap)

dev.off()

gene_info_cfko_lineage_5 <- sort(cutree(lineage_5_info$heatmap$tree_row, k=5))

# Cluster mapping
# 1 - Early, slightly lower and spottier
# 2 - Early more consistent
# 3 - End, less than 50%
# 4 - Early - first half
# 5 - End - Second half

colors <- LaCroixColoR::lacroix_palette(name = "Pamplemousse", n = 2)

names(colors) <- c("cfko_Lineage2", "cfko_Lineage5")

pdf(file.path(image_dir, "example_pseudotime_plots.pdf"))

pseudotime_plots1 <- plotPseudotime(merged_seurat,
                                   lineages = c("cfko_Lineage5",
                                                "cfko_Lineage2"),
                                   gene_list = c("MMP13", "ECRG4", "CLU"),
                                   col_by = "lineage", color = colors,
                                   max_pseudotime = 100)

pseudotime_plots2 <- plotPseudotime(merged_seurat,
                                   lineages = c("cfko_Lineage5"),
                                   gene_list = c("MMP13", "ECRG4", "CLU"),
                                   col_by = "lineage", color = colors,
                                   max_pseudotime = 100)

pseudotime_plots3 <- plotPseudotime(merged_seurat,
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

pseudotime_plots <- plotPseudotime(merged_seurat,
                                   lineages = c("cfko_Lineage2",
                                                "cfko_Lineage7"),
                                   gene_list = c("UPK1A", "LY6D"),
                                   col_by = "lineage", color = NULL)


pseudotime_plot2 <- plotPseudotime(merged_seurat,
                                   lineages = c("cfko_Lineage2"),
                                   gene_list = c("UPK1A", "LY6D"),
                                   col_by = "pseudotime", color = NULL)


pseudotime_plot3 <- plotPseudotime(merged_seurat,
                                   lineages = c("cfko_Lineage2"),
                                   gene_list = c("UPK1A", "LY6D"),
                                   col_by = "lineage", color = NULL)


pseudotime_plot4 <- plotPseudotime(merged_seurat,
                                   lineages = c("cfko_Lineage2"),
                                   gene_list = c("UPK1A", "LY6D"),
                                   col_by = "sample", color = sample_colors)


## WT --------------------------------------------------------------------------

# cfko slingshot
wt_slingshot <- readRDS(file.path(all_sample_dir, "files",
                                    "slingshot", "WT_res.rds"))


wt_pseudotime <- slingPseudotime(wt_slingshot)

wt_fitgam <- readRDS(file.path(all_sample_dir, "files", "slingshot", 
                                 "WT_fitgam_sce.rda"))

# Find converged genes
table(rowData(wt_fitgam)$tradeSeq$converged)

# Testing weather the average gene expression is significantly changing
# across pseudotime
assoRes <- associationTest(wt_fitgam)

lineage_assoRes <- associationTest(wt_fitgam, lineages = TRUE)

write.csv(lineage_assoRes, file.path(all_sample_dir, "files",
                                     "slingshot",
                                     "wt_association_test.csv"))

# I like 3, 4, 5, 7 (1, 2?)

image_dir <- file.path(all_sample_dir, "images", "slingshot")

# Lineage 1 heatmaps
set.seed(100)
lineage_1_info <- make_lineage_heatmaps(seurat_object = merged_seurat,
                                        lineage_name = "wt_Lineage1",
                                        lineage_number = "1",
                                        lineage_asso_res = lineage_assoRes,
                                        image_dir = image_dir,
                                        sample_colors = sample_colors,
                                        waldcutoff = 200, pval_cutoff = 0.05)


png(file.path(all_sample_dir, "images", "slingshot",
              "wt_lineage1_heatmap_wald200_names.png"),
    width = 900, height = 2000)

print(lineage_1_info$heatmap)

dev.off()

gene_info_wt_lineage_1 <- sort(cutree(lineage_1_info$heatmap$tree_row, k=10))

plot_genes <- lapply(1:10, function(x){
  return_genes <- gene_info_wt_lineage_1[gene_info_wt_lineage_1 == x][1:2]
  return(return_genes)
})

plot_genes <- unlist(plot_genes)

pseudotime_plot <- plotPseudotime(merged_seurat,
                                  lineages = c("wt_Lineage1"),
                                  gene_list = names(plot_genes),
                                  col_by = "lineage", color = NULL)


# Cluster mapping
# 1 - High low high low
# 2 - Low High, low
# 3 - Low high low - sharper
# 4 - High low
# 5 - high low high
# 6 - High low (slightly low to start)
# 7 - low high 
# 8 - low high low
# 9 - low high low (less stable)
# 10 - low high

# Lineage 3 heatmaps
set.seed(100)
lineage_3_info <- make_lineage_heatmaps(seurat_object = merged_seurat,
                                        lineage_name = "wt_Lineage3",
                                        lineage_number = "3",
                                        lineage_asso_res = lineage_assoRes,
                                        image_dir = image_dir,
                                        sample_colors = sample_colors,
                                        waldcutoff = 200, pval_cutoff = 0.05)

callback <- function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}


png(file.path(all_sample_dir, "images", "slingshot",
              "wt_lineage3_heatmap_wald200_names.png"),
    width = 900, height = 2000)

print(lineage_3_info$heatmap)

dev.off()

gene_info_wt_lineage_3 <- sort(cutree(lineage_3_info$heatmap$tree_row, k=5))

pseudotime_plot <- plotPseudotime(merged_seurat,
                                   lineages = c("wt_Lineage3"),
                                   gene_list = c("BMP3", "LTF",
                                                 "PROM1", "NT5DC2",
                                                 "SLC34A2", "LIMCH1",
                                                 "INKA1", "IMPDH2",
                                                 "ITPK1", "MALL"),
                                   col_by = "lineage", color = NULL)

# Cluster mapping
# 1 - High late, no expression early
# 2 - Low High, low
# 3 - Low high low - slightly later
# 4 - High low
# 5 - End - second half

# Lineage 4 heatmaps
lineage_4_info <- make_lineage_heatmaps(seurat_object = merged_seurat,
                                        lineage_name = "wt_Lineage4",
                                        lineage_number = "4",
                                        lineage_asso_res = lineage_assoRes,
                                        image_dir = image_dir,
                                        sample_colors = sample_colors,
                                        waldcutoff = 200, pval_cutoff = 0.05)

png(file.path(all_sample_dir, "images", "slingshot",
              "wt_lineage4_heatmap_wald200_names.png"),
    width = 900, height = 1300)

print(lineage_4_info$heatmap)

dev.off()

gene_info_wt_lineage_4 <- sort(cutree(lineage_4_info$heatmap$tree_row, k=5))

pseudotime_plot <- plotPseudotime(merged_seurat,
                                  lineages = c("wt_Lineage4"),
                                  gene_list = c("CEBPD", "RGS12",
                                                "TXN", "DCXR",
                                                "ICAM1", "MUC16",
                                                "LAMB3", "AHNAK",
                                                "MGST3", "EIF4EBP1"),
                                  col_by = "lineage", color = NULL)

# Cluster mapping
# 1 - High to low
# 2 - Low High, low
# 3 - Low high 
# 4 - High low high
# 5 - low high low

# Lineage 5 heatmaps
lineage_5_info <- make_lineage_heatmaps(seurat_object = merged_seurat,
                                        lineage_name = "wt_Lineage5",
                                        lineage_number = "5",
                                        lineage_asso_res = lineage_assoRes,
                                        image_dir = image_dir,
                                        sample_colors = sample_colors,
                                        waldcutoff = 200, pval_cutoff = 0.05)

png(file.path(all_sample_dir, "images", "slingshot",
              "wt_lineage5_heatmap_wald200_names.png"),
    width = 900, height = 1600)

print(lineage_5_info$heatmap)

dev.off()

gene_info_wt_lineage_5 <- sort(cutree(lineage_5_info$heatmap$tree_row, k=8))


# TODO
# pseudotime_plot <- plotPseudotime(merged_seurat,
#                                   lineages = c("wt_Lineage4"),
#                                   gene_list = c("CEBPD", "RGS12",
#                                                 "TXN", "DCXR",
#                                                 "ICAM1", "MUC16",
#                                                 "LAMB3", "AHNAK",
#                                                 "MGST3", "EIF4EBP1"),
#                                   col_by = "lineage", color = NULL)

# Cluster mapping



# Lineage 7 heatmaps
lineage_7_info <- make_lineage_heatmaps(seurat_object = merged_seurat,
                                        lineage_name = "wt_Lineage7",
                                        lineage_number = "7",
                                        lineage_asso_res = lineage_assoRes,
                                        image_dir = image_dir,
                                        sample_colors = sample_colors,
                                        waldcutoff = 200, pval_cutoff = 0.05)

png(file.path(all_sample_dir, "images", "slingshot",
              "wt_lineage7_heatmap_wald200_names.png"),
    width = 900, height = 1600)

print(lineage_7_info$heatmap)

dev.off()

gene_info_wt_lineage_7 <- sort(cutree(lineage_7_info$heatmap$tree_row, k=4))

# TODO
# pseudotime_plot <- plotPseudotime(merged_seurat,
#                                   lineages = c("wt_Lineage4"),
#                                   gene_list = c("CEBPD", "RGS12",
#                                                 "TXN", "DCXR",
#                                                 "ICAM1", "MUC16",
#                                                 "LAMB3", "AHNAK",
#                                                 "MGST3", "EIF4EBP1"),
#                                   col_by = "lineage", color = NULL)

# Cluster mapping

# Plots of lineages alone ------------------------------------------------------
lineage_colors_base <- met.brewer("Egypt", 2)

# Plots of lineages alone
pdf(file.path(all_sample_dir, "images", "slingshot",
              "expression_over_time",
              "wt_lineage1.pdf"),
    width = 4, height = 4)

lineage_colors <- lineage_colors_base[1]
names(lineage_colors) <- "wt_Lineage1"

all_plots <- plotPseudotime(merged_seurat,
                            lineages = c("wt_Lineage1"),
                            gene_list = names(gene_info_wt_lineage_1),
                            col_by = "lineage", color = lineage_colors,
                            line_color = "lineage",
                            alpha = 0.25)

print(all_plots)

dev.off()

pdf(file.path(all_sample_dir, "images", "slingshot",
              "expression_over_time",
              "wt_lineage4.pdf"),
    width = 4, height = 4)

lineage_colors <- lineage_colors_base[1]
names(lineage_colors) <- "wt_Lineage4"

all_plots <- plotPseudotime(merged_seurat,
                            lineages = c("wt_Lineage4"),
                            gene_list = names(gene_info_wt_lineage_4),
                            col_by = "lineage", color = lineage_colors,
                            line_color = "lineage",
                            alpha = 0.25)

print(all_plots)

dev.off()

lineage_colors <- lineage_colors_base[2]
names(lineage_colors) <- "cfko_Lineage2"
pdf(file.path(all_sample_dir, "images", "slingshot",
              "expression_over_time",
              "cfko_lineage2.pdf"),
    width = 4, height = 4)

all_plots <- plotPseudotime(merged_seurat,
                            lineages = c("cfko_Lineage2"),
                            gene_list = names(gene_info_cfko_lineage_2),
                            col_by = "lineage", color = lineage_colors,
                            line_color = "lineage",
                            alpha = 0.25)

print(all_plots)

dev.off()


# Comparing WT and CF ----------------------------------------------------------

# WT lineage 1
# Cluster mapping
# 1 - High low high low
# 2 - Low High, low
# 3 - Low high low - sharper
# 4 - High low
# 5 - high low high
# 6 - High low (slightly low to start)
# 7 - low high 
# 8 - low high low
# 9 - low high low (less stable)
# 10 - low high

# WT lineage 4
# Cluster mapping
# 1 - High to low
# 2 - Low High, low
# 3 - Low high 
# 4 - High low high
# 5 - low high low

# CFKO lineage 2
# Cluster mapping
# 1 = low high
# 2 = high low
# 3 = high for 50% low
# 4 = Low high low

# Lineage colors
lineage_colors_base <- met.brewer("Egypt", 2)

# CFKO lineage 2 and WT lineage 4
combined_lineage <- 
  plot_overlapping_lineages(seurat_object = merged_seurat,
                            cfko_lineage_genes = gene_info_cfko_lineage_2,
                            wt_lineage_genes = gene_info_wt_lineage_4,
                            cfko_lineage_name = "cfko_Lineage2",
                            wt_lineage_name = "wt_Lineage4",
                            save_dir = file.path(all_sample_dir, "images", 
                                                 "slingshot",
                                                 "expression_over_time"))



# 1 cfko + 4 wt look good
# 4 cfko + 3 wt look good
# 1 cfko + 2 wt look good

genes_to_plot <- c("LOC123388619", "LOC123388635",
                   "MUC4", "TXN", "UGGT1")


lineage_colors <- lineage_colors_base
names(lineage_colors) <- c("wt_Lineage4", "cfko_Lineage2")

different_plots <- plotPseudotime(merged_seurat,
                                   lineages = c("cfko_Lineage2",
                                                "wt_Lineage4"),
                                   gene_list = genes_to_plot,
                                   col_by = "lineage", color = lineage_colors,
                                  line_color = "lineage",
                                  alpha = 0.25)


# CFKO lineage 2 and WT lineage 1
combined_lineage <- 
  plot_overlapping_lineages(seurat_object = merged_seurat,
                            cfko_lineage_genes = gene_info_cfko_lineage_2,
                            wt_lineage_genes = gene_info_wt_lineage_1,
                            cfko_lineage_name = "cfko_Lineage2",
                            wt_lineage_name = "wt_Lineage1",
                            save_dir = file.path(all_sample_dir, "images", 
                                                 "slingshot",
                                                 "expression_over_time"))
