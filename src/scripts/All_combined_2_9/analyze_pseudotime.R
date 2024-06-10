library(slingshot)
library(scAnalysisR)
library(Seurat)
library(here)
library(tidyverse)
library(tradeSeq)

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
    data_info <- Seurat::FetchData(seurat_object, vars = c(lineages, gene),
                                   slot = "data")
  } else {
    data_info <- Seurat::FetchData(seurat_object,
                                   vars = c(lineages, gene, col_by),
                                   slot = "data")
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
                                      colors = NULL,
                                      type = "intersection"){
  
  if(is.null(colors)){
    lineage_colors <- met.brewer("Egypt", 2)
    names(lineage_colors) <- c(wt_lineage_name, cfko_lineage_name)
    
  }
  
  if(type == "intersection"){
    plot_genes <- intersect(names(cfko_lineage_genes),
                            names(wt_lineage_genes))
  } else if(type =="difference") {
    plot_genes <- unique(c(setdiff(names(cfko_lineage_genes),
                                   names(wt_lineage_genes)),
                           setdiff(names(wt_lineage_genes),
                                   names(cfko_lineage_genes))))
  }
  
  cfko_lineage_df <- data.frame("cfko_cluster" = cfko_lineage_genes,
                                "gene" = names(cfko_lineage_genes))
  
  wt_lineage_df <- data.frame("wt_cluster" = wt_lineage_genes,
                              "gene" = names(wt_lineage_genes))
  
  if(type == "intersection"){
    combined_lineage <- merge(cfko_lineage_df, wt_lineage_df,
                              by = "gene", all.x = FALSE, all.y = FALSE)    
  } else if(type == "difference"){
    combined_lineage <- merge(cfko_lineage_df, wt_lineage_df,
                              by = "gene", all.x = TRUE, all.y = FALSE)
    
    combined_lineage <- combined_lineage[combined_lineage$gene %in% plot_genes,]
  }
  
  
  
  if(!is.null(save_dir)){
    all_plots <- plotPseudotime(seurat_object,
                                lineages = c(cfko_lineage_name,
                                             wt_lineage_name),
                                gene_list = plot_genes,
                                col_by = "lineage", color = lineage_colors,
                                line_color = "lineage",
                                alpha = 0.25)
    
    
    pdf(file.path(save_dir,
                  paste0(wt_lineage_name, "_and_", cfko_lineage_name,
                         "_", type, ".pdf")),
        width = 4, height = 4)
    
    print(all_plots)
    
    dev.off()    
  }
  
  
  return(combined_lineage)
  
}

# Analysis ---------------------------------------------------------------------
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

all_samples <- "All_combined_2_9"

all_sample_dir <- here("results", all_samples, "R_analysis")

revision_dir <- here("results", all_samples, "R_analysis", 
                     "images", "revision")

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

all_colors <- c("acinar" = "#D4405B",
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


# cfko slingshot
cfko_slingshot <- readRDS(file.path(revision_dir, "CFKO_slingshot.rds"))

cfko_pseudotime <- slingPseudotime(cfko_slingshot)
colnames(cfko_pseudotime) <- paste0("cfko_", colnames(cfko_pseudotime))

cfko_pseudotime <- data.frame(cfko_pseudotime)

# wt slingshot
wt_slingshot <- readRDS(file.path(revision_dir, "WT_slingshot.rds"))

wt_pseudotime <- slingPseudotime(wt_slingshot)
colnames(wt_pseudotime) <- paste0("wt_", colnames(wt_pseudotime))

wt_pseudotime <- data.frame(wt_pseudotime)


# Plot on UMAP
all_umap <- Embeddings(merged_seurat, reduction = "rna.umap")
cfko_umap <- all_umap[rownames(all_umap) %in% rownames(cfko_pseudotime),]

cfko_umap_slingshot <- embedCurves(cfko_slingshot, newDimRed = cfko_umap)


all_curves <- slingCurves(cfko_umap_slingshot)

all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  
  # This line cuts off the long tail... Probably a better option is
  # to just remove unknown before running slingshot.
  base_plot <- plotDimRed(merged_seurat, col_by = "sample", color = sample_colors,
                          plot_type = "rna.umap",
                          ggrastr = TRUE)[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)  +
    ggplot2::ggtitle(paste0("cfko_lineage", x))
  
  return(base_plot)
})

all_umap <- Embeddings(merged_seurat, reduction = "rna.umap")
wt_umap <- all_umap[rownames(all_umap) %in% rownames(wt_pseudotime),]

wt_umap_slingshot <- embedCurves(wt_slingshot, newDimRed = wt_umap)


all_curves <- slingCurves(wt_umap_slingshot)

all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  
  # This line cuts off the long tail... Probably a better option is
  # to just remove unknown before running slingshot.
  base_plot <- plotDimRed(merged_seurat, col_by = "sample", color = sample_colors,
                          plot_type = "rna.umap",
                          ggrastr = TRUE)[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)  +
    ggplot2::ggtitle(paste0("wt_lineage", x))
  
  return(base_plot)
})

# Lineage 2 for CF lineage 7 for WT
colnames(cfko_pseudotime) <- paste0("corrected_", colnames(cfko_pseudotime))
colnames(wt_pseudotime) <- paste0("corrected_", colnames(wt_pseudotime))

merged_seurat <- AddMetaData(merged_seurat, cfko_pseudotime)
merged_seurat <- AddMetaData(merged_seurat, wt_pseudotime)

all_plots <- plotDimRed(merged_seurat, col_by = colnames(cfko_pseudotime),
                        plot_type = "rna.umap")

all_plots <- plotDimRed(merged_seurat, col_by = colnames(wt_pseudotime),
                        plot_type = "rna.umap")

check_correlation <- function(seurat_object, slingshot_df,
                              test_lineage){
  
  old_lineage <- seurat_object[[test_lineage]]
  old_lineage <- old_lineage[rownames(old_lineage) %in%
                               rownames(slingshot_df), , drop = FALSE]
  # Check each lineage
  all_correlations <- lapply(colnames(slingshot_df), function(x){
    return(cor(old_lineage[[1]], slingshot_df[[x]], use = "na.or.complete"))
  })
  names(all_correlations) <- colnames(slingshot_df)
  
  return(all_correlations)
}

# Lineage to test
# WT lineage 1 = lineage 1
# WT lineage 2 = lineage 9
# WT lineage 3 = lineage 7
# CFKO lineage 1 = lineage 1
# CFKO lineage 2 = lineage 4
# CFKO lineage 3 = lineage 2

test_lineage <- c("wt_Lineage1" = "wt", 
                  "wt_Lineage9" = "wt",
                  "wt_Lineage7" = "wt",
                  "cfko_Lineage1" = "cfko",
                  "cfko_Lineage4" = "cfko",
                  "cfko_Lineage2" = "cfko",
                  "cfko_Lineage3" = "cfko")

final_correlations <- lapply(names(test_lineage), function(lineage){
  type <- test_lineage[[lineage]]
  if (type == "wt"){
    slingshot_df <- wt_pseudotime
  } else if (type == "cfko"){
    slingshot_df <- cfko_pseudotime
  } else {
    stop("unrecognized type")
  }
  correlation <- check_correlation(seurat_object = merged_seurat,
                                   slingshot_df = slingshot_df,
                                   test_lineage = lineage)

  return(correlation)
})

names(final_correlations) <- names(test_lineage)

max_values <- lapply(final_correlations, function(correlation){
  max_value <- max(unlist(correlation))
  max_name <- names(correlation[correlation == max_value])
  names(max_value) <- max_name
  return(max_value)
})
# WT lineage 1 = lineage 1
# WT lineage 2 = lineage 9
# WT lineage 3 = lineage 7
# CFKO lineage 1 = lineage 1
# CFKO lineage 2 = lineage 4
# CFKO lineage 3 = lineage 2
lineage_mapping <- c("wt_Lineage1" = "wt_Lineage1",
                     "wt_Lineage9" = "wt_Lineage2",
                     "wt_Lineage7" = "wt_Lineage3",
                     "cfko_Lineage1" = "cfko_Lineage1",
                     "cfko_Lineage4" = "cfko_Lineage2",
                     "cfko_Lineage2" = "cfko_Lineag3")


# Trying genes from the paper with the new values
lineage_colors <- MetBrewer::met.brewer("Egypt", 2)

# Lineage 2 for CF lineage 7 for WT
colnames(cfko_pseudotime) <- paste0("corrected_", colnames(cfko_pseudotime))
colnames(wt_pseudotime) <- paste0("corrected_", colnames(wt_pseudotime))

merged_seurat <- AddMetaData(merged_seurat, cfko_pseudotime)
merged_seurat <- AddMetaData(merged_seurat, wt_pseudotime)

all_plots <- lapply(colnames(cfko_pseudotime), function(x){
  names(lineage_colors) <- c("corrected_wt_Lineage2", x)
  
  different_plots <- plotPseudotime(merged_seurat,
                                    lineages = c(x,
                                                 "corrected_wt_Lineage2"),
                                    gene_list = "CRABP2",
                                    col_by = "lineage", color = lineage_colors,
                                    line_color = "lineage",
                                    alpha = 0.25)[[1]] +
    ggplot2::ggtitle(x)
})

# CFKO lineage 2 looks betst

all_plots <- lapply(colnames(wt_pseudotime), function(x){
  names(lineage_colors) <- c(x, "corrected_cfko_Lineage2")
  
  different_plots <- plotPseudotime(merged_seurat,
                                    lineages = c("corrected_cfko_Lineage2",
                                                 x),
                                    gene_list = "CRABP2",
                                    col_by = "lineage", color = lineage_colors,
                                    line_color = "lineage",
                                    alpha = 0.25)
})

different_plots <- plotPseudotime(merged_seurat,
                                  lineages = c("corrected_cfko_Lineage2",
                                               "corrected_wt_Lineage2"),
                                  gene_list = gene_list,
                                  col_by = "lineage", color = lineage_colors,
                                  line_color = "lineage",
                                  alpha = 0.25)

names(different_plots) <- gene_list


gene_list <- c("ATOX1", "CRABP2", "CFTR", "MUC1")


different_plots <- plotPseudotime(merged_seurat,
                                  lineages = c("cfko_Lineage2",
                                               "wt_Lineage7"),
                                  gene_list = gene_list,
                                  col_by = "lineage", color = lineage_colors,
                                  line_color = "lineage",
                                  alpha = 0.25)

