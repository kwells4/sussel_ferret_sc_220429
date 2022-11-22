library(slingshot)
library(scAnalysisR)
library(Seurat)
library(here)
library(tidyverse)
library(tradeSeq)
library(pheatmap)
library(MetBrewer)

source(here("src/scripts/functions.R"))

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

all_samples <- "All_combined"

all_sample_dir <- here("results", all_samples, "R_analysis")

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

new_colors <- c(all_colors,
                "activated_ductal" = "#D67229")

celltype_levels <- c("acinar", "ductal", "Prolif_acinar", "Prolif_ductal",
                     "progenitor_like_cells", "centroacinar",
                     "transitional_to_acinar")

sample_levels <- names(sample_colors)

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

merged_seurat$celltype_sample <- paste(merged_seurat$RNA_combined_celltype,
                                        merged_seurat$sample, sep = "_")

merged_seurat$RNA_combined_celltype <- factor(merged_seurat$RNA_combined_celltype,
                                              levels = celltype_levels)

merged_seurat$sample <- factor(merged_seurat$sample, levels = sample_levels)

combined_levels <- merged_seurat[[c("sample", "RNA_combined_celltype",
                                   "celltype_sample")]] %>%
  distinct() 

combined_levels <- combined_levels[order(combined_levels$RNA_combined_celltype, 
                                         combined_levels$sample),]
  


merged_seurat$celltype_sample <- factor(merged_seurat$celltype_sample,
                                        levels = combined_levels$celltype_sample)

cytokeritins <- rownames(merged_seurat)[grepl("KRT", rownames(merged_seurat))]

merged_seurat <- AddModuleScore(merged_seurat,
                                features = list("cyto" = cytokeritins),
                                name = "cytokeritins")


best_markers <- c(
  "CTCF", # proliferation
  "HMGB2", # proliferation
  "ID3", # ductal
  "IGFBP3", # activate ductal
  "KRT19", # ductal
  "cytokeritins1", # ductal
  "STAT3", # progenitor
  "FOS", # progenitor
  "TSPAN8", # progenitor
  "SOX4", # early progenitor
  "ALDH1A1" # centroacinar
)


pdf(file.path(all_sample_dir, "images", "cell_type_plots",
              "celltype_heatmaps_sample_celltype.pdf"))



meta_df <- merged_seurat[[c("sample", "RNA_combined_celltype",
                            "celltype_sample")]]


meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$celltype_sample

color_list <- list("sample" = sample_colors,
                   "RNA_combined_celltype" = all_colors)
 
print(plot_heatmap(merged_seurat, gene_list = best_markers,
                   meta_col = "celltype_sample", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   meta_df = meta_ave, color_list = color_list,
                   plot_meta_col = FALSE, cluster_rows = TRUE,
                   cluster_cols = FALSE))

dev.off()

pdf(file.path(all_sample_dir, "images", "cell_type_plots",
              "celltype_heatmaps_celltype.pdf"))

meta_df <- merged_seurat[[c("RNA_combined_celltype")]]


meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$RNA_combined_celltype

color_list <- list("RNA_combined_celltype" = all_colors)

print(plot_heatmap(merged_seurat, gene_list = best_markers,
                   meta_col = "RNA_combined_celltype", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   meta_df = meta_ave, color_list = color_list,
                   plot_meta_col = TRUE, cluster_rows = TRUE,
                   cluster_cols = FALSE))


dev.off()

lineage_colors <- met.brewer("Egypt", 2)
names(lineage_colors) <- c("wt_Lineage4", "cfko_Lineage2")

pseudotime_plot <- plotPseudotime(merged_seurat,
                                  lineages = c("cfko_Lineage2",
                                               "wt_Lineage4"),
                                  gene_list = best_markers,
                                  col_by = "lineage", color = lineage_colors,
                                  line_color = "lineage",
                                  alpha = 0.25)

pdf(file.path(all_sample_dir, "images", "cell_type_plots",
              "celltype_pseudotime_plots.pdf"))

print(pseudotime_plot)

dev.off()




for(ind_sample in names(sample_colors)){
  sample_seurat <- subset(merged_seurat, subset = sample == ind_sample)
  
  violin_plots <- featDistPlot(sample_seurat,
                               geneset = best_markers,
                               color = all_colors,
                               combine = FALSE, sep_by = "sample_cluster",
                               col_by = "RNA_combined_celltype")
  pdf(file.path(all_sample_dir, "images", "cell_type_plots",
                paste0(ind_sample, "_sample_cluster.pdf")))
  
  print(violin_plots)
  
  dev.off()
}


# TODO
# Find clusters that may not match, plot on UMAP
# WT D2 cluster 8 transitional to acinar --> ductal?

# 0, 1, 13, 3 <- activated ductal

merged_seurat$new_celltype <- as.character(merged_seurat$RNA_combined_celltype)

wtd2_activated <- c("WT_D2_0", "WT_D2_1", "WT_D2_13",
                    "WT_D2_3")

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             wtd2_activated] <- "activated_ductal"

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             c("WT_D2_8")] <- "ductal"

# plotDimRed(merged_seurat, col_by = "RNA_combined_celltype", 
#            plot_type = "rna.umap", highlight_group = TRUE,
#            group = "WT_D2_5", meta_data_col = "sample_cluster")
# 
# plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap",
#            color = all_colors)
# 
# plotDimRed(merged_seurat, "new_celltype", plot_type = "rna.umap",
#            color = new_colors)

# WT D5 cluster 10 progenitor like --> ductal

# Cutoff = median 2.5 
wtd5_activated <- c("WT_D5_0", "WT_D5_6", "WT_D5_11")

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             wtd5_activated] <- "activated_ductal"

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             c("WT_D5_10")] <- "ductal"

# plotDimRed(merged_seurat, col_by = "RNA_combined_celltype", 
#            plot_type = "rna.umap", highlight_group = TRUE,
#            group = "WT_D5_12", meta_data_col = "sample_cluster")
# 
# plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap",
#            color = all_colors)
# 
# plotDimRed(merged_seurat, "new_celltype", plot_type = "rna.umap",
#            color = new_colors)

# WT D7 

# Cutoff = median 2.5 
wtd7_activated <- c("WT_D7_0", "WT_D7_4", "WT_D7_7")

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             wtd7_activated] <- "activated_ductal"

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             c("WT_D7_10")] <- "centroacinar"

# plotDimRed(merged_seurat, col_by = "RNA_combined_celltype", 
#            plot_type = "rna.umap", highlight_group = TRUE,
#            group = "WT_D7_10", meta_data_col = "sample_cluster")
# 
# plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap",
#            color = all_colors)
# 
# plotDimRed(merged_seurat, "new_celltype", plot_type = "rna.umap",
#            color = new_colors)


# WT D9

# Cutoff = median 2.5 
wtd9_activated <- c("WT_D9_14", "WT_D9_2")

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             wtd9_activated] <- "activated_ductal"


# plotDimRed(merged_seurat, col_by = "RNA_combined_celltype", 
#            plot_type = "rna.umap", highlight_group = TRUE,
#            group = "WT_D9_1", meta_data_col = "sample_cluster")
# 
# plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap",
#            color = all_colors)
# 
# plotDimRed(merged_seurat, "new_celltype", plot_type = "rna.umap",
#            color = new_colors)

# WT D14

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             c("WT_D14_3",
                               "WT_D14_4")] <- "progenitor_like_cells"

# plotDimRed(merged_seurat, col_by = "RNA_combined_celltype", 
#            plot_type = "rna.umap", highlight_group = TRUE,
#            group = "WT_D14_4", meta_data_col = "sample_cluster")
# 
# plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap",
#            color = all_colors)
# 
# plotDimRed(merged_seurat, "new_celltype", plot_type = "rna.umap",
#            color = new_colors)


# CFKO D2 

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             c("CFKO_D2_11")] <- "progenitor_like_cells"

# plotDimRed(merged_seurat, col_by = "RNA_combined_celltype", 
#            plot_type = "rna.umap", highlight_group = TRUE,
#            group = "CFKO_D2_11", meta_data_col = "sample_cluster")
# 
# plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap",
#            color = all_colors)
# 
# plotDimRed(merged_seurat, "new_celltype", plot_type = "rna.umap",
#            color = new_colors)


# CFKO D5

cfkod5_activated <- c("CFKO_D5_0", "CFKO_D5_2", "CFKO_D5_4",
                      "CFKO_D5_8")

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             cfkod5_activated] <- "activated_ductal"

# plotDimRed(merged_seurat, col_by = "RNA_combined_celltype", 
#            plot_type = "rna.umap", highlight_group = TRUE,
#            group = "CFKO_D5_5", meta_data_col = "sample_cluster")
# 
# plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap",
#            color = all_colors)
# 
# plotDimRed(merged_seurat, "new_celltype", plot_type = "rna.umap",
#            color = new_colors)

# CFKO D7

cfkod7_activated <- c("CFKO_D7_1", "CFKO_D7_2", "CFKO_D7_6")

cfkod7_ductal <- c("CFKO_D7_0", "CFKO_D7_5",
                   "CFKO_D7_8")

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             cfkod7_activated] <- "activated_ductal"


merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             cfkod7_ductal] <- "ductal"

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             c("CFKO_D7_4")] <- "progenitor_like_cells"

# plotDimRed(merged_seurat, col_by = "RNA_combined_celltype", 
#            plot_type = "rna.umap", highlight_group = TRUE,
#            group = "CFKO_D7_8", meta_data_col = "sample_cluster")
# 
# plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap",
#            color = all_colors)
# 
# plotDimRed(merged_seurat, "new_celltype", plot_type = "rna.umap",
#            color = new_colors)


# CFKO D9

cfkod9_activated <- c("CFKO_D9_0", "CFKO_D9_4", "CFKO_D9_7")

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             cfkod9_activated] <- "activated_ductal"


cfkod9_ductal <- c("CFKO_D9_8", "CFKO_D9_9",
                   "CFKO_D9_10")


merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             cfkod9_ductal] <- "ductal"

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             c("CFKO_D9_5")] <- "progenitor_like_cells"

# plotDimRed(merged_seurat, col_by = "RNA_combined_celltype", 
#            plot_type = "rna.umap", highlight_group = TRUE,
#            group = "CFKO_D9_4", meta_data_col = "sample_cluster")
# 
# plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap",
#            color = all_colors)
# 
# plotDimRed(merged_seurat, "new_celltype", plot_type = "rna.umap",
#            color = new_colors)


# CFKO D14

cfkod14_activated <- c("CFKO_D14_5", "CFKO_D14_6")

merged_seurat$new_celltype[merged_seurat$sample_cluster %in%
                             cfkod9_activated] <- "activated_ductal"


# plotDimRed(merged_seurat, col_by = "RNA_combined_celltype", 
#            plot_type = "rna.umap", highlight_group = TRUE,
#            group = "CFKO_D14_5", meta_data_col = "sample_cluster")
# 
# plotDimRed(merged_seurat, "RNA_combined_celltype", plot_type = "rna.umap",
#            color = all_colors)
# 
# plotDimRed(merged_seurat, "new_celltype", plot_type = "rna.umap",
#            color = new_colors)

saveRDS(merged_seurat, file.path(all_sample_dir, "rda_obj",
                                 "seurat_processed.rds"))



for(ind_sample in names(sample_colors)){
  sample_seurat <- subset(merged_seurat, subset = sample == ind_sample)
  
  violin_plots <- featDistPlot(sample_seurat,
                               geneset = best_markers,
                               color = new_colors,
                               combine = FALSE, sep_by = "sample_cluster",
                               col_by = "new_celltype")
  pdf(file.path(all_sample_dir, "images", "cell_type_plots",
                paste0(ind_sample, "_sample_cluster_new_celltype.pdf")))
  
  print(violin_plots)
  
  dev.off()
}

merged_seurat$celltype_sample <- paste(merged_seurat$new_celltype,
                                       merged_seurat$sample, sep = "_")

combined_levels <- merged_seurat[[c("sample", "new_celltype",
                                    "celltype_sample")]] %>%
  distinct() 

combined_levels <- combined_levels[order(combined_levels$new_celltype, 
                                         combined_levels$sample),]

merged_seurat$celltype_sample <- factor(merged_seurat$celltype_sample,
                                        levels = combined_levels$celltype_sample)


pdf(file.path(all_sample_dir, "images", "cell_type_plots",
              "celltype_heatmaps_sample_celltype_updated.pdf"))



meta_df <- merged_seurat[[c("sample", "new_celltype",
                            "celltype_sample")]]


meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$celltype_sample

color_list <- list("sample" = sample_colors,
                   "new_celltype" = new_colors)

print(plot_heatmap(merged_seurat, gene_list = best_markers,
                   meta_col = "celltype_sample", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   meta_df = meta_ave, color_list = color_list,
                   plot_meta_col = FALSE, cluster_rows = TRUE,
                   cluster_cols = FALSE))

dev.off()

pdf(file.path(all_sample_dir, "images", "cell_type_plots",
              "celltype_heatmaps_celltype_updated.pdf"))

meta_df <- merged_seurat[[c("new_celltype")]]


meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$new_celltype

color_list <- list("new_celltype" = new_colors)

print(plot_heatmap(merged_seurat, gene_list = best_markers,
                   meta_col = "new_celltype", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   meta_df = meta_ave, color_list = color_list,
                   plot_meta_col = TRUE, cluster_rows = TRUE,
                   cluster_cols = FALSE))


dev.off()


# UMAP across all samples ------------------------------------------------------
individual_plots <- lapply(unique(merged_seurat$sample), function(x){
  plotDimRed(merged_seurat, col_by = "new_celltype", plot_type = "rna.umap",
             highlight_group = TRUE, group = x, meta_data_col = "sample",
             color = new_colors[x])[[1]]
})

names(individual_plots) <- unique(merged_seurat$sample)

pdf(file.path(all_sample_dir, "images", "cell_type_plots",
              "umap_celltype_updated.pdf"), width = 18,
    height = 18)

print(cowplot::plot_grid(individual_plots$WT_D2, individual_plots$WT_D5,
                         individual_plots$WT_D7, individual_plots$WT_D9,
                         individual_plots$WT_D14, NULL,
                         individual_plots$CFKO_D2,
                         individual_plots$CFKO_D5, individual_plots$CFKO_D7,
                         individual_plots$CFKO_D9, individual_plots$CFKO_D14,
                         NULL, nrow = 4, ncol = 3))

dev.off()

