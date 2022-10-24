library(slingshot)
library(scAnalysisR)
library(Seurat)
library(here)
library(tidyverse)
library(tradeSeq)

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



all_colors <- c("acinar" = "#D4405B",
                "ductal" = "#A5903E",
                "Prolif_acinar" = "#55A470",
                "Prolif_ductal" = "#767FC9",
                "progenitor_like_cells" = "#297878",
                "transitional_to_acinar" = "#874652",
                "centroacinar" = "#78295D")


# cfko slingshot
cfko_slingshot <- readRDS(file.path(all_sample_dir, "files",
                                   "slingshot", "CFKO_res.rds"))


cfko_pseudotime <- slingPseudotime(cfko_slingshot)
colnames(cfko_pseudotime) <- paste0("cfko_", colnames(cfko_pseudotime))

cfko_pseudotime <- data.frame(cfko_pseudotime)

# Check barcodes
dim(cfko_pseudotime)
length(intersect(rownames(cfko_pseudotime), colnames(merged_seurat)))

merged_seurat <- AddMetaData(merged_seurat, metadata = cfko_pseudotime)

pseudotime_plots <- plotDimRed(merged_seurat, colnames(cfko_pseudotime),
                               plot_type= "rna.umap")



# wt slingshot
wt_slingshot <- readRDS(file.path(all_sample_dir, "files",
                                    "slingshot", "WT_res.rds"))


wt_pseudotime <- slingPseudotime(wt_slingshot)
colnames(wt_pseudotime) <- paste0("wt_", colnames(wt_pseudotime))

wt_pseudotime <- data.frame(wt_pseudotime)

# check barcodes
dim(wt_pseudotime)
length(intersect(rownames(wt_pseudotime), colnames(merged_seurat)))

merged_seurat <- AddMetaData(merged_seurat, metadata = wt_pseudotime)

pseudotime_plots <- plotDimRed(merged_seurat, colnames(wt_pseudotime),
                               plot_type= "rna.umap")


# Try out the embedCurves cfko -------------------------------------------------
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
                          plot_type = "rna.umap")[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)  +
    ggplot2::ggtitle(paste0("cfko_lineage", x))
  
  return(base_plot)
})


pdf(file.path(all_sample_dir, "images", "cfko_line_slingshot.pdf"))

print(all_plots)

dev.off()

# I like curves 2, 3, 4, 5, 7, 8

# Plots based on cell type
all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  
  # This line cuts off the long tail... Probably a better option is
  # to just remove unknown before running slingshot.
  base_plot <- plotDimRed(merged_seurat, col_by = "RNA_combined_celltype",
                          color = all_colors,
                          plot_type = "rna.umap")[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)   +
    ggplot2::ggtitle(paste0("cfko_lineage", x))
  
  return(base_plot)
})


pdf(file.path(all_sample_dir, "images", "cfko_line_slingshot_celltype.pdf"))

print(all_plots)

dev.off()


# Plots based on cell type only cells in lineage colored
all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  lineage_name <- paste0("cfko_Lineage", x)
  
  merged_seurat$plot_cells <- ifelse(!is.na(merged_seurat[[lineage_name]][[1]]),
                                     "TRUE", "FALSE")
  
  # This line cuts off the long tail... Probably a better option is
  # to just remove unknown before running slingshot.
  base_plot <- plotDimRed(merged_seurat, col_by = "RNA_combined_celltype",
                          color = all_colors,
                          plot_type = "rna.umap",
                          highlight_group = TRUE,
                          meta_data_col = "plot_cells",
                          group = "TRUE")[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)   +
    ggplot2::ggtitle(paste0("cfko_lineage", x))
  
  return(base_plot)
})


pdf(file.path(all_sample_dir, "images",
              "cfko_line_slingshot_celltype_highlight.pdf"))

print(all_plots)

dev.off()


# Try out the embedCurves wt ---------------------------------------------------
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
                          plot_type = "rna.umap")[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)   +
    ggplot2::ggtitle(paste0("wt_lineage", x))
  
  return(base_plot)
})

# I like curves 2, 3, 4, 5, 7, 8

pdf(file.path(all_sample_dir, "images", "wt_line_slingshot.pdf"))

print(all_plots)

dev.off()



# Plots based on cell type
all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  
  # This line cuts off the long tail... Probably a better option is
  # to just remove unknown before running slingshot.
  base_plot <- plotDimRed(merged_seurat, col_by = "RNA_combined_celltype",
                          color = all_colors,
                          plot_type = "rna.umap")[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)  +
    ggplot2::ggtitle(paste0("wt_lineage", x))
  
  return(base_plot)
})


pdf(file.path(all_sample_dir, "images", "wt_line_slingshot_celltype.pdf"))

print(all_plots)

dev.off()


# Plots based on cell type only cells in lineage colored
all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  lineage_name <- paste0("wt_Lineage", x)
  
  merged_seurat$plot_cells <- ifelse(!is.na(merged_seurat[[lineage_name]][[1]]),
                                     "TRUE", "FALSE")
  
  # This line cuts off the long tail... Probably a better option is
  # to just remove unknown before running slingshot.
  base_plot <- plotDimRed(merged_seurat, col_by = "RNA_combined_celltype",
                          color = all_colors,
                          plot_type = "rna.umap",
                          highlight_group = TRUE,
                          meta_data_col = "plot_cells",
                          group = "TRUE")[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)  +
    ggplot2::ggtitle(paste0("wt_lineage", x))
  
  return(base_plot)
})


pdf(file.path(all_sample_dir, "images",
              "wt_line_slingshot_celltype_highlight.pdf"))

print(all_plots)

dev.off()

# merged_seurat$plot_cells <- NULL

saveRDS(merged_seurat, file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))