# Document information
# This document identifies the activated ductal cells, finds markers
# and makes plots for the supplement. It uses the celltypes that
# were originally identified in 03B_rename_celltype_all.R and were generated
# based on the CFKO and WT samples independently.

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

seurat_wt <- readRDS(file.path(save_dir, "rda_obj/seurat_wt.rds"))
seurat_cfko <- readRDS(file.path(save_dir, "rda_obj/seurat_cfko.rds"))

activated_ductal <- openxlsx::readWorkbook(xlsxFile = here("files/pnas.1918314117.sd01.xlsx"),
                                                 sheet = "ALL-CLUSTERS")

# activated_ductal <- activated_ductal %>%
#   dplyr::filter(p_val_adj < 0.05, avg_logFC < 1) %>%
#   dplyr::group_by(gene) %>%
#   dplyr::add_count(name = "gene_count") %>%
#   dplyr::filter(cluster == 2) %>%
#   dplyr::filter(gene_count == 1)

activated_ductal <- activated_ductal %>%
  dplyr::filter(cluster == 2) 

activated_ductal_genes <- activated_ductal$gene

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

mapping_genes <- mapping_file %>%
  dplyr::select(gene_id, dplyr::contains("Human")) %>%
  dplyr::filter(Human.gene.name %in% activated_ductal_genes)

not_found <- activated_ductal_genes[!activated_ductal_genes
                                    %in% mapping_genes$Human.gene.name]

plot_genes <- c(unique(mapping_genes$gene_id), not_found)

# Pull out the clustering from the two subsets (This is the clustering
# that was used to name populations)
wt_meta <- seurat_wt[[]] %>%
  dplyr::select(RNA_cluster) %>%
  dplyr::rename(ind_cluster = RNA_cluster)

cfko_meta <- seurat_cfko[[]] %>%
  dplyr::select(RNA_cluster) %>%
  dplyr::rename(ind_cluster = RNA_cluster)

full_meta <- rbind(cfko_meta, wt_meta)

seurat_data <- AddMetaData(seurat_data, metadata = full_meta)

seurat_data <- AddModuleScore(seurat_data, features = list(plot_genes),
                              name = "activated_ductal")


all_colors2 <- c("acinar" = "#D4405B",
                 "ductal" = "#A5903E",
                 "Prolif_acinar" = "#55A470",
                 "Prolif_ductal" = "#767FC9",
                 "centroacinar_progenitor" = "#297878",
                 "centroacinar" = "#78295D")

new_colors <- c(all_colors2,
                "activated_ductal" = "#D67229")

celltype_levels <- c("acinar", "ductal", "Prolif_acinar", "Prolif_ductal",
                     "progenitor_like_cells", "centroacinar",
                     "transitional_to_acinar")

sample_colors <- as.character(LaCroixColoR::lacroix_palette("Coconut", 10))
sample_colors[5] <- "#F4E3C7"
names(sample_colors) <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                          "CFKO_D14", "CFKO_D9", "CFKO_D7", "CFKO_D5", "CFKO_D2")
sample_colors <- sample_colors[c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                                 "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")]

genotype_colors <- sample_colors[grepl("D2", names(sample_colors))]
names(genotype_colors) <- gsub("_D2", "", names(genotype_colors))

sample_colors <- sample_colors[!(grepl("D14", names(sample_colors)))]

fig_dir <- file.path(save_dir, "images", "final_figures")

ifelse(!dir.exists(fig_dir), dir.create(fig_dir), FALSE)

seurat_data$celltype <- seurat_data$final_ind_celltype

seurat_data$sample <- factor(seurat_data$sample,
                             levels = names(sample_colors))


# IGFBP3 is the best marker of activated ductal
all_data <- GetAssayData(seurat_data, slot = "data", assay = "RNA")
IGFBP3_data <- all_data["IGFBP3",]

all_meta <- seurat_data[[]] %>%
  dplyr::select(sample, final_ind_celltype, ind_cluster,
                genotype, activated_ductal1) %>%
  dplyr::mutate(genotype_cluster_celltype = 
                  paste(genotype, ind_cluster, final_ind_celltype,
                        sep = "_"))


if(!identical(rownames(all_meta),
              names(IGFBP3_data))){
  all_meta <- all_meta[order(match(rownames(all_meta),
                                   names(IGFBP3_data)))]
}

all_meta$IGFBP3_expression <- IGFBP3_data

# Figure out the cutoff score that's best
all_meta <- all_meta %>%
  dplyr::group_by(genotype_cluster_celltype) %>%
  dplyr::mutate(median_expression = median(IGFBP3_expression),
                median_score = median(activated_ductal1)) %>%
  dplyr::select(-IGFBP3_expression, -sample, -activated_ductal1) %>%
  dplyr::distinct() 

featDistPlot(seurat_data, geneset = "activated_ductal1",
             col_by = "final_ind_celltype", sep_by = "genotype",
             combine = FALSE, color = all_colors2)

featDistPlot(seurat_data, geneset = "activated_ductal1",
             col_by = "genotype", sep_by = "final_ind_celltype",
             combine = FALSE, color = genotype_colors)


seurat_wt <- subset(seurat_data, subset = genotype == "WT")

keep_cells <- seurat_wt[[]] %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::group_by(sample, final_ind_celltype) %>%
  dplyr::add_count(name = "total_cells") %>%
  dplyr::filter(total_cells > 5)

p1 <- featDistPlot(seurat_wt, geneset = "activated_ductal1",
             col_by = "final_ind_celltype", sep_by = "final_ind_celltype",
             combine = FALSE, color = all_colors2)[[1]] +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::ggtitle("WT acitviated ductal by celltype")

seurat_wt <- subset(seurat_wt, cells = keep_cells$barcode)

p2 <- featDistPlot(seurat_wt, geneset = "activated_ductal1",
             col_by = "sample", sep_by = "final_ind_celltype",
             combine = FALSE, color = sample_colors)[[1]] +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::ggtitle("WT acitviated ductal by sample and celltype")


seurat_cfko <- subset(seurat_data, subset = genotype == "CFKO")

keep_cells <- seurat_cfko[[]] %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::group_by(sample, final_ind_celltype) %>%
  dplyr::add_count(name = "total_cells") %>%
  dplyr::filter(total_cells > 5)

p3 <- featDistPlot(seurat_cfko, geneset = "activated_ductal1",
             col_by = "final_ind_celltype", sep_by = "final_ind_celltype",
             combine = FALSE, color = all_colors2)[[1]] +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::ggtitle("CFKO acitviated ductal by celltype")

seurat_cfko <- subset(seurat_cfko, cells = keep_cells$barcode)

p4 <- featDistPlot(seurat_cfko, geneset = "activated_ductal1",
             col_by = "sample", sep_by = "final_ind_celltype",
             combine = FALSE, color = sample_colors)[[1]] +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::ggtitle("CFKO acitviated ductal by sample and celltype")

pdf(file.path(save_dir, "images", "activated_ductal",
              "activated_ductal_module_violin.pdf"), width = 5, height = 4)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

seurat_data$sample_celltype <- paste(seurat_data$sample,
                                     seurat_data$final_ind_celltype, sep = "_")


meta_df <- seurat_data[[c("sample", "final_ind_celltype",
                          "sample_celltype")]]


meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$sample_celltype

color_list <- list("sample" = sample_colors,
                   "final_ind_celltype" = all_colors2)

meta_ave <- meta_ave %>%
  dplyr::arrange( final_ind_celltype)

graphics.off()
pdf(file.path(save_dir, "images", "activated_ductal",
              "activated_ductal_gene_heatmap_sample_celltype.pdf"), width = 10, height = 8)

plot_heatmap(seurat_data, gene_list = plot_genes,
             meta_col = "sample_celltype", color_list = color_list,
             meta_df = meta_ave, plot_meta_col = FALSE,
             average_expression = TRUE, cluster_rows = TRUE)

dev.off()

graphics.off()
pdf(file.path(save_dir, "images", "activated_ductal",
              "activated_ductal_gene_heatmap_sample.pdf"), width = 5, height = 8)


plot_heatmap(seurat_data, gene_list = plot_genes,
             meta_col = "sample", colors = sample_colors,
             average_expression = TRUE, cluster_rows = TRUE)
dev.off()

pdf(file.path(save_dir, "images", "activated_ductal",
              "activated_ductal_module.pdf"), width = 5, height = 4)

# CFKO has higher scores in general than WT
print(ggplot2::ggplot(all_meta, ggplot2::aes(x = median_expression,
                                       y = median_score,
                                       color = genotype)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = genotype_colors))

print(plotDimRed(seurat_data, "activated_ductal1", plot_type = "rna.umap"))
print(plotDimRed(seurat_data, "IGFBP3", plot_type = "rna.umap"))

dev.off()


all_meta <- all_meta %>%
  dplyr::filter(median_expression > 2.5) 
 # dplyr::filter(median_score > 0.075)


seurat_data$genotype_cluster_celltype <-paste(
  seurat_data$genotype, seurat_data$ind_cluster,
  seurat_data$final_ind_celltype, sep = "_"
)

seurat_data$new_celltype <- seurat_data$final_ind_celltype
seurat_data$new_celltype <- as.character(seurat_data$new_celltype)
seurat_data$new_celltype[seurat_data$genotype_cluster_celltype %in%
                           all_meta$genotype_cluster_celltype] <- "activated_ductal"


featDistPlot(seurat_data, geneset = c("activated_ductal1", "IGFBP3"), 
             sep_by = "new_celltype", combine = FALSE, 
             color = new_colors)

scAnalysisR::plot_heatmap(seurat_data, gene_list = plot_genes,
                          meta_col = "new_celltype", 
                          average_expression = TRUE, 
                          colors = new_colors,
                          cluster_rows = TRUE)


scAnalysisR::plot_heatmap(seurat_data, gene_list = plot_genes,
                          meta_col = "final_ind_celltype", 
                          average_expression = TRUE, 
                          colors = all_colors2,
                          cluster_rows = TRUE)

plotDimRed(seurat_data, "new_celltype", plot_type = "rna.umap",
           color = new_colors)

# Find conserved markers -------------------------------------------------------
Idents(seurat_data) <- "new_celltype"
all_markers <- FindConservedMarkers(object = seurat_data,
                                    ident.1 = "activated_ductal",
                                    grouping.var = "sample",
                                    min.cells.group = 15)


sig_markers <- all_markers %>%
  dplyr::filter(max_pval < 0.05)

sig_markers$gene <- rownames(sig_markers)

# conserved markers ------------------------------------------------------------

mean_logfc <- sig_markers %>%
  select(contains("avg_log2FC")) %>%
  dplyr::mutate("logfc_mean" = rowMeans(.), gene = rownames(.)) %>%
  dplyr::filter(logfc_mean > 0)

# START HERE --------------------------------------------------------------------

plot_seurat <- subset(seurat_data,
                      subset = new_celltype %in% c("ductal", "activated_ductal"))

plot_seurat$celltype_sample <- paste0(plot_seurat$new_celltype,
                                      "_", plot_seurat$sample)

meta_df <- plot_seurat[[c("sample", "new_celltype",
                          "celltype_sample")]]


meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$celltype_sample

color_list <- list("sample" = sample_colors,
                   "new_celltype" = new_colors)

meta_ave <- meta_ave %>%
  dplyr::arrange(new_celltype)

ifelse(!dir.exists(file.path(save_dir, "images", "activated_ductal")),
       dir.create(file.path(save_dir, "images", "activated_ductal")),
       FALSE)

graphics.off()

pdf(file.path(save_dir, "images", "activated_ductal",
              "activated_ductal_de_sample_group.pdf"),
    height = 12)

heatmap_data <- plot_heatmap(plot_seurat, gene_list = mean_logfc$gene,
                             meta_col = "celltype_sample", average_expression = TRUE,
                             colors = sample_colors, plot_rownames = TRUE,
                             meta_df = meta_ave, color_list = color_list,
                             plot_meta_col = FALSE, cluster_rows = TRUE,
                             cluster_cols = FALSE, return_data = TRUE)

print(heatmap_data$heatmap)


dev.off()

write.csv(heatmap_data$z_score, file = file.path(save_dir, "images",
                                                 "activated_ductal",
                                                 "heatmap_values.csv"))

plot_data <- plotDimRed(seurat_data, col_by = c("IGFBP3", "CACNA2D1"),
                        plot_type = "rna.umap", ggrastr = TRUE)

names(plot_data) <- c("IGFBP3", "CACNA2D1")

pdf(file.path(save_dir, "images", "activated_ductal",
              "igfpb3.pdf"), width = 5, height = 4)

print(plot_data[[1]])

dev.off()


pdf(file.path(save_dir, "images", "activated_ductal",
              "cacna2d1.pdf"), width = 5, height = 4)

print(plot_data[[2]])

dev.off()

row_order <- rownames(plot_data[[1]]$data)

return_df <- lapply(names(plot_data), function(gene){
  plot <- plot_data[[gene]]$data %>%
    dplyr::select(colour_metric)
  
  if(!identical(rownames(plot), row_order)){
    plot <- plot[order(match(rownames(plot), row_order)),]
  }
  colnames(plot) <- gene
  return(plot)
})

return_df <- do.call(cbind, return_df)
umap_data <- plot_data[[1]]$data %>%
  dplyr::select(dim1, dim2) %>%
  dplyr::rename(UMAP1 = dim1, UMAP2 = dim2)

return_data <- cbind(umap_data, return_df)

write.csv(heatmap_data$z_score, file = file.path(save_dir, "images",
                                                 "activated_ductal",
                                                 "umap_values.csv"))

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))

