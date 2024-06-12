# Document information
# This document makes all of the figures that are seen in the manuscript.

library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(viridis)
library(clustree)
library(harmony)
library(openxlsx)
source("src/scripts/functions.R")

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "All_combined_2_9"

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))
mapping_file <- read.csv(here("files/species_mapping_file.csv"))

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

# DE ---------------------------------------------------------------------------
### Ductal DE ------------------------------------------------------------------
# First find the three populations
seurat_wt <- readRDS(file.path(save_dir, "rda_obj/seurat_wt.rds"))

wt_clusters <- seurat_wt[["RNA_cluster"]]

colnames(wt_clusters) <- c("WT_clusters")

seurat_data <- AddMetaData(seurat_data, wt_clusters)

graphics.off()
plotDimRed(seurat_data, col_by = "WT_clusters", plot_type = "rna.umap", 
           highlight_group = TRUE, group = "ductal",
           meta_data_col = "celltype")

ductal_mapping <- c("0" = "ductal_3",
                    "11" = "ductal_3",
                    "14" = "ductal_1",
                    "15" = "ductal_2",
                    "17" = "ductal_3",
                    "19" = "ductal_3",
                    "2" = "ductal_1",
                    "20" = "ductal_1",
                    "22" = "ductal_3",
                    "23" = "ductal_1",
                    "31" = "ductal_3",
                    "33" = "ductal_3",
                    "34" = "ductal_3",
                    "35" = "ductal_3",
                    "37" = "ductal_2",
                    "39" = "ductal_2",
                    "4" = "ductal_2",
                    "5" = "ductal_1",
                    "6" = "ductal_2",
                    "8" = "ductal_3",
                    "9" = "ductal_2")

seurat_data$ductal_celltype <- seurat_data$celltype
seurat_data$ductal_celltype <- ductal_mapping[seurat_data$WT_clusters]

ductal_colors <- MetBrewer::met.brewer(palette_name = "Greek",
                                       n = 5)[3:5]

names(ductal_colors) <- c("ductal_1", "ductal_2", "ductal_3")

graphics.off()

pdf(file.path(save_dir, "images", "revision", "wt_ductal_umap.pdf"))

print(plotDimRed(seurat_data, col_by = "ductal_celltype", plot_type = "rna.umap",
                 color = ductal_colors, highlight_group = TRUE,
                 group = c("ductal_1", "ductal_2", "ductal_3"),
                 meta_data_col = "ductal_celltype", ggrastr = TRUE))

dev.off()


ductal_seurat <- subset(seurat_data, subset = ductal_celltype %in%
                          c("ductal_1", "ductal_2", "ductal_3"))


Idents(ductal_seurat) <- "ductal_celltype"

ductal_markers <- find_write_markers_orig_orth(seurat_object = ductal_seurat,
                                               save_dir = file.path(save_dir),
                                               meta_col = "ductal_celltype", 
                                               mapping_file = mapping_file,
                                               mapping_gene_col = "gene_id",
                                               mapping_ortholog_col = c("Human.gene.name",
                                                                        "Mouse.gene.name",
                                                                        "Pig.gene.name"),
                                               assay = "RNA", pval = 0.05,
                                               logfc = 0.5, gene_lists = NULL) 

# Save results to excel
save_wb <- openxlsx::createWorkbook()

for (ductal_pop in c("ductal_1", "ductal_2", "ductal_3")){
  keep_markers <- ductal_markers %>%
    dplyr::filter(cluster == ductal_pop & p_val_adj < 0.05)
  
  openxlsx::addWorksheet(wb = save_wb, sheetName = ductal_pop)
  openxlsx::writeData(wb = save_wb, sheet = ductal_pop,
                      x = keep_markers)
}

openxlsx::saveWorkbook(save_wb, file.path(save_dir, "images", "revision",
                                          "WT_ductal_DE.xlsx"),
                       overwrite = TRUE)

# Make a heatmap of the top 10 for each group
plot_genes <- ductal_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 10)

gene_mapping <- c("LOC101680039" = "KRT13*",
                  "LOC101672378" = "CXCL1**",
                  "LOC106003392" = "SAA3*", 
                  "LOC101672959" = "KRT15**",
                  "LOC101672555" = "CDKN2A**",
                  "LOC123390887" = "PRAP1**",
                  "LOC101679740" = "KRT14**",
                  "LOC101686051" = "CCL2*")

graphics.off()

plot_genes$final_gene <- gene_mapping[plot_genes$gene]
plot_genes$final_gene <- ifelse(is.na(plot_genes$final_gene), plot_genes$gene,
                                plot_genes$final_gene)

plot_genes$cluster <- factor(plot_genes$cluster, 
                             levels = names(ductal_colors))

plot_genes <- plot_genes %>%
  dplyr::arrange(cluster)

pdf(file.path(save_dir, "images", "revision",
              "ductal_de_heatmap.pdf"),
    height = 10, width = 8)

heatmap <- plot_heatmap(ductal_seurat, gene_list = plot_genes$gene,
                        meta_col = "ductal_celltype", average_expression = TRUE,
                        colors = ductal_colors, plot_rownames = TRUE,
                        return_data = TRUE, labels_row = plot_genes$final_gene)

print(heatmap$heatmap)

dev.off()

z_score <- data.frame(heatmap$z_score)
z_score$updated_gene <- gene_mapping[rownames(z_score)]
z_score$updated_gene <- ifelse(is.na(z_score$updated_gene), rownames(z_score),
                               z_score$updated_gene)


write.csv(heatmap$z_score, file.path(save_dir, "images", "revision",
                                     "ductal_de_z_score.csv" ))

graphics.off()

### Cell type DE ---------------------------------------------------------------
seurat_data$sample_celltype <- paste(seurat_data$sample,
                                     seurat_data$celltype,
                                     sep = "_")

# Grab out results for that test
de_test <- "D9_combined_celltype"
de_directory <- file.path(save_dir, "files", "DE")

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
save_wb <- openxlsx::createWorkbook()
keep_celltype <- "ductal"

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
pdf(file.path(save_dir, "images", "revision",
              paste0(de_test, "_", keep_celltype,
                     "_heatmap_average.pdf")))
heatmap <- plot_heatmap(celltype_seurat, gene_list = unique(test_de$gene),
                        meta_col = "sample", average_expression = TRUE,
                        colors = sample_colors, plot_rownames = FALSE,
                        return_data = TRUE)

print(heatmap$heatmap)

dev.off()

openxlsx::addWorksheet(wb = save_wb, sheetName = keep_celltype)
openxlsx::writeData(wb = save_wb, sheet = keep_celltype,
                    x = heatmap$z_score,
                    rowNames = TRUE)

fig_height <- round(length(test_de$gene) / 10)

pdf(file.path(save_dir, "images", "revision",
              paste0(de_test, "_", keep_celltype,
                     "_heatmap_average_names.pdf")),
    width = 8, height = fig_height)
print(plot_heatmap(celltype_seurat, gene_list = unique(test_de$gene),
                   meta_col = "sample", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE))

dev.off()

new_gene_list <- test_de %>%
  dplyr::select(gene, Human.gene.name) %>%
  dplyr::distinct() %>%
  dplyr::group_by(gene) %>%
  dplyr::add_count(name = "ferret_count")

one_genes <- new_gene_list %>%
  dplyr::filter(ferret_count == 1)

multi_genes <- new_gene_list %>%
  dplyr::filter(ferret_count > 1) %>%
  dplyr::filter(Human.gene.name != "") %>%
  dplyr::group_by(gene) %>%
  dplyr::add_count(name = "ferret_count")

one_genes_pt_two <- multi_genes %>%
  dplyr::filter(ferret_count == 1)

multi_genes <- multi_genes %>%
  dplyr::filter(ferret_count > 1)

# I'm renaming these by hand
# updated_genes <- c("LOC101682302" = "CEACAM1",
#                    "LOC101683606" = "DDT",
#                    "LOC101676319" = "IFITM3",
#                    "LOC101690213" = "APOL3",
#                    "LOC123394022" = "EIF2S3",
#                    "UPK3BL2" = "UPK3BL2",
#                    "LOC106005009" = "MT1A",
#                    "LOC101680117" = "MT1E",
#                    "EIF5A" = "EIF5A")

updated_genes <- c("")

all_genes <- c(paste(one_genes$gene, one_genes$Human.gene.name, sep = "_"),
               paste(one_genes_pt_two$gene, one_genes_pt_two$Human.gene.name,
                     sep = "_"),
               paste(names(updated_genes), updated_genes, sep = "_"))
names(all_genes) <- c(one_genes$gene, one_genes_pt_two$gene, 
                      names(updated_genes))



pdf(file.path(save_dir, "images", "revision",
              paste0(de_test, "_", keep_celltype,
                     "_heatmap_average_human_names.pdf")),
    width = 8, height = fig_height)
print(plot_heatmap(celltype_seurat, gene_list = unique(test_de$gene),
                   meta_col = "sample", average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   labels_row = all_genes))

dev.off()

openxlsx::saveWorkbook(save_wb, file.path(save_dir, "images", "revision",
                                          "D9_DE_ductal_heatmaps.xlsx"),
                       overwrite = TRUE)


# Batch correction -------------------------------------------------------------
seurat_data <- RunHarmony(seurat_data, c("day"),
                          plot_convergence = TRUE,
                          theta = 3)

harmony_plots <- plotDimRed(seurat_data, col_by = "sample",
                            plot_type = "harmony")

# UMAP -------------------------------------------------------------------------

RNA_pcs <- 30
ADT_pcs <- 8

# Remove previous clustering
remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                               colnames(seurat_data[[]]))]

for (i in remove_cols){
  seurat_data[[i]] <- NULL
}


set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = "RNA", HTO = FALSE,
                         reduction = "harmony")

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

seurat_data <- FindClusters(seurat_data, resolution = c(0.5, 0.8, 1, 1.2,
                                                        1.5, 2.0, 3.0))
clustree(seurat_data)

# UMAP of gene expression
set.seed(0)
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 3.0, assay = "RNA", HTO = FALSE,
                         reduction = "harmony")


seurat_data <- umap_data$object

gene_plots <- umap_data$plots

seurat_data$full_corrected_cluster <- seurat_data$RNA_cluster

cm <- confusionMatrix(seurat_data$full_corrected_cluster, seurat_data$celltype )
cm <- cm / rowSums(cm)

graphics.off()
pdf(file.path(save_dir, "images", "revision",
              "full_batch_corrected_cluster_vs_celltype.pdf"))
print(pheatmap::pheatmap(cm))

dev.off()

pdf(file.path(save_dir, "images", "revision", 
              "celltype_on_corrected_umap.pdf"))
print(plotDimRed(seurat_data, col_by = "celltype", plot_type = "harmony.umap",
                 color = all_colors2))

print(plotDimRed(seurat_data, col_by = "sample", plot_type = "harmony.umap"))

dev.off()
