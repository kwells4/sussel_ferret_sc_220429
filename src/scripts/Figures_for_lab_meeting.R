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

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))


save_dir <- here("results", "figures", "20221110_presentation")
panc_db <-
  read.table("~/Documents/Analysis/references/single_cell_references/gene_lists/PanglaoDB_markers_27_Mar_2020.tsv.gz",
             fill = TRUE, header = TRUE, sep = "\t")

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

development_genes <- c("HIF1A", "ETV6", "ID1", "ID2",
                       "HES1", "TEAD2", "BHLHA15", "RFX6",
                       "RUNX1T1", "ETV5", "NR2F2", "MXI1", "FOXA1", "TCF3",
                       "EGR1", "RBPJ", "SNAI2", "NHLH1", "TEAD2",
                       "PAX6", "RFX6", "RUNX1T1", "HMGN3", "CASZ1",
                       "MLXIPL", "EGR1", "ZBTB7C",
                       "FOXP1", "FEV", "NEUROD1")

development_genes <- development_genes[development_genes %in% 
                                         rownames(merged_seurat)]

development_genes <- unique(development_genes)

# Heatmap of genes found to be developmentally regulated
png(file.path(save_dir,
         "developmental_genes.png"), height = 900, width = 800)
print(plot_heatmap(merged_seurat, gene_list = development_genes,
                   meta_col = "sample",
                   average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   cluster_rows = TRUE))

dev.off()


# Violin plot of genes previously found by Pavana ------------------------------
previous_genes <- c("PAX6", "FOXA3", "SOX9", "PDX1")

violin_plots <- featDistPlot(merged_seurat, geneset = previous_genes,
                             sep_by = "sample",
                             color = sample_colors, combine = FALSE)

pdf(file.path(save_dir, "qpcr_violins.pdf"), width = 8, height = 3)

print(violin_plots)

dev.off()


# UMAP across all samples ------------------------------------------------------
individual_plots <- lapply(unique(merged_seurat$sample), function(x){
  plotDimRed(merged_seurat, col_by = "sample", plot_type = "rna.umap",
             highlight_group = TRUE, group = x, meta_data_col = "sample",
             color = sample_colors[x])[[1]]
})

names(individual_plots) <- unique(merged_seurat$sample)

pdf(file.path(save_dir, "umap_sample.pdf"), width = 18,
    height = 18)

print(cowplot::plot_grid(individual_plots$WT_D2, individual_plots$WT_D5,
                   individual_plots$WT_D7, individual_plots$WT_D9,
                   individual_plots$WT_D15, NULL,
                   individual_plots$CFKO_D2,
                   individual_plots$CFKO_D5, individual_plots$CFKO_D7,
                   individual_plots$CFKO_D9, individual_plots$CFKO_D14,
                   NULL, nrow = 4, ncol = 3))

dev.off()


# UMAP of pseudotime on harmony corrected - these look okay but aren't
# perfect

png(file.path(save_dir, "cfko_Lineage2_umap.png"),
    width = 400, height = 400)

print(plotDimRed(merged_seurat, "cfko_Lineage2", plot_type = "harmony.umap",
           highlight_group = TRUE, group = "CFKO", meta_data_col = "genotype"))

dev.off()

plotDimRed(merged_seurat, "cfko_Lineage2", plot_type = "rna.umap",
           highlight_group = TRUE, group = "CFKO", meta_data_col = "genotype")


png(file.path(save_dir, "wt_Lineage4_umap.png"),
    width = 400, height = 400) 
print(plotDimRed(merged_seurat, "wt_Lineage4", plot_type = "harmony.umap",
           highlight_group = TRUE, group = "WT", meta_data_col = "genotype"))

dev.off()
# Pseudotime density plots
# Make density plots for the two lineages colored by cell type
# Get meta data, select only cell type, sample, and lineage,
# keep only values that aren't NA make plot as below:
wt_plotting_df <- merged_seurat[[]] %>%
  dplyr::select(RNA_combined_celltype, sample, wt_Lineage4) %>%
  dplyr::filter(!is.na(wt_Lineage4))

wt_plotting_df$RNA_combined_celltype <- factor(wt_plotting_df$RNA_combined_celltype,
                                               levels = c("Prolif_ductal",
                                                          "ductal",
                                                          "acinar",
                                                          "transitional_to_acinar",
                                                          "centroacinar",
                                                          "progenitor_like_cells"))

pdf(file.path(save_dir, "wt_pseudotime_ridges.pdf"))
print(ggplot2::ggplot(wt_plotting_df, ggplot2::aes(x = wt_Lineage4,
                                                   y = RNA_combined_celltype,
                                                   group = RNA_combined_celltype,
                                                   fill = RNA_combined_celltype)) + 
        ggridges::geom_density_ridges() +
        ggplot2::scale_fill_manual(values = all_colors) +
        ggplot2::ggtitle("WT lineage 4"))

dev.off()

cfko_plotting_df <- merged_seurat[[]] %>%
  dplyr::select(RNA_combined_celltype, sample, cfko_Lineage2) %>%
  dplyr::filter(!is.na(cfko_Lineage2))

cfko_plotting_df$RNA_combined_celltype <- factor(cfko_plotting_df$RNA_combined_celltype,
                                               levels = c("ductal",
                                                          "acinar",
                                                          "transitional_to_acinar",
                                                          "centroacinar",
                                                          "progenitor_like_cells"))


pdf(file.path(save_dir, "cfko_pseudotime_ridges.pdf"))
print(ggplot2::ggplot(cfko_plotting_df, ggplot2::aes(x = cfko_Lineage2,
                                                   y = RNA_combined_celltype,
                                                   group = RNA_combined_celltype,
                                                   fill = RNA_combined_celltype)) + 
        ggridges::geom_density_ridges() +
        ggplot2::scale_fill_manual(values = all_colors) +
        ggplot2::ggtitle("CFKO lineage 2"))

dev.off()

# Pseudotime plots of interesting genes
lineage_colors <- met.brewer("Egypt", 2)
names(lineage_colors) <- c("wt_Lineage4", "cfko_Lineage2")

genes_plot <- c("KCNN4", "LOC101680039",
                "LOC106005009", "MUC4",
                "SFRP1", "TP63",
                "TRIM29", "TXN")

merged_seurat$lineage_col <- ifelse(!is.na(merged_seurat$cfko_Lineage2),
                                    "CFKO",
                                    ifelse(!is.na(merged_seurat$wt_Lineage4),
                                           "WT", NA))


for(gene in genes_plot){
  png(file.path(save_dir, paste0(gene, "_cfko2_wt4.png")),
      width = 400, height = 400)
  
  print(plotPseudotime(merged_seurat,
                 lineages = c("cfko_Lineage2",
                              "wt_Lineage4"),
                 gene_list = gene,
                 col_by = "lineage", color = lineage_colors,
                 line_color = "lineage",
                 alpha = 0.25))
  
  dev.off()
  
  png(file.path(save_dir, paste0(gene, "_umap_cfko.png")),
      width = 400, height = 400)
  
  print(plotDimRed(merged_seurat, gene, plot_type = "harmony.umap",
                   highlight_group = TRUE, group = "CFKO",
                   meta_data_col = "lineage_col"))
  
  dev.off()
        
  
  png(file.path(save_dir, paste0(gene, "_umap_wt.png")),
      width = 400, height = 400)
  
  print(plotDimRed(merged_seurat, gene, plot_type = "harmony.umap",
                   highlight_group = TRUE, group = "WT",
                   meta_data_col = "lineage_col"))
  
  dev.off()
}

