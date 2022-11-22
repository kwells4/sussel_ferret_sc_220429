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

gene_list <- list("Ductal" = c("RELA", "CTCF", "VTN",
                               "NR3C1", "NFIB", "BHLHA15",
                               "NR2F6", "PTF1A", "CSDA", "HES1",
                               "MYC", "HMGA2", "MEIS1", "RBPJL", "MECOM",
                               "CARHSP1", "GATA4", "SMAD3", "MBD2", "SSRP1",
                               "ID3", "TFDP1", "SMARCC1", "NFATC1",
                               "HMGB2", "HMGB3", "SMARAC5", "TCF7L2"),
                  "Early_progenitor" = c("TEAD2", "NCOR2", "HES6", "KLF10",
                                         "NEUROD2", "SOX4", "NEUROG3",
                                         "PAX4", "TOX3", "NKX2-2",
                                         "SMARCE1", "CSDE1"),
                  "Late_progenitor" = c("MYT1", "RXF6", "NEUROD1", "STAT3",
                                        "FEV", "HBP1", "ST18", "PBX1",
                                        "FOXA2", "FOS"))

gene_list <- lapply(gene_list, function(x){
  x[x %in% rownames(merged_seurat)]
})

gene_plot_one <- featDistPlot(merged_seurat, gene_list$Ductal,
                              col_by = "RNA_combined_celltype",
                              color = all_colors,
                              sep_by = "RNA_combined_celltype",
                              combine = FALSE)


gene_plot <- featDistPlot(merged_seurat, "ALDH1A1",
                              col_by = "RNA_combined_celltype",
                              color = all_colors,
                              sep_by = "RNA_combined_celltype",
                              combine = FALSE)


panc_db_celltypes <- c("Acinar cells", "Alpha cells",
                       "Beta cells", "Delta cells", "Ductal cells",
                       "Gamma (PP) cells")

for(gene_list in panc_db_celltypes){
  list_genes <- panc_db %>%
    dplyr::filter(cell.type == gene_list)
  
  plot_genes <- list_genes$official.gene.symbol
  
  plot_genes <- plot_genes[plot_genes %in% rownames(merged_seurat)]
  
  gene_plot_one <- featDistPlot(merged_seurat, plot_genes,
                                col_by = "RNA_combined_celltype",
                                color = all_colors,
                                sep_by = "RNA_combined_celltype",
                                combine = FALSE)
  
  pdf(file.path(all_sample_dir, "images", "cell_type_plots",
                paste0(gene_list, "_plots.pdf")))
  
  print(gene_plot_one)
  
  dev.off()
  
}


quadir_list <- c("KRTAP2-3", "TFF1", "CRP", "WSB1", "CEL", "SYCN",
                 "SRGN", "AKAP12", "IGFBP3", "SPP1", "OLFM4",
                 "CPA2", "IDO1", "ALB", "REG1B", "PRSS1",
                 "CEL3A", "PNLIP", "CEL", "CTRB2", "CPA2",
                 "ALDH1A3", "CFTR", "CRP", "AQP1", "DEFB1",
                 "KRT19", "SPP1", "TSPAN8")

quadir_list <- quadir_list[quadir_list %in% rownames(merged_seurat)]

gene_plot_one <- featDistPlot(merged_seurat, quadir_list,
                              col_by = "RNA_combined_celltype",
                              color = all_colors,
                              sep_by = "RNA_combined_celltype",
                              combine = FALSE)

pdf(file.path(all_sample_dir, "images", "cell_type_plots",
              paste0("quadir_list", "_plots.pdf")))

print(gene_plot_one)

dev.off()



gene_plot_one <- plotDimRed(merged_seurat, quadir_list,
                            plot_type = "rna.umap")

pdf(file.path(all_sample_dir, "images", "cell_type_plots",
              paste0("quadir_list", "_umap.pdf")))

print(gene_plot_one)

dev.off()

Idents(merged_seurat) <- "RNA_combined_celltype"

DotPlot(merged_seurat, features = unique(quadir_list))


# Genes that look good
# https://www.pnas.org/doi/10.1073/pnas.1918314117
# WSB1 - expressed everywhere but higest in the centroacinar
# ALDH1A1 - centroacinar
# IGFBP3 ductal
# KRT19 ductal
# TSPAN8 - progenitor
# Proliferation - HMGB2, STMN1, MKI67

DotPlot(merged_seurat, features = unique(unlist(gene_list))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


final_list <- c("WSB1", "ALDH1A1", "IGFBP3", "KRT19",
                "TSPAN8", "HMGB2", "STMN1", "MKI67")

DotPlot(merged_seurat, features = final_list) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Note, many of the progenitor markeres are good, use those from this list
print(plot_heatmap(merged_seurat, gene_list = unique(unlist(gene_list)),
                   meta_col = "RNA_combined_celltype",
                   average_expression = TRUE,
                   colors = all_colors, plot_rownames = TRUE,
                   cluster_rows = TRUE))


print(plot_heatmap(merged_seurat, gene_list = unlist(quadir_list),
                   meta_col = "RNA_combined_celltype",
                   average_expression = TRUE,
                   colors = all_colors, plot_rownames = TRUE,
                   cluster_rows = TRUE))

print(plot_heatmap(merged_seurat, gene_list = final_list,
                   meta_col = "RNA_combined_celltype",
                   average_expression = TRUE,
                   colors = all_colors, plot_rownames = TRUE,
                   cluster_rows = TRUE))



# Note, many of the progenitor markeres are good, use those from this list
print(plot_heatmap(merged_seurat, gene_list = unique(unlist(gene_list)),
                   meta_col = "sample",
                   average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   cluster_rows = TRUE))


gene_list_violins <- featDistPlot(merged_seurat,
                                  geneset = unique(unlist(gene_list)),
                                  sep_by = "sample",
                                  col_by = "sample",
                                  color = sample_colors, combine = FALSE)

pdf(file.path(all_sample_dir, "images", "cell_type_plots",
              paste0("gene_list", "_sample.pdf")))

print(gene_list_violins)

dev.off()


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


pdf(file.path(all_sample_dir, "images", "cell_type_plots",
              paste0("gene_list", "_sample.pdf")))

print(plot_heatmap(merged_seurat, gene_list = development_genes,
                   meta_col = "sample",
                   average_expression = TRUE,
                   colors = sample_colors, plot_rownames = TRUE,
                   cluster_rows = TRUE))

dev.off()