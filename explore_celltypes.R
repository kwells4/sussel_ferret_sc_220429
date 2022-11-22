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

#celltype_column <- "RNA_combined_celltype"
celltype_column <- "new_celltype"

all_sample_dir <- here("results", all_samples, "R_analysis")

all_colors <- c("acinar" = "#D4405B",
                "ductal" = "#A5903E",
                "Prolif_acinar" = "#55A470",
                "Prolif_ductal" = "#767FC9",
                "progenitor_like_cells" = "#297878",
                "transitional_to_acinar" = "#874652",
                "centroacinar" = "#78295D")

new_colors <- c(all_colors,
                "activated_ductal" = "#D67229")

color_by <- new_colors

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

# Gene sets
gene_path <- here("files/GSEA_signaling_pathways_with_orthologs.xlsx")

all_sheets <- openxlsx::getSheetNames(gene_path)

gene_lists <- lapply(all_sheets, function(x){
  gene_df <- openxlsx::readWorkbook(gene_path, sheet = x)
  #all_genes <- unique(gene_df$gene_id)
  return(unique(gene_df$gene_id))
})

all_sheets <- sub(" ", "_", all_sheets)

names(gene_lists) <- all_sheets

 # From:
# Single-cell transcriptome and accessible chromatin dynamics during 
# endocrine pancreas development
sc_list <- list("Ductal" = c("RELA", "CTCF", "VTN",
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

sc_list <- lapply(sc_list, function(x){
  x[x %in% rownames(merged_seurat)]
})

# From:
# Single-cell resolution analysis of the human pancreatic ductal progenitor
# cell niche
quadir_list <- c("KRTAP2-3", "TFF1", "CRP", "WSB1", "CEL", "SYCN",
                 "SRGN", "AKAP12", "IGFBP3", "SPP1", "OLFM4",
                 "CPA2", "IDO1", "ALB", "REG1B", "PRSS1",
                 "CEL3A", "PNLIP", "CEL", "CTRB2", "CPA2",
                 "ALDH1A3", "CFTR", "CRP", "AQP1", "DEFB1",
                 "KRT19", "SPP1", "TSPAN8")

quadir_list <- quadir_list[quadir_list %in% rownames(merged_seurat)]


development_genes <- c("HIF1A", "ETV6", "ID1", "ID2",
                       "HES1", "TEAD2", "BHLHA15", "RFX6",
                       "RUNX1T1", "ETV5", "NR2F2", "MXI1", "FOXA1", "TCF3",
                       "EGR1", "RBPJ", "SNAI2", "NHLH1", "TEAD2",
                       "PAX6", "RFX6", "RUNX1T1", "HMGN3", "CASZ1",
                       "MLXIPL", "EGR1", "ZBTB7C",
                       "FOXP1", "FEV", "NEUROD1")

development_genes <- development_genes[development_genes %in% 
                                         rownames(merged_seurat)]

quadir_pseudotime_genes <- c("TFF1", "SPP1", "CPA2", "AKAP12",
                             "CTNND1", "SOX9", "CD24", "CEACAM6",
                             "PROM1", "F3", "GP2")

quadir_pseudotime_genes <- quadir_pseudotime_genes[quadir_pseudotime_genes %in% 
                                         rownames(merged_seurat)]
# From:
# Exocrine ontogenies: On the development of pancreatic acinar,
# ductal and centroacinar cells
exocrine_ontogenies <- c(
  "ALDH1A1", # centroacinar
  "HES1", # centroacinar
  "SOX9", # centroacinar
  "MUC1", # ductal
  "cytokeritins1", # ductal
  "CA2", # ducatal
  "CFTR", # ductal
  "PROM1" # ductal
)

cytokeritins <- rownames(merged_seurat)[grepl("KRT", rownames(merged_seurat))]

merged_seurat <- AddModuleScore(merged_seurat,
                                features = list("cyto" = cytokeritins),
                                name = "cytokeritins")

UMAP_plots <- plotDimRed(merged_seurat, unlist(unique(sc_list)),
                         plot_type = "rna.umap")


violin_plots1 <- featDistPlot(merged_seurat, unlist(unique(sc_list)),
                              sep_by = celltype_column,
                              combine = FALSE, color = color_by)

merged_seurat$sample_orig_cluster <- paste0(merged_seurat$sample, "_",
                                            merged_seurat$original_cluster)

violin_plots2 <- featDistPlot(merged_seurat, unlist(unique(sc_list)),
                              sep_by = "sample_orig_cluster",
                              combine = FALSE) 

pdf(file.path(all_sample_dir, "images", "cell_type_plots", "sc_plots.pdf"))

print(violin_plots1)

dev.off()




# RELA- most places
# CTCF - most places
# NFIB - highest in the WT progenitor population
# NRF2F6 - higher as you move to the progeintor in both. High in the separate
# starting cluster for the WT
# HES1 - high everywhere. Higheset in the separate cluster for WT
# MYC - highest in the starting cluster for WT
# MECOM - higher in WT than KO
# CARHSP1 - very high in lonely WT cluster
# SMAD3 - high in the progentior populations
# SSRP1 - high in ductal starting cells
# ID3 - separates the two WT ductal populations, higher in the starting ductal 
# for CFKO, low in the progenitor
# HMGB2 - in cycling cells
# HES6 - higher in progenitors
# KLF10 - high everywhere?
# SOX4 - high everywhere?
# STAT3 - higher in progenitors?
# HPB1 - higher in progenitors
# PBX1 - higher in all WT (except the weird population)
# FOXA2 - higher in CFKO
# FOS - higher in progenitors?

UMAP_plots <- plotDimRed(merged_seurat, unique(quadir_list),
                         plot_type = "rna.umap")


violin_plots1 <- featDistPlot(merged_seurat, unique(quadir_list),
                              sep_by = celltype_column,
                              combine = FALSE, color = color_by)

merged_seurat$sample_orig_cluster <- paste0(merged_seurat$sample, "_",
                                            merged_seurat$original_cluster)

violin_plots2 <- featDistPlot(merged_seurat, unique(quadir_list),
                              sep_by = "sample_orig_cluster",
                              combine = FALSE)

pdf(file.path(all_sample_dir, "images", "cell_type_plots", "quadir_plots.pdf"))
print(violin_plots1)

dev.off()


UMAP_plots <- plotDimRed(merged_seurat, unique(development_genes),
                         plot_type = "rna.umap")


violin_plots1 <- featDistPlot(merged_seurat, unique(development_genes),
                              sep_by = celltype_column,
                              combine = FALSE, color = color_by)

merged_seurat$sample_orig_cluster <- paste0(merged_seurat$sample, "_",
                                            merged_seurat$original_cluster)

violin_plots2 <- featDistPlot(merged_seurat, unique(development_genes),
                              sep_by = "sample_orig_cluster",
                              combine = FALSE)



violin_plots1 <- featDistPlot(merged_seurat, unique(exocrine_ontogenies),
                              sep_by = celltype_column,
                              combine = FALSE, color = color_by)


pdf(file.path(all_sample_dir, "images", "cell_type_plots",
              "exocrine_ontogenies_plots.pdf"))
print(violin_plots1)

dev.off()

# Good - expected expression
# ID1, ID2, HES1

# Bad - not expected expression
# HIF1A, NRF2

# Make the plots with pseudotime as well

lineage_colors <- met.brewer("Egypt", 2)
names(lineage_colors) <- c("wt_Lineage4", "cfko_Lineage2")

pseudotime_plot <- plotPseudotime(merged_seurat,
               lineages = c("cfko_Lineage2",
                            "wt_Lineage4"),
               gene_list = unique(development_genes),
               col_by = "lineage", color = lineage_colors,
               line_color = "lineage",
               alpha = 0.25)


pseudotime_plot2 <- plotPseudotime(merged_seurat,
                                  lineages = c("cfko_Lineage2",
                                               "wt_Lineage4"),
                                  gene_list = unique(quadir_pseudotime_genes),
                                  col_by = "lineage", color = lineage_colors,
                                  line_color = "lineage",
                                  alpha = 0.25)
