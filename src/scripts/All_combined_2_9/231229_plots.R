library(scAnalysisR)
library(tidyverse)
library(here)
library(Seurat)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "All_combined_2_9"

normalization_method <- "log" # can be SCT or log

HTO <- FALSE
ADT <- FALSE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
base_dir <- here()
base_dir_proj <- file.path(base_dir, "results", sample)

save_dir <- file.path(base_dir_proj, "R_analysis")

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
# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj/seurat_processed.rds"))

gene_list <- c("ECRG4", "LOC123386996", "SPP1", "LY6D", "SPARC",
               "LOC101674231")


plots_1 <- featDistPlot(seurat_data, geneset = gene_list,
                        sep_by = "sample",
                        col_by = "sample",
                        combine = FALSE,
                        color = sample_colors)


plots_2 <- featDistPlot(seurat_data, geneset = gene_list,
                        sep_by = "celltype",
                        col_by = "celltype",
                        combine = FALSE,
                        color = all_colors)

plots_3 <- featDistPlot(seurat_data, geneset = gene_list,
                        sep_by = "sample",
                        col_by = "celltype",
                        combine = FALSE,
                        color = all_colors)

plots_4 <- featDistPlot(seurat_data, geneset = gene_list,
                        sep_by = "celltype",
                        col_by = "sample",
                        combine = FALSE,
                        color = sample_colors)

fig_dir <- file.path(save_dir, "images", "gene_plots")

ifelse(!dir.exists(fig_dir), dir.create(fig_dir),
       FALSE)

pdf(file.path(fig_dir, "sample_sep.pdf"),
    height = 4, width = 8)

print(plots_1)

dev.off()

pdf(file.path(fig_dir, "cell_type_sep.pdf"),
    height = 4, width = 8)

print(plots_2)

dev.off()

pdf(file.path(fig_dir, "sample_sep_celltype_col.pdf"),
    height = 4, width = 8)

print(plots_3)

dev.off()

pdf(file.path(fig_dir, "celltype_sep_sample_col.pdf"),
    height = 4, width = 8)

print(plots_4)

dev.off()

Idents(seurat_data) <- "celltype"
dot_plot <- DotPlot(seurat_data, features = gene_list) +
  ggplot2::theme(axis.text.x = element_text(angle = 90,
                                            vjust = 0.5, hjust=1))

Idents(seurat_data) <- "sample"

dot_plot2 <- DotPlot(seurat_data, features = gene_list) +
  ggplot2::theme(axis.text.x = element_text(angle = 90,
                                            vjust = 0.5, hjust=1))

pdf(file.path(fig_dir, "dot_plot.pdf"),
    height = 8, width = 8)

print(dot_plot)
print(dot_plot2)

dev.off()