library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "All_combined"

wt_samples <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14")

cfko_samples <- c("CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")

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

phase_colors <- LaCroixColoR::lacroix_palette("CranRaspberry", n = 3)
names(phase_colors) <- c("S", "G1", "G2M")

sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))


merged_seurat$sample <- merged_seurat$orig.ident

barplot_res <- scAnalysisR::stacked_barplots(merged_seurat,
                                             meta_col = "RNA_combined_celltype",
                                             split_by = "sample",
                                             color = all_colors,
                                             return_values = TRUE)

pdf(file.path(sample_dir, "images", "cell_type_plots", "cell_type_barplot.pdf"))

print(
barplot_res$barplot +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5, hjust=1))
)
dev.off()

write.csv(barplot_res$data,
          file.path(sample_dir, "files", "celltype_percents.csv"))
