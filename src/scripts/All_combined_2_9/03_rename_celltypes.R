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

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj/seurat_processed.rds"))

# Colors -----------------------------------------------------------------------
all_colors <- c("acinar" = "#D4405B",
                "ductal" = "#A5903E",
                "Prolif_acinar" = "#55A470",
                "Prolif_ductal" = "#767FC9",
                "progenitor_like_cells" = "#297878",
                "transitional_to_acinar" = "#874652",
                "centroacinar" = "#78295D")


all_colors2 <- c("acinar" = "#D4405B",
                 "ductal" = "#A5903E",
                 "Prolif_acinar" = "#55A470",
                 "Prolif_ductal" = "#767FC9",
                 "centroacinar2" = "#297878",
                 "centroacinar1" = "#78295D")

sample_colors <- as.character(LaCroixColoR::lacroix_palette("Coconut", 10))
sample_colors[5] <- "#F4E3C7"
names(sample_colors) <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                          "CFKO_D14", "CFKO_D9", "CFKO_D7", "CFKO_D5", "CFKO_D2")
sample_colors <- sample_colors[c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                                 "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")]

sample_colors <- sample_colors[!(grepl("D14", names(sample_colors)))]

# Figure out cell types --------------------------------------------------------

all_plots <- lapply(names(sample_colors), function(x){
  plotDimRed(seurat_data, col_by = "orig_combined_celltype",
             plot_type = "rna.umap", highlight_group = TRUE,
             group = x, meta_data_col = "sample",
             ggrastr = TRUE)[[1]]
})

names(all_plots) <- names(sample_colors)

seurat_sub <- subset(seurat_data, 
                     subset = RNA_combined_celltype %in%
                       c("progenitor_like_cells",
                         "transitional_to_acinar",
                         "centroacinar"))

featDistPlot(seurat_sub, "ALDH1A1", sep_by = "sample",
             col_by = "orig_combined_celltype", combine = FALSE)

featDistPlot(seurat_sub, "ALDH1A1", col_by = "sample",
             sep_by = "orig_combined_celltype", combine = FALSE,
             color = sample_colors)

featDistPlot(seurat_sub, "ALDH1A1", col_by = "orig_combined_celltype",
             sep_by = "sample", combine = FALSE)


all_plots <- featDistPlot(seurat_data, c("ALDH1A1", "KRT19", "TSPAN8", "ID1"),
                          col_by = "orig_combined_celltype", 
                          combine = FALSE,
                          sep_by = "genotype")

# Rename cell types ------------------------------------------------------------

seurat_name_mapping <- unique(seurat_data$orig_combined_celltype)
names(seurat_name_mapping) <- seurat_name_mapping

seurat_name_mapping[seurat_name_mapping == "progenitor_like_cells"] = 
  "centroacinar-progenitor"

seurat_name_mapping[seurat_name_mapping == "transitional_to_acinar1"] = 
  "centroacinar-progenitor"

seurat_name_mapping[seurat_name_mapping == "centroacinar"] = 
  "centroacinar"

seurat_name_mapping[seurat_name_mapping == "transitional_to_acinar2"] = 
  "centroacinar-progenitor"

seurat_data$final_celltype <- seurat_name_mapping[as.character(seurat_data$orig_combined_celltype)]

seurat_data$final_celltype <- factor(seurat_data$final_celltype,
                                       levels = names(all_colors2))

all_barplots <- scAnalysisR::stacked_barplots(seurat_data,
                                              meta_col = "final_celltype",
                                              split_by = "sample",
                                              return_values = TRUE,
                                              color = all_colors2)

print(
  all_barplots$barplot +
    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 0.5, hjust=1))
)


seurat_sub <- subset(seurat_data, 
                     subset = RNA_combined_celltype %in%
                       c("progenitor_like_cells",
                         "transitional_to_acinar",
                         "centroacinar"))

featDistPlot(seurat_sub, "ALDH1A1", col_by = "sample",
             sep_by = "final_celltype", combine = FALSE,
             color = sample_colors)

saveRDS(seurat_data, file.path(save_dir, "rda_obj/seurat_processed.rds"))