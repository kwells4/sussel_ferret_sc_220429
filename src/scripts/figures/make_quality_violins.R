library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

samples <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
             "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")

merged_seurat <- "All_combined"

# Sample colors
sample_colors <- as.character(LaCroixColoR::lacroix_palette("Coconut", 10))
sample_colors[5] <- "#F4E3C7"
names(sample_colors) <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                          "CFKO_D14", "CFKO_D9", "CFKO_D7", "CFKO_D5", "CFKO_D2")
sample_colors <- sample_colors[c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                                 "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")]


# Set directories
base_dir <- here()

seurat_objs <- lapply(samples, function(x){
  print(x)

  base_dir_proj <- file.path(base_dir, "results", x)
  save_dir <- file.path(base_dir_proj, "R_analysis")
  
  # Read in the data
  seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_unfilt.rds"))
  
  return(seurat_data)
  
})


all_seurat <- merge(seurat_objs[[1]], seurat_objs[2:length(seurat_objs)])

all_seurat$orig.ident <- factor(all_seurat$orig.ident,
                                levels = samples)

pdf(here("results/figures/quality_unfiltered.pdf"), width = 6,
         height = 4)

print(featDistPlot(all_seurat, geneset = "nFeature_RNA",
                   combine = FALSE, sep_by = "orig.ident",
                   color = sample_colors))

print(featDistPlot(all_seurat, geneset = "nCount_RNA",
                   combine = FALSE, sep_by = "orig.ident",
                   color = sample_colors))

print(featDistPlot(all_seurat, geneset = "percent.mt",
                   combine = FALSE, sep_by = "orig.ident",
                   color = sample_colors))

dev.off()

merged_seurat <- readRDS(file.path(here("results/All_combined",
                                        "R_analysis", "rda_obj",
                                        "seurat_processed.rds")))

merged_seurat$orig.ident <- factor(merged_seurat$orig.ident,
                                   levels = samples)

pdf(here("results/figures/quality_filtered.pdf"), width = 6,
         height = 4)

print(featDistPlot(merged_seurat, geneset = "nFeature_RNA",
                   combine = FALSE, sep_by = "orig.ident",
                   color = sample_colors))

print(featDistPlot(merged_seurat, geneset = "nCount_RNA",
                   combine = FALSE, sep_by = "orig.ident",
                   color = sample_colors))

print(featDistPlot(merged_seurat, geneset = "percent.mt",
                   combine = FALSE, sep_by = "orig.ident",
                   color = sample_colors))

dev.off()
