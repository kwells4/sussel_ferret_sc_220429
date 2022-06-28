library(Seurat)
library(here)

sample <- "All_combined"


# Set directories
base_dir <- here()
base_dir_proj <- file.path(base_dir, "results", sample)

save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj/seurat_processed.rds"))

all_genes <- saveRDS(rownames(seurat_data), "src/scripts/shiny/all_genes.rds")