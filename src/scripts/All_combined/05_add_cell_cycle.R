library(Seurat)
library(here)
library(tidyverse)
library(LaCroixColoR)
library(openxlsx)

sample <- "All_combined"

all_samples <- c("CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14", 
                 "WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14")


phase_colors <- LaCroixColoR::lacroix_palette("CranRaspberry", n = 3)
names(phase_colors) <- c("S", "G1", "G2M")

sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

merged_meta <- merged_seurat[[]] %>%
  dplyr::select(sample) %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::mutate(count = gsub(".*_", "", barcode))


individual_meta <- lapply(seq_along(all_samples), function(x){
  print(x)
  print(all_samples[[x]])
  sample_name <- all_samples[[x]]
  seurat_data <- readRDS(file.path("results", sample_name, "R_analysis",
                                   "rda_obj", "seurat_processed.rds"))
  
  seurat_meta <- seurat_data[[]] %>%
    dplyr::select(orig.ident, Phase) %>%
    tibble::rownames_to_column("barcode_old") %>%
    dplyr::mutate(barcode = paste(barcode_old, x, sep = "_"),
                  count = x)
  
  return(seurat_meta)
})


all_meta <- do.call(rbind, individual_meta)

all_meta$barcode_old <- NULL

all.equal(merged_meta$barcode, all_meta$barcode)


all_meta <- all_meta %>%
  tibble::column_to_rownames("barcode")

merged_seurat <- AddMetaData(merged_seurat, all_meta)

plotDimRed(merged_seurat, col_by = "Phase", plot_type = "rna.umap",
           color = c(phase_colors))

# Make excel file --------------------------------------------------------------

write_cell_cycle_to_excel <- function(sample_name, seurat_object,
                                      excel_workbook){
  
  if(sample_name == "all_cells"){
    merged_meta <- seurat_object[[]]
  } else {
    merged_meta <- seurat_object[[]] %>%
      dplyr::filter(orig.ident == sample_name)
  }
  
  merged_meta <- merged_meta %>%
    dplyr::select(Phase, orig.ident, RNA_combined_celltype) %>%
    dplyr::group_by(RNA_combined_celltype) %>%
    dplyr::add_count(name = "total_celltype") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(RNA_combined_celltype, Phase) %>%
    dplyr::add_count(name = "Phase_by_celltype") %>%
    dplyr::mutate(percent_phase = Phase_by_celltype / total_celltype) %>%
    dplyr::distinct(RNA_combined_celltype, Phase, .keep_all = TRUE) %>%
    dplyr::select(Phase, RNA_combined_celltype, percent_phase) %>%
    tidyr::pivot_wider(names_from = RNA_combined_celltype,
                       values_from = percent_phase)
  
  openxlsx::addWorksheet(wb = excel_workbook,
                         sheetName = paste0(sample_name, "_percents"))
  
  openxlsx::writeData(wb = excel_workbook,
                      sheet = paste0(sample_name, "_percents"),
                      x = merged_meta)
  
  merged_meta <- merged_meta %>%
    tibble::column_to_rownames("Phase")
  
  max_meta <- data.frame("top_phase" = rownames(merged_meta)[apply(merged_meta,
                                                                   2, which.max)])
  
  rownames(max_meta) <- colnames(merged_meta)
  
  openxlsx::addWorksheet(wb = excel_workbook,
                         sheetName = paste0(sample_name, "_top_phase"))
  
  openxlsx::writeData(wb = excel_workbook,
                      sheet = paste0(sample_name, "_top_phase"),
                      x = max_meta,
                      rowNames = TRUE)  
}

cell_cycle_percent <- openxlsx::createWorkbook()

# All samples
write_cell_cycle_to_excel(sample = "all_cells", seurat_object = merged_seurat,
                          excel_workbook = cell_cycle_percent)

invisible(lapply(all_samples, function(x){
  write_cell_cycle_to_excel(sample = x, seurat_object = merged_seurat,
                            excel_workbook = cell_cycle_percent)
}))

openxlsx::saveWorkbook(wb = cell_cycle_percent,
                       overwrite = TRUE,
                       file = here("results/All_combined/R_analysis/files/cell_cycle_percents.xlsx"))

saveRDS(seurat_data, file.path(save_dir, "rda_obj/seurat_processed.rds"))
