# Document information
# This document adds the gene description to the DE excel file

library(tidyverse)
library(openxlsx)
library(here)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "All_combined_2_9"

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

gene_description <- read.table(here("files/gene_product_mapping.tsv"),
                               sep = "\t", header = TRUE, fill = TRUE)

gene_description <- rename(gene_description, gene = gene_id)

days <- c("D2", "D5", "D7", "D9")

de_directory <- file.path(save_dir, "files", "DE")

for(day in days){
  # Grab out results for that test
  de_test <- paste0(day, "_combined_celltype")
  
  excel_file <- file.path(de_directory, paste0(de_test, ".xlsx"))
  
  excel_sheets <- openxlsx::getSheetNames(excel_file)
  
  excel_sheets <- excel_sheets[!grepl("gse", excel_sheets)]
  
  original_workbook <- openxlsx::loadWorkbook(excel_file)
  for(sheet in excel_sheets){
    de_df <- openxlsx::readWorkbook(excel_file, sheet = sheet)
    full_de_df <- merge(de_df, gene_description, by = "gene",
                        all.x = TRUE, all.y = FALSE)
    openxlsx::writeData(wb = original_workbook,
                        sheet = sheet, x = full_de_df)
  }
  
  openxlsx::saveWorkbook(wb = original_workbook, file = excel_file, overwrite = TRUE)
}


