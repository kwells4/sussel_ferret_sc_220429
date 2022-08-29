library(here)
library(openxlsx)

base_dir <- here()

all_samples <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                 "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")


red_style <- openxlsx::createStyle(bgFill = "#cc4125")
light_red_style <- openxlsx::createStyle(bgFill = "#e06666")
orange_style <- openxlsx::createStyle(bgFill = "#f6b26b")
yellow_style <- openxlsx::createStyle(bgFill = "#ffd966")
green_style <- openxlsx::createStyle(bgFill = "#93c47d")


excel_wb <- createWorkbook()

write_samples <- lapply(all_samples, function(sample){
  base_dir_proj <- file.path(base_dir, "results", sample)
  save_dir <- file.path(base_dir_proj, "R_analysis")
  
  
  correlations <- read.csv(file.path(save_dir, "files", "all_clustifyr_cors.csv"))
  
  correlations <- correlations %>%
    tibble::column_to_rownames("X") %>%
    dplyr::select(order(colnames(.))) %>%
    tibble::rownames_to_column("cluster")

  
  openxlsx::addWorksheet(wb = excel_wb,
                         sheetName = sample)
  
  openxlsx::writeData(wb = excel_wb,
                      sheet = sample,
                      x = correlations)
  
  openxlsx::conditionalFormatting(wb = excel_wb,
                                  sheet = sample,
                                  cols = 2:ncol(correlations),
                                  rows = 1:nrow(correlations) + 1,
                                  rule = c(0.5, 1),
                                  type = "between",
                                  style = yellow_style)
  
  openxlsx::conditionalFormatting(wb = excel_wb,
                                  sheet = sample,
                                  cols = 2:ncol(correlations),
                                  rows = 1:nrow(correlations) + 1,
                                  rule = c(0.4, 0.5),
                                  type = "between",
                                  style = light_red_style)
  
  openxlsx::conditionalFormatting(wb = excel_wb,
                                  sheet = sample,
                                  cols = 2:ncol(correlations),
                                  rows = 1:nrow(correlations) + 1,
                                  rule = c(0, 0.4),
                                  type = "between",
                                  style = red_style)
  
  
  
  # Top N per column
  invisible(lapply(1:nrow(correlations), function(x){
    openxlsx::conditionalFormatting(wb = excel_wb,
                                    sheet = sample,
                                    cols = 2:ncol(correlations),
                                    rows = x,
                                    rank = 1,
                                    type = "topN",
                                    style = green_style)  
  }))
  
})


saveWorkbook(wb = excel_wb,
             file = here("results/All_combined/R_analysis/files/correlations.xlsx"),
             overwrite = TRUE)


variable_genes <- lapply(seq_along(all_samples), function(x){
  sample_name <- all_samples[[x]]
  seurat_data <- readRDS(file.path("results", sample_name, "R_analysis",
                                   "rda_obj", "seurat_processed.rds"))
  
  variable_genes <- VariableFeatures(seurat_data)
  
  variable_df <- data.frame(variable_genes)
  colnames(variable_df) <- sample_name
  
  return(variable_df)
})

variable_genes <- do.call(cbind, variable_genes)

write.csv(variable_genes,
          file = here("results/All_combined/R_analysis/files/variable_genes.csv"))
