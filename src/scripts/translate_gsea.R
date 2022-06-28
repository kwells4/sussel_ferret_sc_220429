library(openxlsx)
library(here)
library(tidyverse)

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

gene_path <- here("files/GSEA_signalingpathways_EP.xlsx")

new_gene_wb <- openxlsx::createWorkbook()

all_sheets <- openxlsx::getSheetNames(gene_path)

lapply(all_sheets, function(x){
  print(x)
  gene_list <- openxlsx::readWorkbook(gene_path, sheet = x)[[1]]
  total_genes <- length(gene_list)
  ferret_overlap <- length(intersect(gene_list, mapping_file$gene_id))
  human_overlap <- length(intersect(gene_list, mapping_file$Human.gene.name))
  mouse_overlap <- length(intersect(gene_list, mapping_file$Mouse.gene.name))
  
  overlaps <- c("Mouse.gene.name" = mouse_overlap,
                "Human.gene.name" = human_overlap)
  
  best_overlap <- names(overlaps)[overlaps == max(overlaps)]
  
  keep_list <- c("gene_id", best_overlap)
  
  new_list <- mapping_file %>%
    dplyr::select(all_of(keep_list)) %>%
    dplyr::filter(!!as.name(best_overlap) %in% gene_list) %>%
    dplyr::distinct()
  
  colnames(new_list) <- c("gene_id", "ortholog_id")
  
  openxlsx::addWorksheet(wb = new_gene_wb, sheetName = x)
  openxlsx::writeData(wb = new_gene_wb, sheet = x, x = new_list)
  
  
  # diff_list <- setdiff(gene_list, mapping_file[[best_overlap]])
  # 
  # diff_list <- toupper(diff_list)
  # 
  # length(intersect(diff_list, mapping_file$gene_id))
  
})

openxlsx::saveWorkbook(wb = new_gene_wb,
                       file = "files/GSEA_signaling_pathways_with_orthologs.xlsx",
                       overwrite = TRUE)