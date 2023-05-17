library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
library(ggridges)

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


seurat_name_mapping <- names(all_colors)
names(seurat_name_mapping) <- names(all_colors)

seurat_name_mapping[seurat_name_mapping == "progenitor_like_cells"] = 
  "centroacinar2"

seurat_name_mapping[seurat_name_mapping == "transitional_to_acinar"] = 
  "centroacinar2"

merged_seurat$final_celltype <- seurat_name_mapping[as.character(merged_seurat$RNA_combined_celltype)]

table(merged_seurat$final_celltype, merged_seurat$RNA_combined_celltype)

test_barplot <- scAnalysisR::stacked_barplots(merged_seurat,
                                              meta_col = "RNA_combined_celltype",
                                              split_by = "sample",
                                              return_values = TRUE,
                                              color = all_colors)




new_colors <- c("acinar" = "#D4405B",
                "ductal" = "#A5903E",
                "Prolif_acinar" = "#55A470",
                "Prolif_ductal" = "#767FC9",
                "centroacinar" = "#78295D",
                "centroacinar2" = "#297878")

merged_seurat$final_celltype <- factor(merged_seurat$final_celltype,
                                       levels = names(new_colors))

all_barplots <- scAnalysisR::stacked_barplots(merged_seurat,
                                              meta_col = "final_celltype",
                                              split_by = "sample",
                                              return_values = TRUE,
                                              color = new_colors)


pdf(file.path(sample_dir, "images", "cell_type_plots", "new_cell_type_barplot.pdf"))

print(
  all_barplots$barplot +
    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 0.5, hjust=1))
)
dev.off()


# Differential expression analysis
centroacinar_de <- openxlsx::createWorkbook()

# De between centro acinar populations
test_celltypes <- c("CA" = "centroacinar", 
                    "PL" = "progenitor_like_cells",
                    "TA" = "transitional_to_acinar")

all_results <- lapply(names(sample_colors), function(sample_name){
  seurat_sub <- subset(merged_seurat, subset = sample == sample_name)
  all_celltypes <- unique(seurat_sub$RNA_combined_celltype)
  all_celltypes <- all_celltypes[all_celltypes %in% test_celltypes]
  
  Idents(seurat_sub) <- "RNA_combined_celltype"
  
  if(length(all_celltypes) == 1){
    return(NULL)
  } else {
    all_combinations <- combn(all_celltypes, 2)
    all_de <- lapply(1:ncol(all_combinations), function(x){
      test_1 <- all_combinations[1, x]
      test_2 <- all_combinations[2, x]
      all_markers <- FindMarkers(seurat_sub,
                                 ident.1 = test_1,
                                 ident.2 = test_2) %>%
        dplyr::filter(p_val_adj < 0.05) %>%
        dplyr::mutate(sample = sample_name,
                      comparison = paste(test_1, test_2, sep = ":"))
      
      testing_names <- test_celltypes[test_celltypes %in% c(test_1, test_2)]
      
      save_name <- paste(sample_name,
                         names(testing_names)[[1]],
                         names(testing_names)[[2]],
                         sep = "_")
      
      openxlsx::addWorksheet(wb = centroacinar_de,
                             sheetName = save_name)
      
      openxlsx::writeData(wb = centroacinar_de, sheet = save_name,
                          x = all_markers, rowNames = TRUE,
                          colNames = TRUE)
    })
  }
})

openxlsx::saveWorkbook(wb = centroacinar_de, 
                       file = file.path(sample_dir, "files",
                                        "centroacinar_de.xlsx"),
                       overwrite = TRUE)



all_colors <- c("acinar" = "#D4405B",
                "ductal" = "#A5903E",
                "Prolif_acinar" = "#55A470",
                "Prolif_ductal" = "#767FC9",
                "progenitor_like_cells" = "#297878",
                "transitional_to_acinar" = "#874652",
                "centroacinar" = "#78295D")

seurat_sub <- subset(merged_seurat, 
                     subset = RNA_combined_celltype %in%
                       c("progenitor_like_cells",
                         "transitional_to_acinar",
                         "centroacinar"))

featDistPlot(seurat_sub, "ALDH1A1", sep_by = "sample",
            col_by = "RNA_combined_celltype", combine = FALSE,
            color = all_colors)

featDistPlot(seurat_sub, "ALDH1A1", col_by = "sample",
             sep_by = "RNA_combined_celltype", combine = FALSE,
             color = sample_colors)


all_plots <- featDistPlot(merged_seurat, c("KRT19", "TSPAN8", "ID1"),
                          col_by = "RNA_combined_celltype", 
                          combine = FALSE, color = all_colors, 
                          sep_by = "genotype")
pdf(file.path(sample_dir, "images", "interesting_gene_plots.pdf"))
print(all_plots)

dev.off()
