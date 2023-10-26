library(slingshot)
library(scAnalysisR)
library(Seurat)
library(here)
library(tidyverse)
library(tradeSeq)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

all_samples <- "All_combined_2_9"

all_sample_dir <- here("results", all_samples, "R_analysis")

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

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


# cfko slingshot
cfko_slingshot <- readRDS(file.path(all_sample_dir, "files",
                                    "slingshot", "CFKO_slingshot.rds"))


cfko_pseudotime <- slingPseudotime(cfko_slingshot)
colnames(cfko_pseudotime) <- paste0("cfko_", colnames(cfko_pseudotime))

cfko_pseudotime <- data.frame(cfko_pseudotime)

# Check barcodes
dim(cfko_pseudotime)
length(intersect(rownames(cfko_pseudotime), colnames(merged_seurat)))

merged_seurat <- AddMetaData(merged_seurat, metadata = cfko_pseudotime)

pseudotime_plots <- plotDimRed(merged_seurat, colnames(cfko_pseudotime),
                               plot_type= "rna.umap")



# Lineage 3 = lineage 9
# Lineage 2 = lineage 7
# Lineage 1 kind of is lineage 1, but it's not quite right

# wt slingshot
wt_slingshot <- readRDS(file.path(all_sample_dir, "files",
                                  "slingshot", "WT_slingshot.rds"))


wt_pseudotime <- slingPseudotime(wt_slingshot)
colnames(wt_pseudotime) <- paste0("wt_", colnames(wt_pseudotime))

wt_pseudotime <- data.frame(wt_pseudotime)

# check barcodes
dim(wt_pseudotime)
length(intersect(rownames(wt_pseudotime), colnames(merged_seurat)))

merged_seurat <- AddMetaData(merged_seurat, metadata = wt_pseudotime)

pseudotime_plots <- plotDimRed(merged_seurat, colnames(wt_pseudotime),
                               plot_type= "rna.umap")


# Try out the embedCurves cfko -------------------------------------------------
all_umap <- Embeddings(merged_seurat, reduction = "rna.umap")

cfko_umap <- all_umap[rownames(all_umap) %in% rownames(cfko_pseudotime),]

cfko_umap_slingshot <- embedCurves(cfko_slingshot, newDimRed = cfko_umap)


all_curves <- slingCurves(cfko_umap_slingshot)

all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  
  # This line cuts off the long tail... Probably a better option is
  # to just remove unknown before running slingshot.
  base_plot <- plotDimRed(merged_seurat, col_by = "sample", color = sample_colors,
                          plot_type = "rna.umap",
                          ggrastr = TRUE)[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)  +
    ggplot2::ggtitle(paste0("cfko_lineage", x))
  
  return(base_plot)
})


pdf(file.path(all_sample_dir, "images", "cfko_line_slingshot.pdf"))

print(all_plots)

dev.off()

# I like curves 2, 3, 4, 5, 7, 8

# Plots based on cell type
all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  
  # This line cuts off the long tail... Probably a better option is
  # to just remove unknown before running slingshot.
  base_plot <- plotDimRed(merged_seurat, col_by = "final_ind_celltype",
                          color = all_colors,
                          plot_type = "rna.umap", ggrastr = TRUE)[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)   +
    ggplot2::ggtitle(paste0("cfko_lineage", x))
  
  return(base_plot)
})


pdf(file.path(all_sample_dir, "images", "cfko_line_slingshot_celltype.pdf"))

print(all_plots)

dev.off()


# Plots based on cell type only cells in lineage colored
all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  lineage_name <- paste0("cfko_Lineage", x)
  
  merged_seurat$plot_cells <- ifelse(!is.na(merged_seurat[[lineage_name]][[1]]),
                                     "TRUE", "FALSE")
  
  # merged_seurat$plot_cells <- ifelse(rownames(merged_seurat[[]]) %in%
  #                                               rownames(curve1_coord),
  #                                    "TRUE", "FALSE")
  
  # This line cuts off the long tail... Probably a better option is
  # to just remove unknown before running slingshot.
  base_plot <- plotDimRed(merged_seurat, col_by = "final_ind_celltype",
                          color = all_colors,
                          plot_type = "rna.umap",
                          highlight_group = TRUE,
                          meta_data_col = "plot_cells",
                          group = "TRUE", ggrastr = TRUE)[[1]]
  
  # Make sure the cells used in both are the same
  curve1_coord <- curve1_coord[rownames(curve1_coord) %in% 
                                 rownames(base_plot$data),]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)   +
    ggplot2::ggtitle(paste0("cfko_lineage", x))
  
  
  full_plot <- plotDimRed(merged_seurat, col_by = "final_ind_celltype",
                          color = all_colors,
                          plot_type = "rna.umap",
                          ggrastr = TRUE)[[1]]
  
  # Merge all data together
  lineage_info <- curve1_coord %>%
    dplyr::select(rnaUMAP_1, rnaUMAP_2) %>%
    dplyr::rename(line_umap1 = rnaUMAP_1,
                  line_umap2 = rnaUMAP_2) %>%
    tibble::rownames_to_column("barcode")
  
  umap_info <- base_plot$data %>%
    dplyr::select(dim1, dim2, colour_metric) %>%
    dplyr::rename(umap1 = dim1,
                  umap2 = dim2,
                  celltype = colour_metric) %>%
    tibble::rownames_to_column("barcode")
  
  umap_info <- full_join(umap_info, lineage_info, by = "barcode") %>%
    tibble::column_to_rownames("barcode")
  
  full_plot <- full_plot$data %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::filter(!barcode %in% rownames(umap_info)) %>%
    tibble::column_to_rownames("barcode") %>%
    dplyr::rename(umap1 = dim1,
                  umap2 = dim2) %>%
    dplyr::mutate(celltype = NA, line_umap1 = NA,
                  line_umap2 = NA) %>%
    dplyr::select(umap1, umap2, celltype, line_umap1, line_umap2)
  
  # Merge data frames together
  return_data <- rbind(umap_info, full_plot)
  
  return(list(plot = base_plot,
              data = return_data))
})

all_data_cfko <- lapply(all_plots, function(x){
  return(x$data)
})

all_plots <- lapply(all_plots, function(x){
  return(x$plot)
})

pdf(file.path(all_sample_dir, "images",
              "cfko_line_slingshot_celltype_highlight.pdf"))

print(all_plots)

dev.off()


# Try out the embedCurves wt ---------------------------------------------------
all_umap <- Embeddings(merged_seurat, reduction = "rna.umap")

wt_umap <- all_umap[rownames(all_umap) %in% rownames(wt_pseudotime),]

wt_umap_slingshot <- embedCurves(wt_slingshot, newDimRed = wt_umap)


all_curves <- slingCurves(wt_umap_slingshot)

all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  base_plot <- plotDimRed(merged_seurat, col_by = "sample", color = sample_colors,
                          plot_type = "rna.umap", ggrastr = TRUE)[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)   +
    ggplot2::ggtitle(paste0("wt_lineage", x))
  
  return(base_plot)
})

# I like curves 2, 3, 4, 5, 7, 8

pdf(file.path(all_sample_dir, "images", "wt_line_slingshot.pdf"))

print(all_plots)

dev.off()



# Plots based on cell type
all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  
  # This line cuts off the long tail... Probably a better option is
  # to just remove unknown before running slingshot.
  base_plot <- plotDimRed(merged_seurat, col_by = "final_ind_celltype",
                          color = all_colors,
                          plot_type = "rna.umap", ggrastr = TRUE)[[1]]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)  +
    ggplot2::ggtitle(paste0("wt_lineage", x))
  
  return(base_plot)
})


pdf(file.path(all_sample_dir, "images", "wt_line_slingshot_celltype.pdf"))

print(all_plots)

dev.off()


# Plots based on cell type only cells in lineage colored
all_plots <- lapply(1:length(all_curves), function(x){
  c <- all_curves[[x]]
  curve1_coord <- data.frame(c$s[c$ord, c(1,2)])
  curve1_coord$stage <- "line"
  
  lineage_name <- paste0("wt_Lineage", x)
  
  merged_seurat$plot_cells <- ifelse(!is.na(merged_seurat[[lineage_name]][[1]]),
                                     "TRUE", "FALSE")
  
  base_plot <- plotDimRed(merged_seurat, col_by = "final_ind_celltype",
                          color = all_colors,
                          plot_type = "rna.umap",
                          highlight_group = TRUE,
                          meta_data_col = "plot_cells",
                          group = "TRUE", ggrastr = TRUE)[[1]]
  
  # Make sure the cells used in both are the same
  curve1_coord <- curve1_coord[rownames(curve1_coord) %in% 
                                 rownames(base_plot$data),]
  base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
                                              ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
                                              color = "black", size = 1)   +
    ggplot2::ggtitle(paste0("cfko_lineage", x))
  
  
  full_plot <- plotDimRed(merged_seurat, col_by = "final_ind_celltype",
                          color = all_colors,
                          plot_type = "rna.umap",
                          ggrastr = TRUE)[[1]]
  
  # Merge all data together
  lineage_info <- curve1_coord %>%
    dplyr::select(rnaUMAP_1, rnaUMAP_2) %>%
    dplyr::rename(line_umap1 = rnaUMAP_1,
                  line_umap2 = rnaUMAP_2) %>%
    tibble::rownames_to_column("barcode")
  
  umap_info <- base_plot$data %>%
    dplyr::select(dim1, dim2, colour_metric) %>%
    dplyr::rename(umap1 = dim1,
                  umap2 = dim2,
                  celltype = colour_metric) %>%
    tibble::rownames_to_column("barcode")
  
  umap_info <- full_join(umap_info, lineage_info, by = "barcode") %>%
    tibble::column_to_rownames("barcode")
  
  full_plot <- full_plot$data %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::filter(!barcode %in% rownames(umap_info)) %>%
    tibble::column_to_rownames("barcode") %>%
    dplyr::rename(umap1 = dim1,
                  umap2 = dim2) %>%
    dplyr::mutate(celltype = NA, line_umap1 = NA,
                  line_umap2 = NA) %>%
    dplyr::select(umap1, umap2, celltype, line_umap1, line_umap2)
  
  # Merge data frames together
  return_data <- rbind(umap_info, full_plot)
  
  return(list(plot = base_plot,
              data = return_data))
})

all_data_wt <- lapply(all_plots, function(x){
  return(x$data)
})

all_plots <- lapply(all_plots, function(x){
  return(x$plot)
})


pdf(file.path(all_sample_dir, "images",
              "wt_line_slingshot_celltype_highlight.pdf"))

print(all_plots)

dev.off()

# merged_seurat$plot_cells <- NULL
save_data <- openxlsx::createWorkbook()
invisible(lapply(1:length(all_data_cfko), function(x){
  name <- paste0("lineage_", x, "_cfko")
  data <- all_data_cfko[[x]]
  openxlsx::addWorksheet(wb = save_data, sheetName = name)
  openxlsx::writeData(wb = save_data, sheet = name,
                      x = data, rowNames = TRUE)
}))

invisible(lapply(1:length(all_data_wt), function(x){
  name <- paste0("lineage_", x, "_wt")
  data <- all_data_wt[[x]]
  openxlsx::addWorksheet(wb = save_data, sheetName = name)
  openxlsx::writeData(wb = save_data, sheet = name,
                      x = data, rowNames = TRUE)
}))

openxlsx::saveWorkbook(wb = save_data,
                       file = file.path(all_sample_dir, "files",
                                        "slingshot",
                                        "slingshot_umap_plots.xlsx"),
                       overwrite = TRUE)

saveRDS(merged_seurat, file.path(all_sample_dir, "rda_obj",
                                 "seurat_processed.rds"))