library(openxlsx)
library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(viridis)
source("src/scripts/functions.R")

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "All_combined_2_9"

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

all_colors2 <- c("acinar" = "#D4405B",
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

fig_dir <- file.path(save_dir, "images", "final_figures")

ifelse(!dir.exists(fig_dir), dir.create(fig_dir), FALSE)

seurat_data$celltype <- seurat_data$final_ind_celltype

seurat_data$sample <- factor(seurat_data$sample,
                             levels = names(sample_colors))


all_markers <- find_write_markers_orthologs(seurat_object = seurat_data,
                                            save_dir = save_dir,
                                            mapping_file = mapping_file,
                                            mapping_gene_col = "gene_id",
                                            mapping_ortholog_col = c("Mouse.gene.name",
                                                                     "Human.gene.name",
                                                                     "Dog.gene.name",
                                                                     "Pig.gene.name"),
                                            meta_col = "celltype", 
                                            assay = "RNA", pval = 0.05,
                                            logfc = 0.5, gene_lists = NULL,
                                            pairwise = FALSE)


# Fix names
all_markers <- all_markers %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0)
#all_markers_test <- all_markers[51:60,] %>% data.frame()
all_markers_fixed <- lapply(1:nrow(all_markers), function(x){
  row_data <- all_markers[x, ]
  gene_name <- row_data$gene
  if(grepl("LOC", gene_name)){
    # First replace with human name, if that isn't an option, replace with a 
    # different name
    if(!is.na(row_data$Human.gene.name) & row_data$Human.gene.name != ""){
      row_data$gene <- row_data$Human.gene.name
    } else if (!is.na(row_data$Mouse.gene.name) & row_data$Mouse.gene.name != ""){
      row_data$gene <- row_data$Mouse.gene.name
    } else if (!is.na(row_data$Dog.gene.name) & row_data$Dog.gene.name != ""){
      row_data$gene <- row_data$Dog.gene.name
    } else if (!is.na(row_data$Pig.gene.name) & row_data$Pig.gene.name != ""){
      row_data$gene <- row_data$Pig.gene.name
    }
  }
  return(row_data[ , c("gene", "avg_log2FC", "p_val_adj", "cluster")])
})

all_markers_fixed <- do.call(rbind, all_markers_fixed) %>%
  unique()

new_wb <- openxlsx::createWorkbook()

cluster_markers <- lapply(unique(all_markers_fixed$cluster), function(cluster){
  cluster <- as.character(cluster)
  cluster_df <- all_markers_fixed[all_markers_fixed$cluster == cluster,
                                  c("gene", "avg_log2FC", "p_val_adj")]
  colnames(cluster_df) <- c("Gene", "Avg Log2FC", "P val adj")
  openxlsx::addWorksheet(wb = new_wb, sheetName = cluster)
  openxlsx::writeData(wb = new_wb, sheet = cluster, x = cluster_df)
  return(cluster_df)
})

openxlsx::saveWorkbook(wb = new_wb, file = file.path(save_dir, "files",
                                                     "celltype_de.xlsx"),
                       overwrite = TRUE)

names(cluster_markers) <- unique(all_markers_fixed$cluster)

byrnes_markers <- openxlsx::readWorkbook(xlsxFile = file.path("files/byrnes_supp5.xlsx"),
                                         startRow = 2)

acinar_markers <- byrnes_markers %>%
  dplyr::filter(Cluster.ID %in% c("Acinar", "Mature Acinar"))

intersecting_genes_acinar <- intersect(tolower(acinar_markers$Gene.Name),
                                       tolower(cluster_markers$acinar$Gene))

pro_acinar_markers <- byrnes_markers %>%
  dplyr::filter(Cluster.ID %in% c("Proliferating Acinar"))


intersecting_genes_pro_acinar <- intersect(tolower(acinar_markers$Gene.Name),
                                       tolower(cluster_markers$Prolif_acinar$Gene))

intersecting_genes_pro_acinar2 <- intersect(tolower(pro_acinar_markers$Gene.Name),
                                           tolower(cluster_markers$Prolif_acinar$Gene))


ductal_markers <- byrnes_markers %>%
  dplyr::filter(Cluster.ID %in% c("Ductal"))


intersecting_genes_ductal <- intersect(tolower(ductal_markers$Gene.Name),
                                           tolower(cluster_markers$ductal$Gene))


pro_ductal_markers <- byrnes_markers %>%
  dplyr::filter(Cluster.ID %in% c("Proliferating Ductal"))


intersecting_genes_pro_ductal <- intersect(tolower(ductal_markers$Gene.Name),
                                           tolower(cluster_markers$Prolif_ductal$Gene))

intersecting_genes_pro_ductal2 <- intersect(tolower(pro_ductal_markers$Gene.Name),
                                            tolower(cluster_markers$Prolif_ductal$Gene))




# Overlap with the published table 5
# https://www.nature.com/articles/s41467-018-06176-3#MOESM7
# Acinar LDHA, HSPD1, ASNS, NUPR1, FKBP11
# Proliferating acinar LDHA, HSPD1, ASNS, FKBP11
# Ductal ANXA2, HES1, S100A10, CLU, CAPG
# Proliferating ductal ATP1B1, HES1, S100A10

all_markers_fixed[all_markers_fixed$gene %in% 
                    c("LDHA", "HSPD1", "ASNS", "NUPR1", "FKBP11"),] %>%
  data.frame



all_markers_fixed[all_markers_fixed$gene %in% 
                    c("ANXA2", "HES1", "S100A10", "CLU", "CAPG"),] %>%
  data.frame

# Over log fold change 0.5
# Acinar FKBP11, HSPD1, LDHA, NUPR1
# Proliferating acinar LDHA, HSPD1
# Ductal CAPG

test_genes <- c("ADIRF", "CAPG", "S100A10", "KRT5", "KRT15", "KRT19", # Ductal
                "CCNA2", "CCNB1", "CCNB2", "CDC20", "CDCA2", # Cycling
                "HSPD1", "FKBP11", "LDHA", "NUPR1", # Acinar
                "KRT7", "KRT16", "WSB1", "PROM1", "S100A6", "S100A4", # Centroacinar
                "ALDH1A1", "BCL2A1", "KLF4", "TFF3", "TSPAN8", "MUC1", # Centroacinar progenitor
                "MUC13", "MUC16", "MUC4", "MUC5AC", "MUC5B", "MUC20",
                "CEACAM1")


`%notin%` <- Negate(`%in%`)
not_found <- test_genes[test_genes %notin% rownames(seurat_data)]

mapping_file[mapping_file$Human.gene.name %in% not_found,
             c("Human.gene.name", "gene_id")]


test_genes <- c("HSPD1", "FKBP11", "LDHA", "NUPR1", "ASNS", "PRDX2", "PRDX4", 
                "FKBP4", "FKBP11", "NUPR1", # Acinar
                "CCNA2", "CCNB1", "CCNB2", "CDC20", "CDCA2", # Cycling
                "ADIRF", "CAPG", "S100A10", "KRT5", "KRT15", "KRT19", # Ductal
                "KRT7", "LOC101679442", "WSB1", "PROM1", "S100A6", "S100A4", # Centroacinar
                "ALDH1A1", "BCL2A1", "KLF4", "TFF3", "TSPAN8", "MUC1", # Centroacinar progenitor
                "MUC13", "MUC16", "MUC4", "MUC5AC", "MUC5B", "LOC101686708",
                "LOC101682302")

seurat_data$celltype <- factor(seurat_data$celltype, 
                               levels = c("acinar", "Prolif_acinar",
                                          "ductal", "Prolif_ductal",
                                          "centroacinar",
                                          "centroacinar_progenitor"))



# Full heatmap
heatmap_all <- plot_heatmap(seurat_data, gene_list = test_genes,
                            meta_col = "celltype",
                            colors = all_colors2, average_expression = TRUE,
                            return_data = TRUE)

all_z_scores <- heatmap_all$z_score

# Mapping
# LOC101679442 = KRT16
# LOC101686708 = MUC20
# LOC101682302 = CEACAM

rownames(all_z_scores)[rownames(all_z_scores) == "LOC101679442"] <- "KRT16"
rownames(all_z_scores)[rownames(all_z_scores) == "LOC101686708"] <- "MUC20"
rownames(all_z_scores)[rownames(all_z_scores) == "LOC101682302"] <- "CEACAM"

sample_info <- seurat_data[["celltype"]] %>%
  dplyr::distinct() %>%
  dplyr::arrange(celltype)
rownames(sample_info) <- sample_info$celltype

blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", 
                "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")

coloring <- list(celltype = all_colors2)


pdf(file.path(save_dir, "images", "final_figures", "marker_genes.pdf"))
pheatmap::pheatmap(all_z_scores, cluster_rows = FALSE, 
                   cluster_cols = FALSE, show_rownames = TRUE, 
                   show_colnames = FALSE, annotation_col = sample_info, 
                   annotation_colors = coloring, color = blueYellow, border_color = NA, 
                   clustering_method = "complete")

dev.off()

# Make heatmaps separately
seurat_wt <- subset(seurat_data, subset = genotype == "WT")
heatmap_wt <- plot_heatmap(seurat_wt, gene_list = test_genes,
                           meta_col = "celltype",
                           colors = all_colors2, average_expression = TRUE,
                           return_data = TRUE)

wt_z_scores <- heatmap_wt$z_score

seurat_cf <- subset(seurat_data, subset = genotype == "CFKO")
seurat_cf$celltype <- factor(seurat_data$celltype, 
                               levels = c("Prolif_acinar",
                                          "ductal",
                                          "centroacinar",
                                          "centroacinar_progenitor"))

heatmap_cf <- plot_heatmap(seurat_cf, gene_list = test_genes,
                           meta_col = "celltype",
                           colors = all_colors2, average_expression = TRUE,
                           return_data = TRUE)

cf_z_scores <- heatmap_cf$z_score


# Combine heatmaps
colnames(wt_z_scores) <- paste0("wt_", colnames(wt_z_scores))
colnames(cf_z_scores) <- paste0("cf_", colnames(cf_z_scores))

all_z_scores <- cbind(wt_z_scores, cf_z_score)


# Rename the three "LOC"