library(shiny)
library(tidyverse)
library(cowplot)
library(Seurat)
library(scAnalysisR)

gene_list <- readRDS("all_genes.rds")

# Objects ----------------------------------------------------------------------
seurat_object_list <- c("WT_D2" = "WT_D2.rds",
                        "WT_D5" = "WT_D5.rds",
                        "WT_D7" = "WT_D7.rds",
                        "WT_D9" = "WT_D9.rds",
                        "WT_D14" = "WT_D14.rds",
                        "CFKO_D2" = "CFKO_D2.rds",
                        "CFKO_D5" = "CFKO_D5.rds",
                        "CFKO_D7" = "CFKO_D7.rds",
                        "CFKO_D9" = "CFKO_D9.rds",
                        "CFKO_D14" = "CFKO_D14.rds",
                        "WT_combined" = "WT_combined.rds",
                        "CFKO_combined" = "CFKO_combined.rds",
                        "All_combined" = "All_combined.rds")

# Colors -----------------------------------------------------------------------
byrnes_colors <- c("Acinar" = "#D4405B",
                   "Fev_hi" = "#75B140",
                   "Alpha" = "#885FCF",
                   "Ductal" = "#A5903E",
                   "Ngn3_pos" = "#C65CAC",
                   "Prolif_acinar" = "#55A470",
                   "Prolif_ductal" = "#767FC9",
                   "undetermined" = "#D3D3D3")

baron_colors <- c("ductal" = "#A5903E",
                  "acinar" = "#D4405B")

tabula_muris_colors <- c("pancreatic_acinar_cell" = "#D4405B",
                         "pancreatic_stellate_cell" = "#CC6D3A",
                         "pancreatic_ductal_cell" = "#A5903E",
                         "endothelial_cell" = "#45B0CF",
                         "undetermined" = "#D3D3D3")

qadir_colors <- c("immune_cells" = "#D3D3D3",
                  "progenitor_like_cells" = "#297878",
                  "transitional_to_acinar1" = "#874652",
                  "transitional_to_acinar2" = "#CC1B3B",
                  "centroacinar" = "#78295D",
                  "Small_ducts" = "#B37422")

muraro_colors <- c("ductal" = "#A5903E",
                   "acinar" = "#D4405B",
                   "undetermined" = "#D3D3D3")

krentz_colors <- c("AFP" = "#5C6AE0",
                   "Endo1" = "#45B0CF",
                   "undetermined" = "#D3D3D3")

all_colors <- c("Acinar" = "#D4405B",
                "Ductal" = "#A5903E",
                "ductal" = "#A5903E",
                "Prolif_acinar" = "#55A470",
                "Prolif_ductal" = "#767FC9",
                "progenitor_like_cells" = "#297878",
                "transitional_to_acinar1" = "#874652",
                "transitional_to_acinar2" = "#CC1B3B",
                "centroacinar" = "#78295D")

color_mapping <- list("RNA_celltype_byrnes" = byrnes_colors,
                   "RNA_baron_celltype" = baron_colors,
                   "RNA_tabula_muris_celltype" = tabula_muris_colors,
                   "RNA_baron_human_celltype" = baron_colors,
                   "RNA_qadir_celltype" = qadir_colors,
                   "RNA_muraro_celltype" = muraro_colors,
                   "RNA_krentz_celltype" = krentz_colors,
                   "RNA_combined_celltype" = all_colors)

# Sample colors
sample_colors <- as.character(LaCroixColoR::lacroix_palette("Coconut", 10))
sample_colors[5] <- "#F4E3C7"
names(sample_colors) <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                          "CFKO_D14", "CFKO_D9", "CFKO_D7", "CFKO_D5", "CFKO_D2")
sample_colors <- sample_colors[c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14",
                                 "CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")]


#source("functions.R")

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Gene expression shiny ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Gene expression"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Selector for which dataset ----
      selectInput("seurat_object", "Sample to plot",
                  names(seurat_object_list)),
      
      # Input: Selector for which cell type ----
      selectInput("meta_data", "Meta data column:",
                  c("RNA_celltype_byrnes",
                    "RNA_baron_celltype",
                    "RNA_tabula_muris_celltype",
                    "RNA_baron_human_celltype",
                    "RNA_qadir_celltype",
                    "RNA_muraro_celltype",
                    "RNA_krentz_celltype",
                    "RNA_combined_celltype")),
      
      # Input: Selector for which gene ----
      selectInput("gene", "Gene of interest:",
                  gene_list),
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Plots of the gene on a UMAP and violin plots ----
      plotOutput("gene_plot"),
      plotOutput("cluster_plot")
    )
  )
)

# Server -----------------------------------------------------------------------

# Define server logic to plot various variables ----
server <- function(input, output) {
  
  objects <- reactiveValues(seurat_object = NULL)
  
  observeEvent(input$seurat_object, {
    seurat_object <- readRDS(seurat_object_list[[input$seurat_object]])
    objects$celltype_colors <- color_mapping[[input$meta_data]]
    seurat_object$cluster_celltype <- paste(seurat_object$RNA_cluster,
                                            seurat_object[[input$meta_data]][[1]],
                                            sep = "_")
    
    ncolors <- length(unique(seurat_object$cluster_celltype))
    cluster_colors <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(9, "Set1"))(ncolors)
    
    names(cluster_colors) <- unique(seurat_object$cluster_celltype)
    
    objects$cluster_colors <- cluster_colors
    
    objects$seurat_object <- seurat_object
  })
  
  # Generate plots
  output$gene_plot <- renderPlot({
    
    seurat_object <- objects$seurat_object
    
    cluster_colors <- objects$cluster_colors
    
    celltype_colors <- objects$celltype_colors
    
    umap1 <- plotDimRed(seurat_object, col_by = input$gene,
                        plot_type = "rna.umap")[[1]]
    umap2 <- plotDimRed(seurat_object, col_by = input$meta_data,
                        plot_type = "rna.umap",
                        color = celltype_colors)[[1]]
    violin1 <- featDistPlot(seurat_object, input$gene,
                            sep_by = input$meta_data,
                            combine = FALSE,
                            color = celltype_colors)[[1]]
    
    cowplot::plot_grid(umap1, NULL,
                       umap2, violin1,
                       labels = c("A", "", "B", "C"),
                       nrow = 2, ncol = 2,
                       align = "hv",
                       axis = "lr")

    # umap3 <- plotDimRed(seurat_object, col_by = "cluster_celltype",
    #                     plot_type = "rna.umap",
    #                     color = cluster_colors)[[1]]
    # 
    # violin2 <- featDistPlot(seurat_object, input$gene,
    #                         sep_by = "cluster_celltype",
    #                         combine = FALSE,
    #                         color = cluster_colors)[[1]]
    # 
    # 
    # cowplot::plot_grid(umap1, NULL,
    #                    umap2, violin1,
    #                    umap3, violin2,
    #                    labels = c("A", "", "B", "C", "D", "E"),
    #                    nrow = 3, ncol = 2,
    #                    align = "hv",
    #                    axis = "lr")

  })
  
  # Generate plots
  output$cluster_plot <- renderPlot({
    
    seurat_object <- objects$seurat_object
    
    cluster_colors <- objects$cluster_colors
    
    celltype_colors <- objects$celltype_colors
    

    umap3 <- plotDimRed(seurat_object, col_by = "cluster_celltype",
                        plot_type = "rna.umap",
                        color = cluster_colors)[[1]]
    
    violin2 <- featDistPlot(seurat_object, input$gene,
                            sep_by = "cluster_celltype",
                            combine = FALSE,
                            color = cluster_colors)[[1]]
    
    
    cowplot::plot_grid(umap3, violin2,
                       NULL,
                       labels = c("D", "E", ""),
                       nrow = 2, ncol = 2,
                       align = "hv",
                       axis = "lr")
    
  })
}

shinyApp(ui, server)