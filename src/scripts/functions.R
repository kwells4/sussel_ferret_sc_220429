clustifyr_orthologs <- function (seurat_object, ref_mat, save_dir,
                                 save_name, mapping_file,
                                 mapping_gene_col, mapping_ortholog_col,
                                 assay = "RNA", nfeatures = 1000,
                                 clusters = "RNA_cluster", 
                                 plot_type = "rna.umap", cor_cutoff = 0.5) {
  ifelse(!dir.exists(file.path(save_dir, "images")), 
         dir.create(file.path(save_dir, "images")), FALSE)
  ifelse(!dir.exists(file.path(save_dir, "files")), 
         dir.create(file.path(save_dir,  "files")), FALSE)
  seurat_var <- FindVariableFeatures(seurat_object, assay = assay, 
                                     selection.method = "vst",
                                     nfeatures = nfeatures)
  
  # Keep only genes in both
  seurat_mat <- GetAssayData(object = seurat_var, slot = "data")
  seurat_mat <- seurat_mat[rownames(seurat_mat) %in%
                             mapping_file[[mapping_gene_col]],]
  
  # keep only variable genes
  seurat_mat <- seurat_mat[rownames(seurat_mat) %in% 
                             VariableFeatures(seurat_var),]
  
  if(nrow(seurat_mat) < 500){
    warning(paste0("After filtering for genes in your ortholog set only ",
                   nrow(seurat_mat), " gene remain!"))
  }
  
  # Map genes
  all_genes <- data.frame("genes" = rownames(seurat_mat)) %>%
    dplyr::filter(genes %in% dplyr::all_of(mapping_file[[mapping_gene_col]]))
  
  mapping_list <- mapping_file[[mapping_ortholog_col]]
  names(mapping_list) <- mapping_file[[mapping_gene_col]]
  
  all_genes$orthologs <- mapping_list[all_genes$genes]
  
  if(!identical(all_genes$genes, rownames(seurat_mat))){
    stop("Not sure how to handle cases where the row numbers don't add up")
  }
  
  # Change gene names based on the mapping
  rownames(seurat_mat) <- all_genes$orthologs
  
  ref_mat <- ref_mat[rownames(ref_mat) %in% rownames(seurat_mat), ]
  return_list <- list()
  DefaultAssay(seurat_object) <- seurat_assay
  seurat_metadata <- seurat_object[[clusters]]
  seurat_metadata[[clusters]] <- as.character(seurat_metadata[[clusters]])
  seurat_object[[clusters]] <- seurat_metadata[[clusters]]
  
  seurat_genes <- rownames(seurat_mat)
  
  seurat_res <- clustify(input = seurat_mat, metadata = seurat_metadata, 
                         ref_mat = ref_mat, query_genes = seurat_genes,
                         cluster_col = clusters)
  
  return_list$RNA <- seurat_res
  pheatmap::pheatmap(seurat_res, color = viridisLite::viridis(10))
  seurat_cluster <- cor_to_call(seurat_res) %>%
    mutate(type = ifelse(r < cor_cutoff, "undetermined", type))
  new_clusters <- seurat_cluster$type
  names(new_clusters) <- seurat_cluster$cluster
  colname <- paste0(assay, "_", save_name)
  seurat_object[[colname]] <- new_clusters[seurat_object[[clusters]][[1]]]
  plot1 <- plotDimRed(seurat_object, col_by = colname, plot_type = plot_type)
  
  pdf(file.path(save_dir, "images", paste0("RNA_mapping_", 
                                           save_name, ".pdf")))
  print(plot1)
  dev.off()
  return_list$object <- seurat_object
  return(return_list)
}

find_write_markers_orthologs <- function(seurat_object, save_dir,
                                         mapping_file = NULL,
                                         mapping_gene_col = NULL,
                                         mapping_ortholog_col = NULL,
                                         meta_col = "RNA_cluster",
                                         assay = "RNA", pval = 0.05,
                                         logfc = 0.5, gene_lists = NULL,
                                         pairwise = FALSE) {
  
  ifelse(!dir.exists(file.path(save_dir, "files", "DE")),
         dir.create(file.path(save_dir, "files", "DE")), FALSE)
  
  Idents(seurat_object) <- meta_col
  
  if (pairwise) {
    marker_genes <- find_write_markers_pairwise_ortj(seurat_object = seurat_object, 
                                                save_dir = save_dir,
                                                meta_col = meta_col,
                                                assay = assay, 
                                                pval = pval, logfc = logfc,
                                                gene_lists = gene_lists,
                                                mapping_file = mapping_file,
                                                mapping_gene_col = mapping_gene_col,
                                                mapping_ortholog_col = mapping_ortholog_col)
  } else {
    marker_genes <- find_write_markers_orig_orth(seurat_object = seurat_object, 
                                            save_dir = save_dir,
                                            meta_col = meta_col, assay = assay, 
                                            pval = pval, logfc = logfc,
                                            gene_lists = gene_lists,
                                            mapping_file = mapping_file,
                                            mapping_gene_col = mapping_gene_col,
                                            mapping_ortholog_col = mapping_ortholog_col)
  }
  return(marker_genes)
}

find_write_markers_orig_orth <- function(seurat_object, save_dir,
                                         meta_col = "RNA_cluster", 
                                         mapping_file = NULL,
                                         mapping_gene_col = NULL,
                                         mapping_ortholog_col = NULL,
                                         assay = "RNA", pval = 0.05,
                                         logfc = 0.5, gene_lists = NULL) 
{
  marker_genes <- FindAllMarkers(seurat_object, assay = assay, 
                                 only.pos = TRUE)
  if(!is.null(mapping_file)){
    # Add in ortholog gene
    new_mapping <- mapping_file %>%
      dplyr::select(all_of(c(mapping_gene_col, mapping_ortholog_col)))
    
    marker_genes <- merge(marker_genes, new_mapping, by.x = "gene",
                          by.y = mapping_gene_col, all.x = TRUE,
                          all.y = FALSE) %>%
      distinct() %>%
      group_by(cluster) %>%
      arrange(p_val_adj, .by_group = TRUE)
  }
  
  write.csv(marker_genes, file = file.path(save_dir, "files", 
                                           "DE", paste0(assay, "_markers_",
                                                        meta_col, ".csv")))
  gene_wb <- createWorkbook()
  
  full_list <- lapply(unique(marker_genes$cluster), function(x) {
    x <- as.character(x)
    new_df <- marker_genes %>% dplyr::filter(cluster == x & 
                                               p_val_adj < pval &
                                               avg_log2FC > logfc)
    addWorksheet(gene_wb, x)
    writeData(gene_wb, x, new_df)
  })
  
  if (!is.null(gene_lists)) {
    hypergeometric <- hypergeometric_test(seurat_object = seurat_data, 
                                          gene_list = gene_lists,
                                          DE_table = marker_genes, 
                                          DE_p_cutoff = 0.05,
                                          DE_lfc_cutoff = 0.5,
                                          correction_method = "fdr")
    write.csv(hypergeometric, file = file.path(save_dir, 
                                               "files", "DE",
                                               paste0(assay, "_hypergeometric_", 
                                                      meta_col, ".csv")))
    
    full_list <- lapply(unique(hypergeometric$cluster), function(x) {
      x <- as.character(x)
      new_df <- hypergeometric %>% dplyr::filter(cluster == 
                                                   x)
      worksheet_name <- paste0(x, "_gse")
      addWorksheet(gene_wb, worksheet_name)
      writeData(gene_wb, worksheet_name, new_df)
    })
  }
  saveWorkbook(gene_wb, file = file.path(save_dir, "files", 
                                         "DE", paste0(assay, "_markers_",
                                                      meta_col, ".xlsx")), 
               overwrite = TRUE)
  return(marker_genes)
}

find_write_markers_pairwise_ortj <- function(seurat_object, save_dir,
                                             mapping_file = NULL,
                                             mapping_gene_col = NULL,
                                             mapping_ortholog_col = NULL,
                                             meta_col = "RNA_cluster",
                                             assay = "RNA", pval = 0.05,
                                             logfc = 0.5, gene_lists = NULL) {
  marker_genes <- pairwise_markers(seurat_object, assay = seurat_assay, 
                                   meta_col = meta_col)
  
  if(!is.null(mapping_file)){
    # Add in ortholog gene
    new_mapping <- mapping_file %>%
      dplyr::select(all_of(c(mapping_gene_col, mapping_ortholog_col)))
    
    marker_genes <- merge(marker_genes, new_mapping, by.x = "gene",
                          by.y = mapping_gene_col, all.x = TRUE,
                          all.y = FALSE) %>%
      distinct() %>%
      group_by(cluster) %>%
      arrange(p_val_adj, .by_group = TRUE)
  }
  
  write.csv(marker_genes, file = file.path(save_dir, "files", 
                                           "DE",
                                           paste0(assay,
                                                  "_pairwise_markers_",
                                                  meta_col, ".csv")))
  
  gene_wb <- createWorkbook()
  values <- unique(c(marker_genes$cluster_down, marker_genes$cluster_up))
  combinations <- combn(unique(c(marker_genes$cluster_down, 
                                 marker_genes$cluster_up)), m = 2)
  full_list <- lapply(1:ncol(combinations), function(x) {
    ident1 <- combinations[1, x]
    ident2 <- combinations[2, x]
    sheet_name <- paste0(ident1, "_vs_", ident2)
    new_df <- marker_genes %>% 
      dplyr::filter((cluster_up == ident1 & 
                       cluster_down == ident2) | 
                      (cluster_up == ident2 & cluster_down == ident1)) %>% 
      dplyr::filter(p_val_adj < pval & avg_log2FC > logfc)
    addWorksheet(gene_wb, sheet_name)
    writeData(gene_wb, sheet_name, new_df)
  })
  marker_genes$cluster <- paste0(marker_genes$cluster_up, "_vs_", 
                                 marker_genes$cluster_down)
  if (!is.null(gene_lists)) {
    hypergeometric <- hypergeometric_test(seurat_object = seurat_data, 
                                          gene_list = gene_lists,
                                          DE_table = marker_genes, 
                                          DE_p_cutoff = 0.05,
                                          DE_lfc_cutoff = 0.5, correction_method = "fdr")
    write.csv(hypergeometric, file = file.path(save_dir, 
                                               "files", "DE", paste0(assay,
                                                                     "_hypergeometric_", 
                                                                     meta_col,
                                                                     ".csv")))
    full_list <- lapply(unique(hypergeometric$cluster), function(x) {
      x <- as.character(x)
      new_df <- hypergeometric %>% dplyr::filter(cluster == 
                                                   x)
      worksheet_name <- paste0(x, "_gse")
      addWorksheet(gene_wb, worksheet_name)
      writeData(gene_wb, worksheet_name, new_df)
    })
  }
  saveWorkbook(gene_wb, file = file.path(save_dir, "files", 
                                         "DE", paste0(assay,
                                                      "_pairwise_markers_",
                                                      meta_col, ".xlsx")), 
               overwrite = TRUE)
  
  return(marker_genes)
}

make_plots <- function(seurat_object, cluster_name, celltype_name,
                       cluster_colors, celltype_colors, gene,
                       assay = "RNA", plot_type = "rna.umap",
                       save_name = NULL){
  
  violin1 <- featDistPlot(seurat_object, gene, sep_by = celltype_name,
                          col_by = celltype_name, color = celltype_colors,
                          assay = assay)
  
  
  umap1 <- plotDimRed(seurat_object, col_by = gene, plot_type = plot_type,
                      assay = assay)[[1]]
  
  umap3 <- plotDimRed(seurat_object, col_by = celltype_name,
                      plot_type = plot_type, color = celltype_colors)[[1]]
  
  umap2 <- NULL
  umap4 <- plotDimRed(seurat_object, col_by = cluster_name,
                      plot_type = plot_type, color = cluster_colors)[[1]]
  violin2 <- featDistPlot(seurat_object, gene, sep_by = cluster_name,
                          col_by = cluster_name, color = cluster_colors,
                          assay = assay)
  text_labels <- c("A", "", "B", "C", "D", "E")
  
  
  full_plot <- plot_grid(umap1, umap2,
                         umap3, violin1,
                         umap4, violin2,
                         labels = text_labels,
                         nrow = 3,
                         ncol = 2,
                         align = "hv",
                         axis = "l")
  
  if(!is.null(save_name)){
    pdf(save_name, height = 10, width = 10)
    print(full_plot)
    dev.off()
  }
  return(full_plot)
  
}
