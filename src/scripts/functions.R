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
