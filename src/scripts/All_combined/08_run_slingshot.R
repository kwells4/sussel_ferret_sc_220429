library(slingshot)
library(tidyverse)
library(Seurat)
library(tradeSeq)

# Follows two tutorials
# https://statomics.github.io/tradeSeq/articles/tradeSeq.html
# https://statomics.github.io/tradeSeq/articles/fitGAM.html

print("Done laoding packages!!!")

input_args <- commandArgs(trailingOnly = TRUE)

dimensional_reductions <- input_args[1]
clusters <- input_args[2]
start_cluster <- str_split(input_args[3], pattern = ",")[[1]]
seurat_object_path <- input_args[4]
save_name <- input_args[5]
nknots <- input_args[6]

nknots <- as.integer(nknots)

slingshot_file <- paste0(save_name, "_slingshot.rds")

# If slingshot was already run for the file, just read in the object
# otherwise run slingshot
if(file.exists(slingshot_file)){
	slingshot_res <- readRDS(slingshot_file)
} else {
	print("Using starting cluster: ")
	print(start_cluster)

	dim_red_mat <- read.table(dimensional_reductions, row.names = 1,
		                      sep = "\t")

	cluster_mat <- read.table(clusters, row.names = 1,
		                      sep = "\t")

	cluster_vector <- as.character(cluster_mat[ , 1])

	names(cluster_vector) <- rownames(cluster_mat)

	cluster_vector <- factor(cluster_vector)

	print("Starting slingshot analysis!!!")

	slingshot_res <- slingshot(data = dim_red_mat, clusterLabels = cluster_vector,
		                       start.clus = start_cluster)

	saveRDS(slingshot_res, slingshot_file)
}

icmat_file <- paste0(save_name, "_icmat.rds")

seurat_object <- readRDS(seurat_object_path)

slingshot_cells <- cellnames(slingshot_res)

seurat_object <- subset(seurat_object, cells = slingshot_cells)

# Find variable genes, we will only use these
seurat_object <- FindVariableFeatures(seurat_object, nfeatures = 100)

counts <- GetAssayData(seurat_object, slot = "data")

counts <- counts[rownames(counts) %in% VariableFeatures(seurat_object),]

if(file.exists(icmat_file)){
	icMat <- readRDS(icmat_file)

} else {
	print("Finding genes!")

	# Fit negative binomial
	# Figure out knots
	set.seed(0)
	pdf(paste0(save_name, "evaluateK_plots.pdf"))
	icMat <- evaluateK(counts = counts, sds = slingshot_res, k = 3:10, 
	                   nGenes = 200, verbose = T)

	dev.off()
	
	saveRDS(icMat, icmat_file)
}

# Run fitGAM - note for future implementations (ie when accounting for batch,
# covariates can be included here.)
print("Running fitGAM!")
set.seed(0)

pseudotime <- slingPseudotime(slingshot_res, na = FALSE)
cellWeights <- slingCurveWeights(slingshot_res)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
	          nknots = nknots, verbose = FALSE)

fitgam_save <- paste0(save_name, "_fitgam_sce.rda")
saveRDS(sce, fitgam_save)
