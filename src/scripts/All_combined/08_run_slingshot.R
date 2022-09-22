library(slingshot)

input_args <- commandArgs(trailingOnly = TRUE)

dimensional_reductions <- input_args[1]
clusters <- input_args[2]
start_cluster <- input_args[3]
save_name <- input_args[4]

dim_red_mat <- read.table(dimensional_reductions, row.names = 1,
	                      sep = "\t")

cluster_mat <- read.table(clusters, row.names = 1,
	                      sep = "\t")


cluster_vector <- cluster_mat[ , 1]
names(cluster_vector) <- rownames(cluster_mat)

slingshot_res <- slingshot(data = dim_red_mat, clusterLabels = cluster_vector,
	                       start.clus = start_cluster)

saveRDS(slingshot_res, save_name)