library(BayesSpace)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(mclust)
library(aricode)
library(clevr)  # Homogeneity, Completeness, V-Measure
library(cluster)  # ASW
library(pryr)  # Memory usage
library(microbenchmark)  # Timing

# Define the number of clusters for the BRCA dataset
batch_cluster_map <- list(
  'V1_Human_Breast_Cancer_Block_A_Section_1' = 20
)

# Paths
data_path <- file.path("/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/BRCA1")
save_path <- file.path("/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/RESULTS/BRCA1/BayesSpace")

# Function to calculate additional metrics
calculate_metrics <- function(ground_truth, clusters) {
  tryCatch({
    valid_indices <- !is.na(clusters) & !is.na(ground_truth)
    clusters <- clusters[valid_indices]
    ground_truth <- ground_truth[valid_indices]

    ARI <- adjustedRandIndex(ground_truth, clusters)
    AMI <- AMI(ground_truth, clusters)
    homogeneity <- clevr::homogeneity(ground_truth, clusters)
    completeness <- clevr::completeness(ground_truth, clusters)
    v_measure <- clevr::v_measure(ground_truth, clusters)

    asw <- silhouette(clusters, dist(reducedDim(dlpfc, "PCA")))[, 3]  # Average Silhouette Width
   

    return(list(ARI = ARI, AMI = AMI, Homogeneity = homogeneity,
                Completeness = completeness, V_Measure = v_measure,
                ASW = mean(asw)))
  }, error = function(e) {
    warning("Error calculating metrics: ", e$message)
    return(rep(NA, 6))
  })
}

# Process each sample
for (sample.name in names(batch_cluster_map)) {
  cat("Processing batch:", sample.name, "\n")
  n_clusters <- batch_cluster_map[[sample.name]]

  dir.input <- file.path(data_path, sample.name)
  dir.output <- file.path(save_path)
  if (!dir.exists(dir.output)) {
    dir.create(dir.output, recursive = TRUE)
  }

  # Load data
  dlpfc <- read10Xh5(dir.input)
  dlpfc <- scuttle::logNormCounts(dlpfc)

  set.seed(101)
  dec <- scran::modelGeneVar(dlpfc)
  top <- scran::getTopHVGs(dec, n = 2000)

  set.seed(102)
  dlpfc <- scater::runPCA(dlpfc, subset_row=top)

  # Preprocessing for BayesSpace
  dlpfc <- spatialPreprocess(dlpfc, platform="Visium", skip.PCA=TRUE)

  q <- n_clusters
  d <- 15

  # Run BayesSpace clustering with timing
  start_time <- Sys.time()
  set.seed(104)
  dlpfc <- spatialCluster(dlpfc, q=q, d=d, platform='Visium', 
                          nrep=10000, gamma=3, save.chain=FALSE)
  end_time <- Sys.time()
  
  elapsed_time <- as.numeric(difftime(end_time, start_time, units="secs"))
  memory_usage <- mem_used()

  labels <- dlpfc$spatial.cluster
  gt <- dlpfc$layer_guess

  metrics <- calculate_metrics(gt, labels)
  metrics$Time <- elapsed_time
  metrics$Memory <- memory_usage / (1024^2)  # Convert to MB

  # Save metrics
  write.csv(as.data.frame(metrics), file = file.path(dir.output, 'clustering_metrics.csv'), row.names = FALSE)

  # Spatial Clustering Plot
  cluster_plot <- clusterPlot(dlpfc, label=labels, palette=NULL, size=0.05) +
    scale_fill_viridis_d(option = "A", labels = as.character(1:max(labels))) +
    labs(title="BayesSpace") +
    theme(plot.title = element_text(hjust = 0.5, size = 16))

  ggsave(file.path(dir.output, 'spatial_clustering.pdf'), plot = cluster_plot,
         width = 6, height = 6, dpi = 300, device = "pdf")

  # Save data as CSV
  write.csv(reducedDim(dlpfc, "PCA"), file = file.path(dir.output, 'low_dim_data.csv'), row.names = TRUE)
  write.csv(colData(dlpfc), file = file.path(dir.output, 'cell_metadata.csv'), row.names = TRUE)

  # UMAP
  umap_coords <- as.data.frame(reducedDim(dlpfc, "UMAP_neighbors15"))
  umap_coords$spot_id <- rownames(umap_coords)
  write.csv(umap_coords, file = file.path(dir.output, "spatial_umap_coords.csv"), row.names = FALSE)

  # UMAP Plot
  umap_plot <- ggplot(umap_coords, aes(x = V1, y = V2, color = as.factor(labels))) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_brewer(palette = "Set1", labels = as.character(1:max(labels))) +
    labs(title = "BayesSpace", x = "UMAP 1", y = "UMAP 2", color = 'Cluster') +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"))

  ggsave(file.path(dir.output, 'umap.pdf'), plot = umap_plot,
         width = 6, height = 6, dpi = 300, device = "pdf")

  # Expression Matrix
  expression_data <- as.data.frame(as.matrix(assay(dlpfc, "counts")))
  write.csv(t(expression_data), file = file.path(dir.output, "expression_matrix.csv"), row.names = TRUE)
}
