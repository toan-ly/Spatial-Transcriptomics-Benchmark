library(BayesSpace)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(mclust)
library(aricode)
library(clevr)  # For homogeneity, completeness, v-measure
library(cluster) # For silhouette score
library(bench)   # For benchmarking execution time
library(readxl)
library(pryr)
library(peakRAM)
library(scater)  # For runUMAP

batch_cluster_map <- list('-0.04' = 8, '-0.09' = 8, '-0.14' = 8, '-0.24' = 8, '-0.19' = 8)

calculate_metrics <- function(ground_truth, clusters, data_matrix) {
  tryCatch({
    # Remove NA values
    valid_indices <- !is.na(clusters) & !is.na(ground_truth)
    clusters <- clusters[valid_indices]
    ground_truth <- ground_truth[valid_indices]
    data_matrix <- data_matrix[valid_indices, ]

    if (length(ground_truth) != length(clusters)) {
      stop("Mismatch in length between ground_truth and clusters.")
    }
    if (any(is.na(ground_truth)) || any(is.na(clusters))) {
      stop("Missing values found in ground_truth or clusters.")
    }

    # # Convert to factors
    # ground_truth <- as.factor(ground_truth)
    # clusters <- as.factor(clusters)

    ARI <- adjustedRandIndex(ground_truth, clusters)
    AMI <- AMI(ground_truth, clusters)
    homogeneity <- clevr::homogeneity(ground_truth, clusters)
    completeness <- clevr::completeness(ground_truth, clusters)
    v_measure <- clevr::v_measure(ground_truth, clusters)

    # Calculate Silhouette Score (ASW)
    dist_matrix <- dist(data_matrix)
    sil <- silhouette(clusters, dist_matrix)
    ASW <- mean(sil[, 3])

    # Placeholder for CHAOS and PAS
    CHAOS <- NA  
    PAS <- NA    

    return(list(ARI = ARI, 
                AMI = AMI, 
                Homogeneity = homogeneity, 
                Completeness = completeness, 
                V_Measure = v_measure, 
                ASW = ASW, 
                CHAOS = CHAOS, 
                PAS = PAS))
  }, error = function(e) {
    warning("Error calculating metrics: ", e$message)
    return(list(ARI = NA, 
                AMI = NA, 
                Homogeneity = NA, 
                Completeness = NA, 
                V_Measure = NA, 
                ASW = NA, 
                CHAOS = NA, 
                PAS = NA))
  })
}

save_results <- function(data, labels, metrics_df, dir.output) {
  # Save metrics
  write.csv(metrics_df, file = file.path(dir.output, 'metrics.csv'), row.names = FALSE)
  
  # Save spatial clustering plot
  # cluster_plot <- clusterPlot(data, label=labels, palette=NULL, size=0.05) +
  #   scale_fill_viridis_d(option = "A", labels = 1:7, name=NULL) +
  #   labs(title=paste("BayesSpace (ARI =", round(metrics_df$ARI, 2), ")", sep="")) +
  #   theme(plot.title = element_text(hjust = 0.5, size = 16))
  
  # ggsave(file.path(dir.output, 'clustering.pdf'), plot = cluster_plot,
  #        width = 6, height = 6, dpi = 300, device = "pdf")
  
  # Save reduced dimension data
  write.csv(reducedDim(data, "PCA"),
            file = file.path(dir.output, 'low_dim_data.csv'),
            row.names = TRUE)
  
  # Save metadata
  write.csv(colData(data),
            file=file.path(dir.output, 'cell_metadata.csv'),
            row.names = TRUE)
  
  # Save UMAP coordinates and plot
  # umap_coords <- as.data.frame(reducedDim(data, "UMAP_neighbors15"))
  umap_coords <- as.data.frame(reducedDim(data, "UMAP"))
  umap_coords$spot_id <- rownames(umap_coords)
  write.csv(umap_coords,
            file = file.path(dir.output, "spatial_umap_coords.csv"),
            row.names = TRUE)
  
  umap_plot <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = as.factor(labels))) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "BayesSpace", x = "UMAP 1", y = "UMAP 2", color = 'Cluster') +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"))
  
  ggsave(file.path(dir.output, 'umap.pdf'), plot = umap_plot,
         width = 6, height = 6, dpi = 300, device = "pdf")
}


data_path <- file.path("/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/mHypothalamus")
save_path <- file.path("/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/Results/MERFISH/mHypothalamus/BayesSpace")

metrics_list <- list()

for (sample.name in names(batch_cluster_map)) {
  cat("Processing batch:", sample.name, "\n")
  n_clusters <- batch_cluster_map[[sample.name]]

  # dir.input <- file.path(data_path, sample.name)
  dir.output <- file.path(save_path, sample.name)

  if (!dir.exists(file.path(dir.output))) {
    dir.create(file.path(dir.output), recursive = TRUE)
  }


  # benchmark <- mark({
  memory_usage <- peakRAM({
    start_time <- Sys.time()

    set.seed(102)

    # Load MERFISH data

    # load(sprintf("%s/MERFISH_input.RData", data_path))
    # count <- t(RNA)
    # info <- as.data.frame(xyz)
    # colnames(info) <- c("array_row", "array_col", "z")
    # colnames(count) <- NULL

    filename <- paste0(data_path, '/MERFISH_Animal1_cnts.xlsx')
    cnts <- as.data.frame(read_excel(filename, sheet=sample.name))
    row.names(cnts) <- cnts[[1]]
    cnts <- cnts[-1]

    infoname <- paste0(data_path, '/MERFISH_Animal1_info.xlsx')
    info <- as.data.frame(read_excel(infoname, sheet=sample.name))
    row.names(info) <- info[[1]]
    gt_labels <- list(info$z)
    info <- info[-1]
    

    count <- as.matrix(cnts)
    info <- info[c('x', 'y')]
    colnames(info) <- c('array_row', 'array_col')
    colData <- info

    rownames(colData) <- colnames(count)

    sce <- SingleCellExperiment(assays=list(counts=as(count, "dgCMatrix")), colData=colData)

    # colnames(sce) <- make.names(colnames(sce), unique=TRUE)

    sce <- spatialPreprocess(sce, n.PCs = 20, n.HVGs = 2000, log.normalize = TRUE)
  
    q <- n_clusters  
    d <- 20
    
    # colData(sce)$array_row <- colData(sce)$row
    # colData(sce)$array_col <- colData(sce)$col

    data <- spatialCluster(
      sce, 
      q=q, d=d,
      init.method="mclust", model="t", 
      nrep=10000, 
      # gamma=3,
      # save.chain = TRUE,
    )

    data <- runUMAP(data, dimred="PCA", name="UMAP")

    labels <- data$spatial.cluster
    gt <- data$cluster.init
    pca_data <- reducedDim(data, "PCA")

    metrics <- calculate_metrics(gt, labels, pca_data)
    cat("ARI for batch", sample.name, ":", metrics$ARI, "\n")    

    end_time <- Sys.time()
  })
  # }, iterations = 1L)

  execution_time <- as.numeric(end_time - start_time)
  memory_usage <- memory_usage$Peak_RAM_Used_MiB

  metrics_df <- data.frame(
    Sample = sample.name,
    ARI = metrics$ARI,
    AMI = metrics$AMI,
    Homogeneity = metrics$Homogeneity,
    Completeness = metrics$Completeness,
    V_Measure = metrics$V_Measure,
    ASW = metrics$ASW,
    Time = execution_time,
    Memory = memory_usage
  )

  save_results(data, labels, metrics_df, dir.output)
  metrics_list[[sample.name]] <- metrics_df
}

metrics_df <- do.call(rbind, metrics_list)
write.csv(metrics_df, file = file.path(save_path, "metrics.csv"), row.names = TRUE)
