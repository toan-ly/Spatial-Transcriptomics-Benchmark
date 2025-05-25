# Import necessary libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(mclust)
library(aricode)
library(clevr)  # For homogeneity, completeness, v-measure
options(bitmapType = 'cairo')
library(cluster) # For silhouette score
library(bench)   # For benchmarking execution time
library(peakRAM()) # For memory usage tracking
library(readxl)


batch_cluster_map <- list("-0.04" = 8, "-0.09" = 8, "-0.14" = 8, "-0.24" = 8, "-0.19" = 8)

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

    # Convert factors to numeric for clustering metrics
    if (is.factor(clusters)) {
      clusters <- as.numeric(as.character(clusters))
    }
    if (is.factor(ground_truth)) {
      ground_truth <- as.numeric(as.character(ground_truth))
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

load_dataset <- function(input_path, sample.name, cluster.number) {
    dir.input <- file.path(input_path, sample.name)
    
    filename = paste0(input_path, '/MERFISH_Animal1_cnts.xlsx')
    cnts <- as.data.frame(read_excel(filename, sheet = sample.name))
    row.names(cnts) <- cnts[[1]]
    cnts <- cnts[-1]

    infoname <- paste0(input_path, '/MERFISH_Animal1_info.xlsx')
    info <- as.data.frame(read_excel(infoname, sheet = sample.name))
    row.names(info) <- info[[1]]
    info <- info[-1]

    sp_data <- CreateSeuratObject(
      counts = cnts, 
      project = 'MERFISH',
      min.cells = 3,
      names.delim = "-",
      names.field = 2)

    sp_data <- AddMetaData(sp_data, metadata = info$x, col.name = 'row')
    sp_data <- AddMetaData(sp_data, metadata = info$y, col.name = 'col')
    sp_data <- AddMetaData(sp_data, metadata = info$z, col.name = 'layer_guess')

    # sp_data$orig.ident <- 1

    return(sp_data)
}

save_results <- function(sp_data, metrics_df, dir.output) {
  # Save metrics
  write.csv(metrics_df, file = file.path(dir.output, 'metrics.csv'), row.names = FALSE)
  
  # Generate and save clustering plots separately
  p1 <- DimPlot(sp_data, reduction = "umap", label = TRUE) + 
    ggtitle(paste("Seurat")) + 
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          legend.title = element_blank())  
          
  # p2 <- SpatialDimPlot(sp_data, label = TRUE, label.size = 3) + 
  #   ggtitle(paste("Seurat (ARI =", round(metrics_df$ARI[1], 2), ")")) + 
  #   theme(plot.title = element_text(hjust = 0.5, size = 16),
  #         legend.title = element_blank())
  
  
  ggsave(file.path(dir.output, 'umap.pdf'), 
         plot = p1, width = 6, height = 6, dpi = 300)
  # ggsave(file.path(dir.output, 'clustering.pdf'), 
  #        plot = p2, width = 6, height = 6, dpi = 300)
  
  # Save reduced dimension data
  write.csv(sp_data@reductions$pca@cell.embeddings,
            file = file.path(dir.output, 'low_dim_data.csv'),
            row.names = TRUE)
  
  # Save metadata
  write.csv(sp_data@meta.data,
            file = file.path(dir.output, 'cell_metadata.csv'),
            row.names = TRUE)
  
  # Save UMAP coordinates
  umap_coords <- as.data.frame(sp_data@reductions$umap@cell.embeddings)
  umap_coords$spot_id <- rownames(umap_coords)
  write.csv(umap_coords,
            file = file.path(dir.output, "spatial_umap_coords.csv"),
            row.names = TRUE)
  
}

data_path <- file.path("/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/mHypothalamus")
save_path <- file.path("/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/Results/MERFISH/mHypothalamus/Seurat")
metrics_list <- list()

# Loop through all batches
for (sample.name in names(batch_cluster_map)) {
  cat("Processing batch:", sample.name, "\n")
  n_clusters <- batch_cluster_map[[sample.name]]

  # Define input and output directories
  # dir.input <- file.path(data_path, sample.name)
  dir.output <- file.path(save_path, sample.name)
  if (!dir.exists(dir.output)) {
    dir.create(dir.output, recursive = TRUE)
  }

  benchmark <- mark({
    # # Load the dataset
    # load(sprintf("%s/MERFISH_Animal1.RData", data_path))
    
    # cnts <- cnts_mult[sample.name] # a list of gene expression count matrices
    # cnts_all <- do.call(cbind, cnts)

    # sp_data <- CreateSeuratObject(
    #   counts = cnts_all,
    # )

    sp_data <- load_dataset(data_path, sample.name, n_clusters)

   
    # # Data normalization using SCTransform
    # sp_data <- SCTransform(sp_data, assay = "Spatial", verbose = FALSE)

    sp_data <- NormalizeData(sp_data)
    sp_data <- ScaleData(sp_data, features = rownames(sp_data))

    # Dimensionality reduction and clustering
    # sp_data <- RunPCA(sp_data, assay = "SCT", verbose = FALSE, npcs = 50)
    sp_data <- RunPCA(sp_data, features = rownames(sp_data), verbose = FALSE)
    sp_data <- FindNeighbors(sp_data, reduction = "pca", dims = 1:30)

    # Find optimal resolution for clustering
    sp_data <- tryCatch({
      for (resolution in seq(1, 0.01, by = -0.01)) {
        sp_data <- FindClusters(sp_data, verbose = FALSE, resolution = resolution)
        if (length(levels(sp_data@meta.data$seurat_clusters)) == n_clusters) {
          cat("Optimal resolution found for batch", sample.name, ": ", resolution, "\n")
          break
        }
      }
      sp_data
    }, error = function(e) {
      stop("Error during clustering for batch ", sample.name, ": ", e$message)
    })

    # Run UMAP for visualization
    sp_data <- RunUMAP(sp_data, reduction = "pca", dims = 1:30)

    # Evaluate clustering performance
    gt <- sp_data@meta.data$layer_guess
    pred <- sp_data@meta.data$seurat_clusters
    pca_data <- sp_data@reductions$pca@cell.embeddings
    metrics <- calculate_metrics(gt, pred, pca_data)
    cat("ARI for batch", sample.name, ":", metrics$ARI, "\n")
  }, iterations = 1L)

  execution_time <- as.numeric(benchmark$time[[1]])  # Time in seconds
  memory_usage <- as.numeric(benchmark$mem_alloc[[1]]) / (1024^2)  # Memory in MB

  # Create metrics dataframe
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

  # Save results
  save_results(sp_data, metrics_df, dir.output)
  metrics_list[[sample.name]] <- metrics_df
}

metrics_df <- do.call(rbind, metrics_list)
write.csv(metrics_df, file = file.path(save_path, "metrics.csv"), row.names = TRUE)
cat("All batches processed successfully.\n")
