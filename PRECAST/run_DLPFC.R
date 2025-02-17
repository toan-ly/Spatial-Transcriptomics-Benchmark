library(PRECAST)
library(Seurat)
library(ggplot2)
library(bench)
library(clevr)  # For homogeneity, completeness, v-measure
library(cluster) # For silhouette score
library(aricode)
library(mclust)
library(dplyr)
library(BayesSpace)

batch_cluster_map <- list(
  '151669' = 5, '151670' = 5, '151671' = 5, '151672' = 5,
  '151673' = 7, '151674' = 7, '151675' = 7, '151676' = 7,
  '151507' = 7, '151508' = 7, '151509' = 7, '151510' = 7
)

#' Calculate clustering evaluation metrics
#' @param ground_truth Vector of ground truth labels
#' @param clusters Vector of predicted cluster labels
#' @param data_matrix Matrix of data points
#' @return List of evaluation metrics
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
    output_path <- '/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/results3/DLPFC/PRECAST'
    dir.input <- file.path(input_path, sample.name)
    dir.output <- file.path(output_path, sample.name, '/')
    meta.input <- file.path(input_path, sample.name, 'gt')

    layer.input <- file.path(input_path, sample.name, 'gt/layered')

    if(!dir.exists(file.path(dir.output))){
        dir.create(file.path(dir.output), recursive = TRUE)
    }

    filename <- paste0(sample.name, "_filtered_feature_bc_matrix.h5")
    sp_data <- Load10X_Spatial(dir.input, filename = filename)


    df_meta <- read.table(file.path(meta.input, 'tissue_positions_list_GTs.txt'))


    original_row_names <- row.names(df_meta)
    split_data <- strsplit(df_meta$V1, split = ",")
    df_meta <- do.call(rbind, lapply(split_data, function(x) {
    data.frame(V1=x[1], V2=x[2], V3=x[3], V4=x[4], V5=x[5], V6=x[6], V7=x[7])
    }))
    row.names(df_meta) <- df_meta$V1
    df_meta$V3 <- as.numeric(df_meta$V3)
    df_meta$V4 <- as.numeric(df_meta$V4)

    # Set the row names of df_meta_matched to be V1
    # Identify the cells that are in both sp_data and df_meta
    common_cells <- colnames(sp_data[["Spatial"]]) %in% rownames(df_meta)

    # Subset sp_data to keep only these cells
    sp_data <- sp_data[, common_cells]

    # Initialize an empty dataframe to hold the final results
    layer.data <- data.frame()

    if(as.numeric(cluster.number) == 5) {
        for(i in 3:6){
            file.name <- paste0(sample.name, "_L", i, "_barcodes.txt")
            file.path <- file.path(layer.input, file.name)

            data.temp <- read.table(file.path, header = FALSE, stringsAsFactors = FALSE) # assuming the file has no header
            data.temp <- data.frame(barcode = data.temp[,1], layer = paste0("layer", i), row.names = data.temp[,1])

            # Append to the final dataframe
            layer.data <- rbind(layer.data, data.temp)
        }
    } else {
        for(i in 1:6){
            file.name <- paste0(sample.name, "_L", i, "_barcodes.txt")
            file.path <- file.path(layer.input, file.name)

            data.temp <- read.table(file.path, header = FALSE, stringsAsFactors = FALSE) # assuming the file has no header
            data.temp <- data.frame(barcode = data.temp[,1], layer = paste0("layer", i), row.names = data.temp[,1])

            # Append to the final dataframe
            layer.data <- rbind(layer.data, data.temp)
        }
    }


    # For the WM file
    file.name <- paste0(sample.name, "_WM_barcodes.txt")
    file.path <- file.path(layer.input, file.name)

    data.temp <- read.table(file.path, header = FALSE, stringsAsFactors = FALSE) # assuming the file has no header
    data.temp <- data.frame(barcode = data.temp[,1], layer = "WM", row.names = data.temp[,1])

    # Append to the final dataframe
    layer.data <- rbind(layer.data, data.temp)

    sp_data <- AddMetaData(sp_data,
                        metadata = df_meta['V3'],
                        col.name = 'row')
    sp_data <- AddMetaData(sp_data,
                        metadata = df_meta['V4'],
                        col.name = 'col')
    sp_data <- AddMetaData(sp_data,
                        metadata = layer.data['layer'],
                        col.name = 'layer_guess_reordered')
    
    return(sp_data)
}


save_results <- function(sp_data, metrics_df, pca_data, dir.output) {
  # Save metrics
  write.csv(metrics_df, file = file.path(dir.output, 'metrics.csv'), row.names = FALSE)
          
  cluster_plot <- SpatialDimPlot(sp_data, label = TRUE, label.size = 3) + 
    labs(title=paste("PRECAST (ARI = ", round(metrics_df$ARI, 2), ")", sep="")) +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          legend.title = element_blank())
  
  ggsave(file.path(dir.output, 'clustering.pdf'), 
         plot = cluster_plot, width = 6, height = 6, dpi = 300)
  
  # Save reduced dimension data
  write.csv(pca_data,
            file = file.path(dir.output, 'low_dim_data.csv'),
            row.names = TRUE)
  
  # Save metadata
  write.csv(sp_data@meta.data,
            file = file.path(dir.output, 'cell_metadata.csv'),
            row.names = TRUE)
  
  # Save UMAP coordinates
  umap_coords <- as.data.frame(sp_data@reductions$UMAP@cell.embeddings)
  umap_coords$spot_id <- rownames(umap_coords)
  write.csv(umap_coords,
            file = file.path(dir.output, "spatial_umap_coords.csv"),
            row.names = FALSE)

  umap_plot <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = as.factor(sp_data@meta.data$cluster))) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "PRECAST", x = "UMAP 1", y = "UMAP 2", color = 'Cluster') +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.title = element_blank())
  ggsave(file.path(dir.output, 'umap.pdf'), plot = umap_plot,
         width = 6, height = 6, dpi = 300, device = "pdf")
}

run_sample <- function(input_path, sample.name, cluster.number) {
    dir.output <- file.path('/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/results3/DLPFC/PRECAST', sample.name)
    if (!dir.exists(dir.output)) {
        dir.create(dir.output, recursive = TRUE)
    }

    benchmark <- mark({
        sp_data <- load_dataset(input_path, sample.name, cluster.number)

        # Run PRECAST
        preobj <- CreatePRECASTObject(seuList = list(sp_data),
                          selectGenesMethod = "HVGs", 
                          gene.number = 2000,)
        preobj@seulist              
        PRECASTObj <- AddAdjList(preobj, platform = "Visium")
        PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 1, maxIter = 30, verbose = TRUE)

        PRECASTObj <- PRECAST(PRECASTObj, K = as.numeric(cluster.number))

        resList <- PRECASTObj@resList
        PRECASTObj <- SelectModel(PRECASTObj)

        seuInt <- PRECASTObj@seulist[[1]]
        seuInt@meta.data$cluster <- factor(unlist(PRECASTObj@resList$cluster))
        seuInt@meta.data$batch <- 1
        seuInt <- Add_embed(PRECASTObj@resList$hZ[[1]], seuInt, embed_name = "PRECAST")
        posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
        seuInt <- Add_embed(posList[[1]], seuInt, embed_name = "position")
        Idents(seuInt) <- factor(seuInt@meta.data$cluster)

        seuInt <- AddUMAP(seuInt, n_comp=2)

        gt <- seuInt@meta.data$layer_guess
        pred <- PRECASTObj@resList$cluster[[1]]
        pca_data <- PRECASTObj@resList$hZ[[1]]     
          
        metrics <- calculate_metrics(gt, pred, pca_data)

    }, iterations = 1L)

    execution_time <- as.numeric(benchmark$time[[1]]) 
    memory_usage <- as.numeric(benchmark$mem_alloc[[1]]) / (1024^2)  # Memory in MBsample.name <- '151673'

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
    save_results(seuInt, metrics_df, pca_data, dir.output)
}


input_path <- "/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/DLPFC12"
for (sample.name in names(batch_cluster_map)) {
    cluster.number <- batch_cluster_map[[sample.name]]
    run_sample(input_path, sample.name, cluster.number)
}

