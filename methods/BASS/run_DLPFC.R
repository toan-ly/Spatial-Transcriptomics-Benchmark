library(BASS)
library(Matrix)
library(Seurat)
library(ggplot2)
library(bench)
library(clevr)  # For homogeneity, completeness, v-measure
library(cluster) # For silhouette score
library(aricode)
library(mclust)
library(dplyr)
library(BayesSpace)
library(SeuratData)


batch_cluster_map <- list(
#   '151669' = 5, '151670' = 5, '151671' = 5, '151672' = 5,
#   '151673' = 7, '151674' = 7, '151675' = 7, '151676' = 7,
  '151507' = 7, '151508' = 7, '151509' = 7, '151510' = 7
)

calculate_metrics <- function(ground_truth, clusters) {
  tryCatch({
    # Remove NA values
    valid_indices <- !is.na(clusters) & !is.na(ground_truth)
    clusters <- clusters[valid_indices]
    ground_truth <- ground_truth[valid_indices]

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

    return(list(ARI = ARI, 
                AMI = AMI, 
                Homogeneity = homogeneity, 
                Completeness = completeness, 
                V_Measure = v_measure))
  }, error = function(e) {
    warning("Error calculating metrics: ", e$message)
    return(list(ARI = NA, 
                AMI = NA, 
                Homogeneity = NA, 
                Completeness = NA, 
                V_Measure = NA))
  })
}

save_results <- function(sp_data, metrics_df, pca_data, dir.output) {
  # Save metrics
  write.csv(metrics_df, file = file.path(dir.output, 'metrics.csv'), row.names = FALSE)
          
  cluster_plot <- SpatialDimPlot(sp_data, label = TRUE, label.size = 3, group.by = 'cluster') + 
    labs(title=paste("BASS (ARI = ", round(metrics_df$ARI, 2), ")", sep="")) +
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

}

run_sample <- function(input_path, sample.name, cluster.number) {
    dir.output <- file.path('/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/results3/DLPFC/BASS', sample.name)
    if (!dir.exists(dir.output)) {
        dir.create(dir.output, recursive = TRUE)
    }
    dir.input <- file.path(input_path, sample.name)
    sp_data <- Load10X_Spatial(dir.input, filename = 'filtered_feature_bc_matrix.h5')
    df_meta <- read.table(file.path(dir.input, 'metadata.tsv'))
    sp_data <- AddMetaData(sp_data, metadata = df_meta$layer_guess, col.name = 'layer_guess')
    
    smp <- sample.name
    if (smp %in% c('151669', '151670', '151671', '151672')) {
        index <- 2
    } else if (smp %in% c('151673', '151674', '151675', '151676')) {
        index <- 3
    } else if (smp %in% c('151507', '151508', '151509', '151510')) {
        index <- 1
    }
    benchmark <- mark({
        load(sprintf('/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/spatialLIBD_p%d.RData', index))
        C = 20
        R = as.numeric(cluster.number)

        # Run BASS
        BASS <- createBASSObject(cntm[smp], xym[smp], C = C, R = R, 
            beta_method = "SW", init_method = "mclust",
            nsample = 1000)

        BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE,
            geneSelect = "sparkx", nSE = 3000, doPCA = TRUE,
            scaleFeature = FALSE, nPC = 20)

        BASS <- BASS.run(BASS)
        BASS <- BASS.postprocess(BASS)

        pred <- BASS@results$z[[1]]
        gt <- sp_data@meta.data$layer_guess
          
        metrics <- calculate_metrics(gt, pred)
        sp_data <- AddMetaData(sp_data, metadata = pred, col.name = 'cluster')

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
      ASW = NA,
      Time = execution_time,
      Memory = memory_usage
    )

    pca_data <- as.data.frame(t(BASS@X_run))
    # change index to row names
    rownames(pca_data) <- colnames(sp_data)
    save_results(sp_data, metrics_df, pca_data, dir.output)
}


input_path <- "/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/DLPFC_new"
for (sample.name in names(batch_cluster_map)) {
    cat("Processing batch:", sample.name, "\n")
    cluster.number <- batch_cluster_map[[sample.name]]
    run_sample(input_path, sample.name, cluster.number)
}
