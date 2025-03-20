library(BASS)
library(Matrix)
library(Seurat)
library(ggplot2)
library(bench)
library(clevr) # For homogeneity, completeness, v-measure
library(cluster) # For silhouette score
library(aricode)
library(mclust)
library(dplyr)
library(BayesSpace)
library(SeuratData)
library(pryr)

batch_cluster_map <- list("-0.04" = 8, "-0.09" = 8, "-0.14" = 8, "-0.24" = 8, "-0.29" = 8)

calculate_metrics <- function(ground_truth, clusters) {
  tryCatch(
    {
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

      return(list(
        ARI = ARI,
        AMI = AMI,
        Homogeneity = homogeneity,
        Completeness = completeness,
        V_Measure = v_measure
      ))
    },
    error = function(e) {
      warning("Error calculating metrics: ", e$message)
      return(list(
        ARI = NA,
        AMI = NA,
        Homogeneity = NA,
        Completeness = NA,
        V_Measure = NA
      ))
    }
  )
}

save_results <- function(BASS, metrics_df, pca_data, dir.output) {
  # Save metrics
  write.csv(metrics_df, file = file.path(dir.output, "metrics.csv"), row.names = FALSE)

  # cluster_plot <- SpatialDimPlot(sp_data, label = TRUE, label.size = 3, group.by = 'cluster') +
  #   labs(title=paste("BASS (ARI = ", round(metrics_df$ARI, 2), ")", sep="")) +
  #   theme(plot.title = element_text(hjust = 0.5, size = 16),
  #         legend.title = element_blank())

  # ggsave(file.path(dir.output, 'clustering.pdf'),
  #        plot = cluster_plot, width = 6, height = 6, dpi = 300)

  # Save reduced dimension data
  write.csv(pca_data,
    file = file.path(dir.output, "low_dim_data.csv"),
    row.names = TRUE
  )

  meta_df <- data.frame(
    BASS@xy,
    unlist(BASS@results$init_c), 
    unlist(BASS@results$init_z), 
    unlist(BASS@results$c), 
    unlist(BASS@results$z) 
  )
  colnames(meta_df) <- c("x", "y", "init_c", "init_z", "c", "pred")
  write.csv(meta_df,
    file = file.path(dir.output, "cell_metadata.csv"),
    row.names = TRUE
  )
}

load_dataset <- function(input_path, sample.name, cluster.number) {
  dir.input <- file.path(input_path, sample.name)

  filename <- paste0(input_path, "/MERFISH_Animal1_cnts.xlsx")
  cnts <- as.data.frame(read_excel(filename, sheet = sample.name))
  row.names(cnts) <- cnts[[1]]
  cnts <- cnts[-1]

  infoname <- paste0(input_path, "/MERFISH_Animal1_info.xlsx")
  info <- as.data.frame(read_excel(infoname, sheet = sample.name))
  row.names(info) <- info[[1]]
  info <- info[-1]

  C <- 20
  R <- as.numeric(cluster.number)

  # Run BASS
  BASS <- createBASSObject(
    list(cnts), list(info),
    C = C, R = R,
    beta_method = "SW"
  )

  return(BASS)
}

load_dataset2 <- function(input_path, sample.name, cluster.number) {
  load(sprintf("%s/MERFISH_Animal1.RData", input_path))
  # smps <- c("-0.04", "-0.09", "-0.14", "-0.19", "-0.24")
  cnts <- cnts_mult[sample.name] # a list of gene expression count matrices
  xys <- lapply(info_mult[sample.name], function(info.i){
    info.i$x <- info.i$x - min(info.i$x)
    info.i$y <- info.i$y - min(info.i$y)
    as.matrix(info.i[, c("x", "y")])
  }) # a list of spatial coordinates matrices
  # hyper-parameters
  C <- 20 # number of cell types
  R <- 8 # number of spatial domains

  BASS <- createBASSObject(cnts, xys, C = C, R = R, beta_method = "SW")

  return(BASS)
}

run_sample <- function(input_path, sample.name, cluster.number) {
  dir.output <- file.path("/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/Results/MERFISH/mHypothalamus/BASS", sample.name)
  if (!dir.exists(dir.output)) {
    dir.create(dir.output, recursive = TRUE)
  }
  dir.input <- file.path(input_path, sample.name)
  BASS <- load_dataset2(input_path, sample.name, cluster.number)

  start_time <- Sys.time()
  start_mem <- pryr::mem_used()
  # benchmark <- mark(
  #   {
      BASS <- BASS.preprocess(
        BASS,
        doLogNormalize = TRUE,
        doPCA = TRUE,
        scaleFeature = TRUE,
        nPC = 20
      )

      BASS <- BASS.run(BASS)
      # plot(1:BASS@burnin, BASS@samples$beta, xlab = "Iteration", 
      #   ylab = expression(beta), type = "l")

      BASS <- BASS.postprocess(BASS)

      pred <- BASS@results$z[[1]]
      gt <- BASS@results$c[[1]]

      metrics <- calculate_metrics(gt, pred)
  #   },
  #   iterations = 1L
  # )

  # execution_time <- as.numeric(benchmark$time[[1]])
  # memory_usage <- as.numeric(benchmark$mem_alloc[[1]]) / (1024^2) # Memory in MBsample.name <- '151673'
  end_time <- Sys.time()
  execution_time <- as.numeric(end_time - start_time)

  end_mem <- pryr::mem_used()
  memory_usage <- as.numeric(end_mem - start_mem) / (1024^2)  # Memory in MB

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
  save_results(BASS, metrics_df, pca_data, dir.output)
}


input_path <- "/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/mHypothalamus"
for (sample.name in names(batch_cluster_map)) {
  cat("Processing batch:", sample.name, "\n")
  cluster.number <- batch_cluster_map[[sample.name]]
  run_sample(input_path, sample.name, cluster.number)
}
