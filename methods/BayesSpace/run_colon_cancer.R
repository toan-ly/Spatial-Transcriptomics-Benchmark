library(Matrix)

load_colon_visium_hd <- function(dir_input) {
  dir_input <- normalizePath(dir_input, winslash = "/")
  data_dir <- file.path(dir_input, "qc")
  if (!file.exists(file.path(data_dir, "counts.mtx"))) {
    data_dir <- dir_input
  }
  spatial_dir <- file.path(dir_input, "spatial")

  obs <- read.delim(file.path(data_dir, "observations.tsv"),
                    row.names = 1, check.names = FALSE)
  features <- read.delim(file.path(data_dir, "features.tsv"),
                         row.names = 1, check.names = FALSE)
  coords <- read.delim(file.path(data_dir, "coordinates.tsv"),
                       row.names = 1, check.names = FALSE)
  labels <- read.delim(file.path(dir_input, "labels.tsv"),
                       row.names = 1, check.names = FALSE)

  mat <- Matrix::readMM(file.path(data_dir, "counts.mtx"))
  mat <- as(mat, "CsparseMatrix")
  stopifnot(nrow(mat) == nrow(obs), nrow(mat) == nrow(coords))

  rownames(mat) <- rownames(obs)
  colnames(mat) <- rownames(features)

  selected <- tolower(as.character(obs$selected)) == "true"
  mat   <- mat[selected, , drop = FALSE]
  obs   <- obs[selected, , drop = FALSE]
  coords <- coords[rownames(obs), , drop = FALSE]

  ground_truth <- labels[rownames(obs), "label", drop = TRUE]
  keep <- !is.na(ground_truth) & ground_truth != "Outside"
  mat   <- mat[keep, , drop = FALSE]
  obs   <- obs[keep, , drop = FALSE]
  coords <- coords[rownames(obs), , drop = FALSE]
  ground_truth <- ground_truth[keep]

  coord_x <- coords[[1]]
  coord_y <- coords[[2]]

  cat("Loaded", nrow(mat), "spots x", ncol(mat), "genes\n")
  print(table(ground_truth))

  list(
    counts = t(mat),                 # genes x cells (Seurat / SCE convention)
    obs = obs,
    coord_x = coord_x,
    coord_y = coord_y,
    ground_truth = ground_truth,
    spatial_dir = spatial_dir
  )
}

calculate_metrics <- function(ground_truth, clusters, data_matrix) {
  valid <- !is.na(clusters) & !is.na(ground_truth)
  clusters <- clusters[valid]
  ground_truth <- ground_truth[valid]
  data_matrix <- data_matrix[valid, , drop = FALSE]

  list(
    ARI = adjustedRandIndex(ground_truth, clusters),
    AMI = aricode::AMI(ground_truth, clusters),
    Homogeneity = clevr::homogeneity(ground_truth, clusters),
    Completeness = clevr::completeness(ground_truth, clusters),
  # V_Measure = clevr::v_measure(ground_truth, clusters),
    ASW = mean(cluster::silhouette(
      as.numeric(factor(clusters)), dist(data_matrix)
    )[, 3]),
    CHAOS = NA,
    PAS = NA
  )
}

suppressPackageStartupMessages({
  library(BayesSpace)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scran)
  library(scater)
  library(ggplot2)
  library(aricode)
  library(clevr)
  library(cluster)
  library(Matrix)
})

SEEDS <- c(42, 123, 456, 789, 2024)
SAMPLE_NAME <- "visium_hd_cancer_colon_square_016um"
N_CLUSTERS <- 6
N_PCS <- 15
NREP <- 10000   # use 1000 for a quick test on ~50k spots

DATA_PATH <- "/home/lytq/Spatial-Transcriptomics-Benchmark/data/colon_cancer"
OUTPUT_PATH <- "/home/lytq/Spatial-Transcriptomics-Benchmark/results/colon_cancer"

# --- include load_colon_visium_hd() and calculate_metrics() from above ---

dir_input <- file.path(DATA_PATH, SAMPLE_NAME)
loaded <- load_colon_visium_hd(dir_input)

for (seed in SEEDS) {
  cat("\n==============================\nRUNNING SEED:", seed, "\n==============================\n")
  set.seed(seed)

  dir_out <- file.path(OUTPUT_PATH, as.character(seed), "BayesSpace")
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

  start_time <- Sys.time()
  gc(reset = TRUE)

  sce <- SingleCellExperiment(
    assays = list(counts = loaded$counts),
    colData = DataFrame(
      loaded$obs,
      fine_annot_type = loaded$ground_truth,
      array_row = as.integer(loaded$obs$row),
      array_col = as.integer(loaded$obs$col),
      pxl_row_in_fullres = loaded$coord_y,
      pxl_col_in_fullres = loaded$coord_x
    )
  )

  sce <- scuttle::logNormCounts(sce)

  set.seed(seed)
  dec <- scran::modelGeneVar(sce)
  top <- scran::getTopHVGs(dec, n = 2000)

  set.seed(seed)
  sce <- scater::runPCA(sce, subset_row = top, ncomponents = N_PCS)

  set.seed(seed)
  sce <- spatialPreprocess(sce, platform = "Visium", skip.PCA = TRUE)

  set.seed(seed)
  sce <- spatialCluster(
    sce, q = N_CLUSTERS, d = N_PCS, platform = "Visium",
    nrep = NREP, gamma = 3, save.chain = FALSE
  )

  set.seed(seed)
  sce <- runUMAP(sce, dimred = "PCA", name = "UMAP")

  time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  mem_before <- gc()
  gc(reset = TRUE)
  mem_after <- gc()
  memory_used <- mem_before["Vcells", "max used"] - mem_after["Vcells", "max used"]

  labels <- sce$spatial.cluster
  metrics <- calculate_metrics(
    sce$fine_annot_type,
    labels,
    reducedDim(sce, "PCA")
  )
  metrics$Time <- time_taken
  metrics$Memory <- memory_used
  write.csv(as.data.frame(metrics), file.path(dir_out, "metrics.csv"), row.names = FALSE)
  cat("ARI =", metrics$ARI, "\n")

  cluster_plot <- clusterPlot(sce, label = labels, size = 0.05) +
    labs(title = sprintf("BayesSpace (ARI=%.4f)", metrics$ARI)) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(dir_out, "clustering.pdf"), cluster_plot,
         width = 6, height = 6, dpi = 300)

  umap_coords <- as.data.frame(reducedDim(sce, "UMAP"))
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  umap_coords$spot_id <- rownames(umap_coords)

  umap_plot <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = factor(labels))) +
    geom_point(size = 0.1, alpha = 0.8) +
    labs(title = "BayesSpace", color = "Cluster") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(dir_out, "umap.pdf"), umap_plot, width = 10, height = 4, dpi = 300)

  write.csv(reducedDim(sce, "PCA"),
            file.path(dir_out, "low_dim_data.csv"), row.names = TRUE)
  write.csv(as.data.frame(colData(sce)),
            file.path(dir_out, "cell_metadata.csv"), row.names = TRUE)
  write.csv(umap_coords[, c("spot_id", "UMAP1", "UMAP2")],
            file.path(dir_out, "spatial_umap_coords.csv"), row.names = FALSE)

  rm(sce); gc()
}