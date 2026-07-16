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
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(aricode)
  library(clevr)
  library(cluster)
  library(Matrix)
})

options(bitmapType = "cairo")
options(future.globals.maxSize = 8 * 1024^3)  # SCT on ~50k spots needs RAM

SEEDS <- c(42, 123, 456, 789, 2024)
SAMPLE_NAME <- "visium_hd_cancer_colon_square_016um"
N_CLUSTERS <- 6

DATA_PATH <- "/home/lytq/Spatial-Transcriptomics-Benchmark/data/colon_cancer"
OUTPUT_PATH <- "/home/lytq/Spatial-Transcriptomics-Benchmark/results/colon_cancer"

# --- include load_colon_visium_hd() and calculate_metrics() from above ---

search_resolution <- function(sp_data, n_clusters) {
  for (resolution in seq(1.0, 0.05, by = -0.02)) {
    sp_data <- FindClusters(sp_data, resolution = resolution, verbose = FALSE)
    if (length(unique(sp_data$seurat_clusters)) == n_clusters) {
      cat("Resolution:", resolution, "\n")
      return(sp_data)
    }
  }
  warning("Exact cluster count not found; using last resolution.")
  sp_data
}

plot_spatial_clusters <- function(meta, pred_col, title, out_file) {
  p <- ggplot(meta, aes(x = coord_x, y = -coord_y, color = .data[[pred_col]])) +
    geom_point(size = 0.05) +
    coord_equal() +
    labs(title = title, color = "Cluster") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(out_file, p, width = 6, height = 6, dpi = 300)
}

dir_input <- file.path(DATA_PATH, SAMPLE_NAME)
loaded <- load_colon_visium_hd(dir_input)

for (seed in SEEDS) {
  cat("\n==============================\nRUNNING SEED:", seed, "\n==============================\n")
  set.seed(seed)

  dir_out <- file.path(OUTPUT_PATH, as.character(seed), "Seurat")
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

  start_time <- Sys.time()
  gc(reset = TRUE)

  sp_data <- CreateSeuratObject(
    counts = loaded$counts,
    assay = "Spatial",
    project = SAMPLE_NAME
  )
  sp_data <- AddMetaData(sp_data, loaded$obs)
  sp_data$fine_annot_type <- loaded$ground_truth
  sp_data$coord_x <- loaded$coord_x
  sp_data$coord_y <- loaded$coord_y

  sp_data <- SCTransform(sp_data, assay = "Spatial", verbose = FALSE)
  sp_data <- RunPCA(sp_data, assay = "SCT", npcs = 50, verbose = FALSE)
  sp_data <- FindNeighbors(sp_data, reduction = "pca", dims = 1:30)
  sp_data <- search_resolution(sp_data, N_CLUSTERS)
  sp_data <- RunUMAP(sp_data, reduction = "pca", dims = 1:30)

  time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  mem_before <- gc()
  gc(reset = TRUE)
  mem_after <- gc()
  memory_used <- mem_before["Vcells", "max used"] - mem_after["Vcells", "max used"]

  metrics <- calculate_metrics(
    sp_data$fine_annot_type,
    sp_data$seurat_clusters,
    Embeddings(sp_data, "pca")
  )
  metrics$Time <- time_taken
  metrics$Memory <- memory_used
  write.csv(as.data.frame(metrics), file.path(dir_out, "metrics.csv"), row.names = FALSE)
  cat("ARI =", metrics$ARI, "\n")

  sp_data$pred_shift <- as.numeric(as.character(sp_data$seurat_clusters)) + 1

  meta <- sp_data@meta.data
  plot_spatial_clusters(
    meta, "pred_shift",
    sprintf("Seurat (ARI=%.4f)", metrics$ARI),
    file.path(dir_out, "clustering.pdf")
  )

  p_umap <- DimPlot(sp_data, reduction = "umap", group.by = "fine_annot_type") +
    ggtitle("Manual Annotation")
  p_umap2 <- DimPlot(sp_data, reduction = "umap", group.by = "pred_shift") +
    ggtitle("Seurat")
  ggsave(file.path(dir_out, "umap.pdf"), p_umap | p_umap2, width = 10, height = 4, dpi = 300)

  write.csv(Embeddings(sp_data, "pca"),
            file.path(dir_out, "low_dim_data.csv"), row.names = TRUE)
  write.csv(sp_data@meta.data,
            file.path(dir_out, "cell_metadata.csv"), row.names = TRUE)

  umap_df <- as.data.frame(Embeddings(sp_data, "umap"))
  umap_df$spot_id <- rownames(umap_df)
  write.csv(umap_df[, c("spot_id", "UMAP_1", "UMAP_2")],
            file.path(dir_out, "spatial_umap_coords.csv"), row.names = FALSE)

  rm(sp_data); gc()
}