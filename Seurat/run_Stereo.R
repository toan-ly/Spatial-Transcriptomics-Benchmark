# Import necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(clevr)
library(SingleCellExperiment)

# options(bitmapType = 'cairo')

# Define clustering parameters
n_domains <- 7

# Define input and output directories
data_path <- "Spatial-Transcriptomics-Benchmark/data/Mouse_Olfactory_Bulb/"
save_path <- "Spatial-Transcriptomics-Benchmark/results3/Mouse_Olfactory_Bulb/Seurat/"
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}


# Load raw.h5ad file
# cat("Loading raw.h5ad file...\n")
# raw_h5ad <- file.path(data_path)

# Convert .h5ad to Seurat-compatible .h5Seurat format
# h5seurat_file <- sub(".h5ad$", ".h5Seurat", raw_h5ad)
# if (!file.exists(h5seurat_file)) {
#   Convert(raw_h5ad, dest = "h5Seurat", overwrite = TRUE)
# }

# Load the converted file into Seurat
# sp_data <- LoadH5Seurat(h5seurat_file)

counts <- read.csv(file.path(data_path, "RNA_counts.tsv"), sep = "\t", row.names = 1)
position <- read.csv(file.path(data_path, "position.tsv"), sep = "\t")
used_barcodes <- read.csv(file.path(data_path, "used_barcodes.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(counts) <- gsub("^X", "Spot_", colnames(counts))
rownames(position) <- paste0('Spot_', position$label)
position <- position %>% select(x, y)

if (all(used_barcodes$V1 %in% colnames(counts))) {
  counts <- counts[, used_barcodes$V1]
} else {
  stop("Missing barcodes in counts")
}

counts <- counts[, used_barcodes$V1, drop = FALSE]
position <- position[used_barcodes$V1, ]
# counts[is.na(counts)] <- 0


sp_data <- CreateSeuratObject(counts = counts, assay = "Spatial")
sp_data <- AddMetaData(sp_data, metadata = position, col.name = c("x", "y"))

sp_data <- subset(sp_data, nCount_Spatial > 0)
# Data normalization using SCTransform
cat("Running SCTransform...\n")
sp_data <- SCTransform(sp_data, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction and clustering
cat("Performing dimensionality reduction and clustering...\n")
sp_data <- RunPCA(sp_data, assay = "SCT", verbose = FALSE, npcs = 50)
sp_data <- FindNeighbors(sp_data, reduction = "pca", dims = 1:30)

# Find optimal resolution for clustering
sp_data <- tryCatch({
  for (resolution in seq(1, 0.01, by = -0.01)) {
    sp_data <- FindClusters(sp_data, verbose = FALSE, resolution = resolution)
    if (length(unique(sp_data@meta.data$seurat_clusters)) == n_domains) {
      cat("Optimal resolution found: ", resolution, "\n")
      break
    }
  }
  sp_data
}, error = function(e) {
  stop("Error during clustering: ", e$message)
})

# Run UMAP for visualization
sp_data <- RunUMAP(sp_data, reduction = "pca", dims = 1:30)

# Visualization of clustering results
p1 <- DimPlot(sp_data, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(sp_data, label = TRUE, label.size = 3) + ggtitle("Seurat")
clusters_plot <- p1 + p2
ggsave(file.path(save_path, "clustering.pdf"), plot = p2, 
      width = 6, height = 6)

# Save results
cat("Saving results...\n")
write.table(sp_data@meta.data,
            file = file.path(save_path, "cell_metadata.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

