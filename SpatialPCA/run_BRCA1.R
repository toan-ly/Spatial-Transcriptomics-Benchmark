library(SpatialPCA)
library(ggplot2)
library(Matrix)
library(Seurat)

input_path <- "/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/BRCA1/V1_Human_Breast_Cancer_Block_A_Section_1"
adata <- Load10X_Spatial(input_path)
df_meta <- read.table(file.path(input_path, 'metadata.tsv'), sep = '\t', header=TRUE)
adata <- AddMetaData(adata, metadata = df_meta$fine_annot_type, col.name = 'fine_annot_type')
adata <- adata[, df_meta[, 1]]

gene_counts <- rowSums(GetAssayData(adata, slot = "counts") > 0)
keep_genes <- names(gene_counts[gene_counts >= 10])
adata <- subset(adata, features = keep_genes)


df_coords <- read.table(file.path(input_path, 'spatial/tissue_positions_list.csv'), sep = '\t', header=TRUE)
coords_rownames <- row.names(df_coords)
split_data <- strsplit(df_coords$V1, split = ",")


xy_coords <- adata@images$slice1@coordinates
xy_coords <- xy_coords[c('imagerow', 'imagecol')]
colnames(xy_coords) <- c('x_coord', 'y_coord')