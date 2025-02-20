library(SpatialPCA)
library(Seurat)
library(ggplot2)
library(Matrix)
library(SPARK)
library(Seurat)
library(peakRAM)


batch_cluster_map <- list(
  '151669' = 5, '151670' = 5, '151671' = 5, '151672' = 5,
  '151673' = 7, '151674' = 7, '151675' = 7, '151676' = 7,
  '151507' = 7, '151508' = 7, '151509' = 7, '151510' = 7
)

input_path <- "/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/DLPFC_spatialpca"
sample.name = '151676'
dir.output <- file.path('/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/Results/DLPFC/SpatialPCA', sample.name)
if (!dir.exists(dir.output)) {
  dir.create(dir.output, recursive = TRUE)
}

load(file.path(input_path, paste0('LIBD_sample', sample.name, ".RData")))

xy_coords = as.matrix(xy_coords)
rownames(xy_coords) = colnames(count_sub)

LIBD = CreateSpatialPCAObject(counts=count_sub, location=xy_coords, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions = 20, min.features=20)
mem <- peakRAM({
    start_time <- Sys.time()
    LIBD = SpatialPCA_buildKernel(LIBD, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
    LIBD = SpatialPCA_EstimateLoading(LIBD,fast=TRUE,SpatialPCnum=20)
    LIBD = SpatialPCA_SpatialPCs(LIBD, fast=TRUE)
    end_time <- Sys.time()
    T = end_time - start_time
})
T




# dir.input <- file.path('/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/DLPFC12', sample.name)
# dir.output <- file.path('/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/Results/DLPFC/SpatialPCA', sample.name, '/')
# meta.input <- file.path(dir.input, 'gt')
# layer.input <- file.path(dir.input, 'gt/layered')


# if(!dir.exists(file.path(dir.output))){
# dir.create(file.path(dir.output), recursive = TRUE)
# }

# filename <- paste0(sample.name, "_filtered_feature_bc_matrix.h5")
# sp_data <- Load10X_Spatial(dir.input, filename = filename)

# df_meta <- read.table(file.path(meta.input, 'tissue_positions_list_GTs.txt'))


# original_row_names <- row.names(df_meta)
# split_data <- strsplit(df_meta$V1, split = ",")
# df_meta <- do.call(rbind, lapply(split_data, function(x) {
# data.frame(V1=x[1], V2=x[2], V3=x[3], V4=x[4], V5=x[5], V6=x[6], V7=x[7])
# }))
# row.names(df_meta) <- df_meta$V1
# df_meta$V3 <- as.numeric(df_meta$V3)
# df_meta$V4 <- as.numeric(df_meta$V4)
# #df_meta_matched <- df_meta[df_meta$V1 %in% row.names(sp_data@meta.data),]
# # Set the row names of df_meta_matched to be V1
# # Identify the cells that are in both sp_data and df_meta
# common_cells <- colnames(sp_data[["Spatial"]]) %in% rownames(df_meta)

# # Subset sp_data to keep only these cells
# sp_data <- sp_data[, common_cells]

# # Initialize an empty dataframe to hold the final results
# layer.data <- data.frame()

# if(as.numeric(cluster.number) == 5) {
# for(i in 3:6){
#     file.name <- paste0(sample.name, "_L", i, "_barcodes.txt")
#     file.path <- file.path(layer.input, file.name)

#     data.temp <- read.table(file.path, header = FALSE, stringsAsFactors = FALSE) # assuming the file has no header
#     data.temp <- data.frame(barcode = data.temp[,1], layer = paste0("layer", i), row.names = data.temp[,1])

#     # Append to the final dataframe
#     layer.data <- rbind(layer.data, data.temp)
# }
# } else {
# for(i in 1:6){
#     file.name <- paste0(sample.name, "_L", i, "_barcodes.txt")
#     file.path <- file.path(layer.input, file.name)

#     data.temp <- read.table(file.path, header = FALSE, stringsAsFactors = FALSE) # assuming the file has no header
#     data.temp <- data.frame(barcode = data.temp[,1], layer = paste0("layer", i), row.names = data.temp[,1])

#     # Append to the final dataframe
#     layer.data <- rbind(layer.data, data.temp)
# }
# }


# # For the WM file
# file.name <- paste0(sample.name, "_WM_barcodes.txt")
# file.path <- file.path(layer.input, file.name)

# data.temp <- read.table(file.path, header = FALSE, stringsAsFactors = FALSE) # assuming the file has no header
# data.temp <- data.frame(barcode = data.temp[,1], layer = "WM", row.names = data.temp[,1])

# # Append to the final dataframe
# layer.data <- rbind(layer.data, data.temp)



# sp_data <- AddMetaData(sp_data,
#                     metadata = df_meta['V3'],
#                     col.name = 'row')
# sp_data <- AddMetaData(sp_data,
#                     metadata = df_meta['V4'],
#                     col.name = 'col')
# sp_data <- AddMetaData(sp_data,
#                     metadata = layer.data['layer'],
#                     col.name = 'layer_guess_reordered')
# count <- sp_data@assays$Spatial@layers$counts

# # get coordinates
# #gtlabels <- list(sp_data@meta.data$layer_guess_reordered)
# coord <- data.frame(row=sp_data@meta.data$row, col=sp_data@meta.data$col)
# row.names(coord) <- row.names(sp_data@meta.data)

# LIBD = CreateSpatialPCAObject(counts=count, location=as.matrix(coord), project = "SpatialPCA",gene.type="spatial",sparkversion="spark",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions = 20, min.features=20)

# LIBD = SpatialPCA_buildKernel(LIBD, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
# LIBD = SpatialPCA_EstimateLoading(LIBD,fast=FALSE,SpatialPCnum=20)
# LIBD = SpatialPCA_SpatialPCs(LIBD, fast=FALSE)

# clusterlabel <- walktrap_clustering(clusternum=as.numeric(cluster.number),latent_dat=LIBD@SpatialPCs,knearest=70 )