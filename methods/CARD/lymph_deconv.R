library(CARD)
library(Matrix)
library(matrixStats)


setwd('/home/lytq/Spatial-Transcriptomics-Benchmark/data/lymph')

org_st_count = read.csv('Out_gene_expressions_10000genes.csv',header = T, row.names = 1)
sc_exp = read.table('mapped_sc_exp.txt',header = T,row.names = 1)
sc_anno = read.table('somatosensory_sc_labels.txt',header = F, sep='\n')
spatial_location = read.csv('Out_rect_locations.csv',header = T, row.names = 1)

colnames(spatial_location) = c('x','y')

sc_count = Matrix(as.matrix(sc_exp), sparse = TRUE)
sc_meta = data.frame(matrix(ncol = 3, nrow = ncol(sc_exp)))
colnames(sc_meta) = c('cellID','cellType','sampleInfo')
sc_meta$sampleInfo = "sample1"
sc_meta$cellType = sc_anno$V1
sc_meta$cellID = colnames(sc_count)
rownames(sc_meta) = sc_meta$cellID

CARD_obj = createCARDObject(
  sc_count = as.matrix(sc_exp),
  sc_meta = sc_meta,
  spatial_count = t(as.matrix(org_st_count)),
  spatial_location = spatial_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 

CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

out_path = '/home/lytq/Spatial-Transcriptomics-Benchmark/Results/Deconvolution/lymph/CARD'
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
write.csv(CARD_obj@Proportion_CARD, file = paste0(out_path, '/CARD_lymph_10000.csv'), row.names = TRUE)

