library(spacexr)
library(Matrix)

input_dir = "/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/lymph/"
load_seqFISH=function(n_genes){
    st_counts_fp=sprintf("%sOut_gene_expressions_%dgenes.csv",input_dir,n_genes)
    st_locations_fp=sprintf("%sOut_rect_locations.csv",input_dir)
    sc_counts_fp=sprintf("%sraw_somatosensory_sc_exp.txt",input_dir)
    sc_labels_fp=sprintf("%ssomatosensory_sc_labels.txt",input_dir)
    
    st_counts=read.csv(st_counts_fp,sep=",",row.names=1)
    st_counts=t(st_counts)
    st_locations=read.csv(st_locations_fp,sep=",",row.name=1)
    st_locations=st_locations[,c("x","y")]
    
    sc_counts=read.csv(sc_counts_fp,sep="\t",row.names=1)
    sc_labels=read.csv(sc_labels_fp,header=FALSE)$V1
    names(sc_labels)=colnames(sc_counts)
    
    ret_list=list(
            st_counts=st_counts,
            st_locations=st_locations,
            sc_counts=sc_counts,
            sc_labels=sc_labels
    )
    return(ret_list)
}

data = load_seqFISH(10000)
out_dir="/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/Results/Deconvolution/lymph/RCTD"
dir.create(out_dir,recursive = TRUE, showWarnings = FALSE)
out_matrix_norm_fp=file.path(out_dir,sprintf("RCTD_lymph_10000.csv"))

sc_reference=Reference(
    counts=data$sc_counts,
    cell_types=factor(data$sc_labels)
)

st_data=SpatialRNA(
    counts=data$st_counts,
    coords=data$st_location,
    require_int=FALSE
)

start_time <- Sys.time()

myRCTD <- create.RCTD(st_data, sc_reference, max_cores = 1, CELL_MIN_INSTANCE = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

end_time <- Sys.time()

weights=myRCTD@results$weights
norm_weights=normalize_weights(weights)
print(end_time-start_time)

write.csv(as.matrix(norm_weights),out_matrix_norm_fp)