setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE\\data1")
write.table(test_spot_list[[1]],file = "st_count.txt",sep = "\t")

setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE\\data1\\result")
stride_deconv=read.table(file = "result_spot_celltype_frac.txt",sep = "\t")
colnames(stride_deconv)=stride_deconv[1,]
stride_deconv=stride_deconv[-1,]
rownames(stride_deconv)=stride_deconv[,1]
stride_deconv=stride_deconv[,-1]
test_spot_metadata=test_spot_list[[2]]

colnames(stride_deconv) <-
  gsub(
    pattern = "[[:punct:]]|[[:blank:]]",
    ".",
    x = colnames(stride_deconv),
    perl = TRUE
  )
stride_deconv=as.matrix(stride_deconv)
stride_deconv=apply(stride_deconv, 2, as.numeric)
ct_cols <- colnames(stride_deconv)
#stride_deconv=t(stride_deconv)
spatial_decon_syn <-
  test_synthetic_performance(
    test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols]),
    spot_composition_mtrx = spot_composition_mtrx[, ct_cols]
  )
test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols])
stride_deconv=STRIDE_deconv
spot_composition_mtrx = stride_deconv[, ct_cols]
spot_composition_mtrx[which(spot_composition_mtrx<0.09)]=0
spot_composition_mtrx=spot_composition_mtrx/apply(spot_composition_mtrx, 1, sum)
#test_spots_metadata_mtrx = as.matrix(synthe_com1[, ct_cols])
#spot_composition_mtrx = synthe_res_norm[, ct_cols]

#test_spots_metadata_mtrx=as.matrix(test_spot_metadata[, ct_cols])
#spot_composition_mtrx=decon_mtrx[,ct_cols]


pearson_correlation_parameter=rep(0,1000)
for (i in 1:1000) {
  pearson_correlation_parameter[i]=cor(spot_composition_mtrx[i,],test_spots_metadata_mtrx[i,])
}
pearson_STRIDE=pearson_correlation_parameter
setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set\\data1\\STRIDE")
saveRDS(pearson_STRIDE,"pearson_STRIDE.rds")
saveRDS(stride_deconv,file="STRIDE_deconv.rds")










saveRDS()

