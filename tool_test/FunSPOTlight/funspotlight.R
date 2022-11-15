library(magrittr)
library(FunSPOTlight)

setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set")
marker_genes=readRDS("marker_genes_sc5rJUQ026")

data <- Seurat::SCTransform(object = data,
                            assay = "RNA",
                            verbose = FALSE)
Seurat::Idents(data) <- data$nnet2
marker_genes <- Seurat::FindAllMarkers(
  object = data,
  assay = "RNA",
  slot = "data",
  min.pct = 0,
  only.pos = TRUE,
  logfc.threshold = 0
)

marker_genes %>% dplyr::count(cluster)

marker_genes_filt <- marker_genes %>%
  dplyr::filter(pct.1 > 0.1)
marker_genes_filt <- marker_genes %>%
  dplyr::filter(pct.2 < 0.02)
marker_genes_filt <- marker_genes %>%
  dplyr::filter(p_val_adj<0.01)

marker_genes_filt %>% dplyr::count(cluster)

se_sc_down <- SPOTlight::downsample_se_obj(
  se_obj = data,
  clust_vr = "nnet2",
  cluster_markers = marker_genes_filt,
  cl_n = 200,
  hvg = 0
)

geno = as.matrix(se_sc_down@assays[["RNA"]]@counts)
ST.matrix = as.matrix(test_spot_list[["topic_profiles"]])
ST.matrix = ST.matrix[rownames(geno),]
cell.type.factor = se_sc_down$nnet2

fun_decon=fun_spotlight2(geno,
                         ST.matrix,
                         cell.type.factor,
                         nkonts = 400,
                         nfpca = 25,
                         Upper_Bound = 100000,
                         Lower_Bound = 0.00001,
                         min_count = 0.09)

#fun_decon=fun_spotlight(geno,ST.matrix,cell.type.factor,min_count = 0.09,)
decon_mtrx = t(fun_decon)
fun_spotlight_deconv <-
  decon_mtrx[, colnames(decon_mtrx) != "ress_ss"]
test_spot_metadata=test_spot_list[[2]]
ct_cols <- colnames(fun_spotlight_deconv)
spatial_decon_syn <-
  SPOTlight::test_synthetic_performance(
    test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols]),
    spot_composition_mtrx = fun_spotlight_deconv[, ct_cols]
  )


test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols])
spot_composition_mtrx = fun_spotlight_deconv[, ct_cols]
#test_spots_metadata_mtrx=as.matrix(test_spot_metadata[, ct_cols])
spot_composition_mtrx=decon_mtrx[,ct_cols]
#spot_composition_mtrx=decon_mtrx[,ct_cols]

pearson_correlation_parameter=rep(0,1000)
for (i in 1:1000) {
  pearson_correlation_parameter[i]=cor(spot_composition_mtrx[i,],test_spots_metadata_mtrx[i,])
}
median(pearson_correlation_parameter)
quantile(pearson_correlation_parameter,0.25)
quantile(pearson_correlation_parameter,0.75)
mean(pearson_correlation_parameter)
var(pearson_correlation_parameter)
sd(pearson_correlation_parameter)

test_spots_metadata_mtrx[which(pearson_correlation_parameter<0.2),]

setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set\\data1\\plot_data")
saveRDS(pearson_FUNSPOTlight,"pearson_FunSPOTlight.rds")
saveRDS(fun_spotlight_deconv,"deconv_FunSPOTlight.rds")
saveRDS(root_square,"RMSE_FunSPOTlight.rds")
#test_spots_metadata_mtrx[which(pearson_correlation_parameter>0.8),]

#setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set\\data1\\rctd")
#saveRDS(pearson_RCTD,"pearson_RCTD.rds")
obs_test=test_spot_list[[2]]
obs_test=obs_test/apply(obs_test, 1, sum)
colnames(obs_test)
root_square=matrix(rep(0,2000),ncol=2)
for(i in 1:1000){
  cell_present=colnames(test_spot_list[[2]])[test_spot_list[[2]][i,]>0]
  root_present=0
  for (j in 1:length(cell_present)) {
    root_present=root_present+
      (decon_mtrx[i,cell_present[j]]-obs_test[i,cell_present[j]])^2
  }
  root_square[i,1]=sqrt(root_present/length(cell_present))
  cell_absent=colnames(test_spot_list[[2]])[test_spot_list[[2]][i,]==0]
  root_absent=0
  for (j in 1:length(cell_absent)) {
    root_absent=root_absent+
      (decon_mtrx[i,cell_absent[j]]-obs_test[i,cell_absent[j]])^2
  }
  root_square[i,2]=sqrt(root_absent/length(cell_absent))
  
}

quantile(root_square[,1],0.25)
quantile(root_square[,1],0.5)
quantile(root_square[,1],0.75)

decon_mtrx=rctd_deconv
decon_mtrx=cell2location_deconv
table(cell2location_deconv>0.09)
cell2location_deconv
apply(rctd_deconv, 1, sum)



