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
  dplyr::filter(pct.1 > 0.9)

marker_genes_filt %>% dplyr::count(cluster)

se_sc_down <- SPOTlight::downsample_se_obj(
  se_obj = data,
  clust_vr = "nnet2",
  cluster_markers = marker_genes_filt,
  cl_n = 200,
  hvg = 0
)
data=data[]
spotlight_ls <- SPOTlight::spotlight_deconvolution(se_sc = data,
                                                   counts_spatial = test_spot_list[[1]],
                                                   clust_vr = "nnet2",
                                                   cluster_markers = marker_genes,
                                                   cl_n = 200,
                                                   hvg = 3000,
                                                   ntop = NULL,
                                                   transf = "uv",
                                                   method = "nsNMF",
                                                   min_cont = 0.09,
                                                   assay = "RNA",
                                                   slot = "counts")

spotlight_deconv <- spotlight_ls[[2]][, colnames(spotlight_ls[[2]]) != "res_ss"]
synthetic_comp <- as.matrix(synthetic_mixtures[[2]] /
                              rowSums(synthetic_mixtures[[2]]))
synthetic_comp=t(synthetic_comp)

spatial_decon_syn <-
  SPOTlight::test_synthetic_performance(
    test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols]),
    spot_composition_mtrx = spotlight_deconv[, ct_cols]
  )

test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols])
spot_composition_mtrx = spotlight_deconv[, ct_cols]
#test_spots_metadata_mtrx=as.matrix(test_spot_metadata[, ct_cols])
#spot_composition_mtrx=decon_mtrx[,ct_cols]


pearson_correlation_parameter=rep(0,1000)
for (i in 1:1000) {
  pearson_correlation_parameter[i]=cor(spot_composition_mtrx[i,],test_spots_metadata_mtrx[i,])
}
spotlight_deconv[which(spotlight_deconv<0.09)]=0
spotlight_deconv=spotlight_deconv/apply(spotlight_deconv,1,sum)
setwd("D:\\research\\undergraduate\\fda_ST\\SPOTlight")
saveRDS(pearson_SPOTlight,"pearson_SPOTlight.rds")
pearson_SPOTlight=pearson_correlation_parameter
saveRDS(spotlight_deconv,"spotlight_deconv.rds")


