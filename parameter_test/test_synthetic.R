test_spot_ls <- test_spot_fun(se_obj = pdac_A,
                              clust_vr = "annotation",
                              n = 1000)
test_spot_counts <- as.matrix(test_spot_ls[[1]])
colnames(test_spot_counts) <- paste("mixt", 1:ncol(test_spot_counts), sep = "_")
test_spot_metadata <- test_spot_ls[[2]]
marker_genes=cluster_markers_a
marker_genes_filt <- marker_genes %>%
  dplyr::filter(pct.2 < 0.2)

marker_genes_filt %>% dplyr::count(cluster)



# Downsample scRNAseq to select gene set and number of cells to train the model
se_sc_down <- downsample_se_obj(se_obj = pdac_A,
                                clust_vr = "annotation",
                                cluster_markers = marker_genes_filt,
                                cl_n = 100,
                                hvg = 0)

geno = as.matrix(se_sc_down@assays$RNA@counts)
ST.matrix = test_spot_counts
ST.matrix = ST.matrix[rownames(geno),]
cell.type.factor = se_sc_down$annotation

decon_mtrx = fun_spotlight2(
  geno = geno,
  ST.matrix = ST.matrix,
  nfpca = 50,
  cell.type.factor = cell.type.factor,
  min_count = 0.01,
  Upper_Bound=10,
  Lower_Bound=0.1
)

decon_mtrx = t(decon_mtrx)

fun_spotlight_deconv <-
  decon_mtrx[, colnames(decon_mtrx) != "ress_ss"]
ct_cols <- colnames(fun_spotlight_deconv)
raw_statistics_ls <- test_synthetic_performance(test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols]),
                                                spot_composition_mtrx = fun_spotlight_deconv[, ct_cols])
test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols])
spot_composition_mtrx = fun_spotlight_deconv[, ct_cols]
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


