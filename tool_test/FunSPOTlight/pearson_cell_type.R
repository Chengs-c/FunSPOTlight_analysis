#pearson_for_cell_type
n_cell_type=8
pearson_for_cell_type=matrix(rep(0,8),nrow = 8)
test_spot_metadata=test_spot_metadata/apply(test_spots_metadata_mtrx,1,sum)
test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols])
spot_composition_mtrx = funspot_pct002[, ct_cols]

rownames(pearson_for_cell_type)=colnames(test_spots_metadata_mtrx)
for (i in 1:8) {
  pearson_for_cell_type[i]=cor(spot_composition_mtrx[,i],test_spots_metadata_mtrx[,i])
}
