library(Seurat)
setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set\\data1")
sc_data=read.csv("sc5rJUQ026.counts.csv")
setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set\\data1")
sc_metadata=read.csv("2103-Breastcancer_metadata.csv")

rownames(sc_data)=sc_data[,1]
sc_data=sc_data[,-1]


sc_metadata_sub=sc_metadata[which(sc_metadata$Cell%in%colnames(sc_data)),]
rownames(sc_metadata_sub)=sc_metadata_sub$Cell
sc_metadata_sub=sc_metadata_sub[colnames(sc_data),]
cell_type_factor = sc_metadata_sub$CellType
nnet2 = data.frame(cell_type_factor)
colnames(nnet2) = "nnet2"
rownames(nnet2) = sc_metadata_sub$Cell
data = CreateSeuratObject(counts = sc_data,
                          assay = "RNA",
                          meta.data = nnet2)

synthetic_mixtures <- SPOTlight::test_spot_fun(se_obj = data,
                                               clust_vr = "nnet2",
                                               n = 1000,
                                               verbose = TRUE)

setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set\\data2")
saveRDS(data,"data_sct.rds")
saveRDS(synthetic_mixtures,"test_spot_list.rds")
saveRDS(data,"data_sct.rds")

#setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set")
#marker_genes=readRDS("marker_genes_sc5rJUQ026")