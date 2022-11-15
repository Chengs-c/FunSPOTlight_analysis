library(SPOTlight)
library(fda)
library(magrittr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
#setwd("D:\\study\\undergraduate\\spolight\\code\\fda_ST\\fun_spotlight")
#source("fun_spotlight.r")
library(FunSPOTlight)

#parameter
n_syn_mixt = 1000
data_name = "cerebellum_singlecell"
data_name = "se_Quartz"
marker_gene_access = 1
hvg = 0
num_fpca=10
num_repeat = 5
cl_n=40

#load data

if (data_name == "se_Quartz") {
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\SPOTlight_deconvolution_analysis-master\\analysis\\tool_benchmarking"
  )
  data = readRDS("se_quartz.rds")
  clust_vr = "nnet2"
  table(data$nnet2)
}

if (data_name == "cerebellum_singlecell") {
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\data_RCTD\\RCTD_data\\mouse cerebellum"
  )
  data = readRDS("1000cellsSubsampled_cerebellum_singlecell.rds")
  clust_vr = "liger_ident_coarse"
  table(data$liger_ident_coarse)
  ##cell_type = names(table(data$liger_ident_coarse))[which(table(data$liger_ident_coarse) >
  ##                                                          200)]
  #data=data[,data$liger_ident_coarse%in%cell_type]
  
  data_counts = data@assays$RNA@counts
  #data_counts <- data_counts[!duplicated(rownames(data_counts)), ]
  #data_counts=data_counts[sample(nrow(data_counts),round(0.5*nrow(data_counts))),]
  cell_type_factor = data$liger_ident_coarse
  nnet2 = data.frame(cell_type_factor)
  colnames(nnet2) = "nnet2"
  rownames(nnet2) = colnames(data_counts)
  
  data = CreateSeuratObject(counts = data_counts,
                            assay = "RNA",
                            meta.data = nnet2)
  #data = data[, data$nnet2 %in% cell_type]
  clust_vr = "nnet2"
}

#marker gene
if (marker_gene_access == 0) {
  se_sc <- Seurat::SCTransform(object = data,
                               assay = "RNA")
  # se_sc <- Seurat::SCTransform(object = data,
  #                              assay = "RNA",
  #                              verbose = FALSE)
  Seurat::Idents(se_sc) <- se_sc$nnet2
  marker_genes <- Seurat::FindAllMarkers(
    object = se_sc,
    assay = "SCT",
    slot = "data",
    min.pct = 0,
    only.pos = TRUE,
    logfc.threshold = 0
  )
  
  marker_genes %>% dplyr::count(cluster)
  
  marker_genes_filt <- marker_genes %>%
    dplyr::filter(pct.1 > 0.9)
  
  marker_genes_filt %>% dplyr::count(cluster)
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\fda_ST\\parameter_test\\cl_n\\markergene"
  )
  saveRDS(marker_genes, file = sprintf("%s_markergene.rds", data_name))
} else{
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\fda_ST\\parameter_test\\cl_n\\markergene"
  )
  marker_genes = readRDS(sprintf("%s_markergene.rds", data_name))
  
  marker_genes %>% dplyr::count(cluster)
  
  marker_genes_filt <- marker_genes %>%
    dplyr::filter(pct.1 > 0.9)
  
  marker_genes_filt = marker_genes %>%
    dplyr::filter(pct.2<0.1)
  
  marker_genes_filt %>% dplyr::count(cluster)
}

#test_spot=FALSE
#generate st data
for (i in 1:num_repeat) {
  print(sprintf("Generating %s synthetic test mixtures", n_syn_mixt))
  set.seed(i)
  
  test_spot_ls <- test_spot_fun(se_obj = data,
                                clust_vr = clust_vr,
                                n = n_syn_mixt)
  setwd(
    "D:\\research\\undergraduate\\fda_ST\\parameter_test\\synthetic_list"
  )
  saveRDS(test_spot_ls, file = sprintf("test_spot_ls%s.rds", i))
  
}



##
#data_performance = data.frame(matrix(0, nrow = (4 * num_fpca), ncol = 3))
#colnames(data_performance) = c("cl_n", "value", "metric")
#data_performance[, 3] = rep(c("accuracy", "sensitivity", "specificity", "F1"), num_fpca)

data_performance_pearson = data.frame(matrix(0, 
                                             nrow = (3 * num_repeat * num_fpca*n_syn_mixt),
                                             ncol = 3))
colnames(data_performance_pearson) = c("cl_n", "value", "num_repeat")
#data_performance_JSD[, 3] = rep(c("0.25", "0.50", "0.75"), num_fpca * num_repeat)

data_performance_RMSE = data.frame(matrix(0, nrow = (2 * num_repeat * num_fpca*n_syn_mixt),
                                          ncol = 4))
colnames(data_performance_RMSE) = c("cl_n", "value", "num_repeat", "metric")

for (i in 1:num_fpca) {
  nfpca = i * 10+10
  data_performance_pearson[(1 + (i - 1) * 1000 * num_repeat):(1000*i*num_repeat), 1] = rep(nfpca, 1000*num_repeat)
  data_performance_RMSE[(1 + 2*(i - 1) * 1000 * num_repeat):(2*1000*i*num_repeat),1]=rep(nfpca, 2000*num_repeat)
  
  for (j in 1:num_repeat) {
    data_performance_pearson[(1 + (i - 1) * 1000 * num_repeat+(j-1)*1000):((i - 1) * 1000 * num_repeat+j*1000),3]=j
    data_performance_RMSE[(1 + 2*(i - 1) * 1000 * num_repeat+(j-1)*2000):(2*(i - 1) * 1000 * num_repeat+j*2000),3]=j
    setwd(
      "D:\\research\\undergraduate\\fda_ST\\parameter_test\\synthetic_list"
    )
    test_spot_ls = readRDS(file = sprintf("test_spot_ls%s.rds", j))
    test_spot_counts <- as.matrix(test_spot_ls[[1]])
    colnames(test_spot_counts) <-
      paste("mixt", 1:ncol(test_spot_counts), sep = "_")
    test_spot_metadata <- test_spot_ls[[2]]
    #test_spot_metadata=test_spot_metadata/rowSums(test_spot_metadata)
    
    se_sc_down <- downsample_se_obj(
      se_obj = data,
      clust_vr = clust_vr,
      cluster_markers = marker_genes,
      cl_n = cl_n,
      hvg = hvg
    )
    
    print("Deconvolute synthetic spots")
    
    geno = as.matrix(se_sc_down@assays$RNA@counts)
    ST.matrix = test_spot_counts
    ST.matrix = ST.matrix[rownames(geno), ]
    cell.type.factor = se_sc_down$nnet2
    
    decon_mtrx = fun_spotlight2(
      geno = geno,
      ST.matrix = ST.matrix,
      cell.type.factor = cell.type.factor,
      nkonts = 200,
      nfpca = nfpca,
      ct_mode = "median",
      cluster_marker = marker_genes_filt
    )
    
    decon_mtrx = t(decon_mtrx)
    
    fun_spotlight_deconv <-
      decon_mtrx[, colnames(decon_mtrx) != "ress_ss"]
    ct_cols <- colnames(fun_spotlight_deconv)
    spatial_decon_syn <-
      test_synthetic_performance(
        test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols]),
        spot_composition_mtrx = fun_spotlight_deconv[, ct_cols]
      )
    test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols])
    spot_composition_mtrx = fun_spotlight_deconv[, ct_cols]
    #test_spots_metadata_mtrx=as.matrix(test_spot_metadata[, ct_cols])
    #spot_composition_mtrx=decon_mtrx[,ct_cols]
    pearson_correlation_parameter=rep(0,1000)
    for (k in 1:1000) {
      pearson_correlation_parameter[k]=cor(spot_composition_mtrx[k,],test_spots_metadata_mtrx[k,])
    }
    data_performance_pearson[(1 + (i - 1) * 1000 * num_repeat+(j-1)*1000):((i - 1) * 1000 * num_repeat+j*1000),2]=pearson_correlation_parameter
    
    obs_test=test_spot_ls[[2]]
    obs_test=obs_test/apply(obs_test, 1, sum)
    colnames(obs_test)
    root_square=matrix(rep(0,2000),ncol=2)
    for(k in 1:1000){
      cell_present=colnames(test_spot_ls[[2]])[test_spot_ls[[2]][k,]>0]
      root_present=0
      for (l in 1:length(cell_present)) {
        root_present=root_present+
          (decon_mtrx[k,cell_present[l]]-obs_test[k,cell_present[l]])^2
      }
      root_square[k,1]=sqrt(root_present/length(cell_present))
      cell_absent=colnames(test_spot_ls[[2]])[test_spot_ls[[2]][k,]==0]
      root_absent=0
      if(length(cell_absent!=0)){
        for (l in 1:length(cell_absent)) {
          root_absent=root_absent+
            (decon_mtrx[k,cell_absent[l]]-obs_test[k,cell_absent[l]])^2
        }
      }
      root_square[k,2]=sqrt(root_absent/length(cell_absent))
    }
    data_performance_RMSE[(1 + 2*(i - 1) * 1000 * num_repeat+(j-1)*2000):(2*(i - 1) * 1000 * num_repeat+j*2000-1000),2]=root_square[,1]
    data_performance_RMSE[(1 + 2*(i - 1) * 1000 * num_repeat+(j-1)*2000):(2*(i - 1) * 1000 * num_repeat+j*2000-1000),4]="present"
    data_performance_RMSE[(1001 + 2*(i - 1) * 1000 * num_repeat+(j-1)*2000):(2*(i - 1) * 1000 * num_repeat+j*2000),2]=root_square[,2]
    data_performance_RMSE[(1001 + 2*(i - 1) * 1000 * num_repeat+(j-1)*2000):(2*(i - 1) * 1000 * num_repeat+j*2000),4]="absent"
  }
  
}
setwd("D:\\research\\undergraduate\\fda_ST\\parameter_test\\nfpca")
saveRDS(data_performance_pearson,"data_performance_pearson.rds")
saveRDS(data_performance_RMSE,"data_performance_RMSE.rds")
saveRDS(data_performance_RMSE1,"data_performance_RMSE1.rds")
data_performance_RMSE1=data_performance_RMSE[data_performance_RMSE[,4]==0,]

cl_n_list=seq(1,20)*10
colourCount=20
#plot
performance_plt <- data_performance %>%
  dplyr::mutate(
    cl_n = factor(x = cl_n,
                  levels = cl_n_list)
  ) %>%
  ggplot() +
  geom_point(aes(x = cl_n,
                 y = value,
                 color = cl_n),
             size = 5,
             alpha = 0.9) +
  ylim(0.25,1)+
  facet_wrap(. ~ metric) +
  labs(title = "Deconvolution parameter benchmarking(cell number)",
       x = "cl_n/num_cell_per_celltype",
       y = "Metric Value") +
  theme_classic()  +
  scale_fill_manual(values = colorRampPalette(brewer.pal(20, "Set3"))(colourCount))+
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(size = 18, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 18),
    #axis.ticks = element_blank(),
    #axis.text = element_blank(),
    axis.title = element_text(size = 20),
    #axis.line = element_blank(),
    legend.title = element_text(size=20),
    legend.text = element_text(size=20)
  ) +
  guides(color=guide_legend(override.aes = list(size=7)))

performance_plt
colnames(data_performance_JSD)=c("cl_n","value","quantile")
quan_list=c("0.25","0.50","0.75")
#scale_color_brewer(palette = "Set3")
# Tool = factor(x = Tool,
#               levels = tech_list)

performance_JSD_plt <- data_performance_JSD %>%
  dplyr::mutate(
    quantile=factor(x=quantile,levels = quan_list)
  ) %>%
  ggplot() +
  geom_point(aes(x = cl_n,
                 y = value,
                 color = quantile),
             size = 1,
             alpha = 0.9) +
  ylim(0,0.5)+
  labs(title = "Deconvolution parameter benchmarking(cell number)",
       x = "cl_n",
       y = "JSD") +
  theme_classic()  +
  scale_fill_manual(values = colorRampPalette(brewer.pal(13, "Set3"))(colourCount))+
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(size = 18, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 18),
    #axis.ticks = element_blank(),
    #axis.text = element_blank(),
    axis.title = element_text(size = 20),
    #axis.line = element_blank(),
    legend.title = element_text(size=20),
    legend.text = element_text(size=20)
  ) +
  guides(color=guide_legend(override.aes = list(size=7)))


performance_JSD_plt
#facet_wrap(. ~ metric)

setwd("D:\\study\\undergraduate\\spolight\\code\\data_RCTD\\RCTD_data\\truedata\\plot_4.1\\data1\\cl_n")
data_performance=readRDS("data_performance.rds")
data_performance_JSD=readRDS("data_performance_JSD.rds")
saveRDS(data_performance,"data_performance.rds")
saveRDS(data_performance_JSD,"data_performance_JSD.rds")

