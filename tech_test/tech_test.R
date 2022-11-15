setwd("D:\\study\\undergraduate\\spolight\\code\\fda_ST\\fun_spotlight")

library(Seurat)
library(magrittr)
library(Biobase)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(SPOTlight)
library(FunSPOTlight)

setwd("D:\\study\\undergraduate\\spolight\\code\\fda_ST\\fun_spotlight")

tech_list = c(
  "Chromium",
  "inDrop",
  "C1HT-medium",
  "C1HT-small",
  "CEL-Seq2",
  "ddSEQ",
  "Drop-Seq",
  "ICELL8",
  "MARS-Seq",
  "Chromium (sn)",
  "Quartz-Seq2",
  "mcSCRB-Seq",
  "Smart-Seq2"
)
num_tech=length(tech_list)
num_repeat=1
n_syn_mixt=1000
test_tech=rep(0,length(tech_list))
syn_con=TRUE
marker_gene_access=0

data_performance_pearson = data.frame(matrix(rep(0,(3 * num_repeat * num_tech*n_syn_mixt)*3), 
                                             nrow = ( num_repeat * num_tech*n_syn_mixt),
                                             ncol = 3))
colnames(data_performance_pearson) = c("technology", "value", "num_repeat")
#data_performance_JSD[, 3] = rep(c("0.25", "0.50", "0.75"), num_cl_n * num_repeat)

data_performance_RMSE = data.frame(matrix(rep(0,(2 * num_repeat * num_tech*n_syn_mixt)*4),
                                          nrow = (2*num_repeat * num_tech*n_syn_mixt),
                                          ncol = 4))
colnames(data_performance_RMSE) = c("technology", "value", "num_repeat", "metric")

for (i in 1:num_tech) {
  data_performance_pearson[(1+(i-1)*num_repeat*n_syn_mixt):(i*num_repeat*n_syn_mixt),1]=tech_list[i]
  data_performance_RMSE[(1+(i-1)*num_repeat*n_syn_mixt*2):(i*num_repeat*n_syn_mixt*2),1]=tech_list[i]
  
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\fda_ST\\parameter_test\\data\\sce"
  )
  sc_obj <- readRDS(sprintf("sce%s.rds", i))
  
  if(sc_obj$batch[1]==tech_list[i]){
    test_tech[i]=1
  }
  
  colnames(sc_obj@assays@data$counts)[2]
  rownames(sc_obj@assays@data$counts)[1]
  dim(sc_obj@assays@data$counts)
  
  
  counts = sc_obj@assays@data$counts
  colname_counts=colnames(counts)
  #counts=cbind(rownames(counts),counts)
  #counts=data.frame(counts)
  #colnames(counts)=c("Genes", colname_counts)
  #counts=tibble::column_to_rownames(counts,"Genes")
  #rownames(counts)[duplicated(rownames(counts))]
  
  nnet2 = data.frame(sc_obj$nnet2)
  colnames(nnet2) = "nnet2"
  rownames(nnet2)=colname_counts
  
  sc_obj_seurat = CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    meta.data = nnet2
  )
  
  sc_obj_seurat <- sc_obj_seurat[, sc_obj_seurat$nnet2 != "Megakaryocytes"]
  sc_obj_seurat <- sc_obj_seurat[, sc_obj_seurat$nnet2 != "unclassified"]
  
  if(marker_gene_access==0){
    sc_obj_seurat <- Seurat::SCTransform(object = sc_obj_seurat,
                                 assay = "RNA",
                                 verbose = FALSE)
    Seurat::Idents(sc_obj_seurat) <- sc_obj_seurat$nnet2
    marker_genes <- Seurat::FindAllMarkers(
      object = sc_obj_seurat,
      assay = "SCT",
      slot = "data",
      min.pct = 0,
      only.pos = TRUE,
      logfc.threshold = 0
    )
    
    marker_genes %>% dplyr::count(cluster)
    
    marker_genes_filt <- marker_genes %>%
      dplyr::filter(pct.2 < 0.2)
    
    marker_genes_filt %>% dplyr::count(cluster)
    saveRDS(marker_genes,sprintf("marker_gene_%s",tech_list[i]))
  }else{
    marker_genes=readRDS(sprintf("marker_gene_%s",tech_list[i]))
    
    marker_genes_filt <- marker_genes %>%
      dplyr::filter(pct.2 < 0.2)
    
    marker_genes_filt %>% dplyr::count(cluster)
  }
  
  for (j in num_repeat) {
    
    if(syn_con){
      print(sprintf("Generating %s synthetic test mixtures", n_syn_mixt))
      set.seed(i*100)
      test_spot_ls <- test_spot_fun(se_obj = sc_obj_seurat,
                                    clust_vr = "nnet2",
                                    n = n_syn_mixt)
      setwd("D:\\research\\undergraduate\\fda_ST\\tech_test\\stdata")
      saveRDS(test_spot_ls,sprintf("test_spot_list_%s_%s",tech_list[i],1))
    }else{
      setwd("D:\\research\\undergraduate\\fda_ST\\tech_test\\stdata")
      test_spot_ls=readRDS(sprintf("test_spot_list_%s_%s",tech_list[i],1))
    }
    
    # Downsample scRNAseq to select gene set and number of cells to train the model
    sc_obj_seurat_down <- downsample_se_obj(se_obj = sc_obj_seurat,
                                    clust_vr = "nnet2",
                                    cluster_markers = marker_genes_filt,
                                    cl_n = 100,
                                    hvg = 0)
    
    
    print("Deconvolute synthetic spots")
    test_spot_counts=test_spot_ls[[1]]
    geno = as.matrix(sc_obj_seurat_down@assays$RNA@counts)
    ST.matrix = as.matrix(test_spot_counts)
    ST.matrix = ST.matrix[rownames(geno),]
    cell.type.factor = sc_obj_seurat_down$nnet2
    
    decon_mtrx = fun_spotlight2(geno = geno,
                               ST.matrix = ST.matrix,
                               cell.type.factor = cell.type.factor)
    
    decon_mtrx = t(decon_mtrx)
    
    fun_spotlight_deconv <-
      decon_mtrx[, colnames(decon_mtrx) != "ress_ss"]
    ct_cols <- colnames(fun_spotlight_deconv)
    test_spot_metadata=test_spot_ls[[2]]
    
    test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols])
    spot_composition_mtrx = fun_spotlight_deconv[, ct_cols]
    #test_spots_metadata_mtrx=as.matrix(test_spot_metadata[, ct_cols])
    #spot_composition_mtrx=decon_mtrx[,ct_cols]
    pearson_correlation_parameter=rep(0,1000)
    for (k in 1:1000) {
      pearson_correlation_parameter[k]=cor(spot_composition_mtrx[k,],test_spots_metadata_mtrx[k,])
    }
    data_performance_pearson[(1 + (i - 1) * 1000 * num_repeat+(j-1)*1000):((i - 1) * 1000 * num_repeat+j*1000),2]=pearson_correlation_parameter
    data_performance_pearson[(1 + (i - 1) * 1000 * num_repeat+(j-1)*1000):((i - 1) * 1000 * num_repeat+j*1000),3]=j
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
        root_square[k,2]=sqrt(root_absent/length(cell_absent))
      }
      
    }
    data_performance_RMSE[(1 + 2*(i - 1) * 1000 * num_repeat+(j-1)*2000):(2*(i - 1) * 1000 * num_repeat+j*2000-1000),2]=root_square[,1]
    data_performance_RMSE[(1 + 2*(i - 1) * 1000 * num_repeat+(j-1)*2000):(2*(i - 1) * 1000 * num_repeat+j*2000-1000),4]="present"
    data_performance_RMSE[(1001 + 2*(i - 1) * 1000 * num_repeat+(j-1)*2000):(2*(i - 1) * 1000 * num_repeat+j*2000),2]=root_square[,2]
    data_performance_RMSE[(1001 + 2*(i - 1) * 1000 * num_repeat+(j-1)*2000):(2*(i - 1) * 1000 * num_repeat+j*2000),4]="absent"
  }
  
}

setwd("D:\\research\\undergraduate\\fda_ST\\tech_test")
saveRDS(data_performance_pearson,"data_performance_pearson.rds")
saveRDS(data_performance_RMSE,"data_performance_RMSE")


for (i in 1:num_tech) {
  data_performance_pearson[(1+(i-1)*num_repeat*n_syn_mixt):(i*num_repeat*n_syn_mixt),1]=tech_list[i]
  data_performance_RMSE[(1+(i-1)*num_repeat*n_syn_mixt*2):(i*num_repeat*n_syn_mixt*2),1]=tech_list[i]
  for (j in 1:num_repeat) {
    
    setwd(sprintf("D:\\study\\undergraduate\\spolight\\code\\fda_ST\\parameter_test\\result\\data\\multip_tech%s",j))
    fun_spotlight_deconv=readRDS(sprintf("tech_%s.rds",tech_list[i]))[[1]]
    ct_cols <- colnames(fun_spotlight_deconv)
    
    test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols])
    spot_composition_mtrx = fun_spotlight_deconv[, ct_cols]
    #test_spots_metadata_mtrx=as.matrix(test_spot_metadata[, ct_cols])
    #spot_composition_mtrx=decon_mtrx[,ct_cols]
    pearson_correlation_parameter=rep(0,1000)
    for (k in 1:1000) {
      pearson_correlation_parameter[k]=cor(spot_composition_mtrx[k,],test_spots_metadata_mtrx[k,])
    }
    data_performance_pearson[(1 + (i - 1) * 1000 * num_repeat+(j-1)*1000):((i - 1) * 1000 * num_repeat+j*1000),2]=pearson_correlation_parameter
    
    
  }
}

library(ggplot2)
library(ggsci)
  
p1=ggplot(data_performance_pearson,aes(x=technology,y=value,fill=technology))+
  geom_boxplot()+
  theme_classic()+
  #ylim(0,1)+
  #scale_fill_locuszoom(alpha = 0.6)
  scale_fill_ucscgb(alpha = 0.4)

colnames(data_performance_RMSE)=c("Technology","RMSE","Num_repeat","Group")
plot_RMSE=ggplot2::ggplot(data = data_performance_RMSE,
                          aes(x=Technology,
                              y=RMSE,
                              fill=Group))+
  geom_boxplot()+
  theme_classic()+
  ylim(0,1)+
  #scale_fill_locuszoom(alpha = 0.6)
  scale_fill_npg(alpha = 0.6)+
  theme(
    #plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(angle = 60, hjust = 1),
    #axis.text.y = element_text(size = 18),
    #axis.ticks = element_blank(),
    #axis.text = element_blank(),
    #axis.title = element_text(size = 20),
    #axis.line = element_blank(),
    #legend.title = element_text(size=20),
    #legend.text = element_text(size=20)
  ) 
    
  guides(color=guide_legend(override.aes = list(size=7)))
