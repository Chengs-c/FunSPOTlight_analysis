



library(spacexr)
library(Seurat)
library(dplyr)

## Parameters
set.seed(321)
num_sample_cell_type = 2
num_syn = 5


setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set\\data2")
se_quartz <-readRDS(file = "data.rds")
synthetic_mixtures <-readRDS(file = "test_spot_list.rds")

cell_type = as.factor(se_quartz$nnet2)
cell_type=as.factor(gsub("/",".",cell_type))
names(cell_type)=colnames(se_quartz@assays[["RNA"]]@counts)
count = as.matrix(se_quartz@assays$RNA@counts)
reference = Reference(counts = count, cell_types = cell_type)
    
setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set\\data2\\RCTD")
   tmp <- synthetic_mixtures[[1]] %>%
      data.frame() %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::mutate(spot_0 = spot_1) %>%
      dplyr::select(gene, spot_0, everything()) %>%
      readr::write_csv(x = ,
                       file = "MappedDGEForR.csv",
                       col_names = TRUE)
    
    # Create aritificial coordinates for the synthetic spots
    #nrow(synthetic_mixtures[[2]])
    
    # Since we have 1000 spots we can create an array of 20 * 50 matrix
    coord <- expand.grid(1:20, 1:50)
    colnames(coord) <- c("xcoord", "ycoord")
    df_coord <-
      data.frame("barcodes" = paste("spot", 1:1000, sep = "_"), coord)
    
    setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set\\data2\\RCTD")
    readr::write_csv(x = df_coord,
                     file = "BeadLocationsForR.csv")
    
    
    ## RCTD deconvolution
    ### Read data in RCTD
    # setwd("D:/study/undergraduate/spolight/code/SPOTlight_deconvolution_analysis-master")
    # reference <- spacexr::Reference(refdir = "analysis/tool_benchmarking/RCTD_data/reference")
    # #reference <- RCTD::dgeToSeurat(refdir = "analysis/tool_benchmarking/RCTD_data/reference")
    
    #setwd("D:/study/undergraduate/spolight/code/SPOTlight_deconvolution_analysis-master")
    #puck <- spacexr::read.SpatialRNA(datadir = "analysis/tool_benchmarking/RCTD_data/spatial")
    
    setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set\\data2")
    puck = spacexr::read.SpatialRNA(datadir = "RCTD",
                                    count_file = "MappedDGEForR.csv",
                                    coords_file = "BeadLocationsForR.csv")
    
    
    ### Creating and running RCTD
    
    # myRCTD <- RCTD::create.RCTD(spatialRNA = puck,
    #                             reference = reference,
    #                             max_cores = 1,
    #                             CELL_MIN = 18)
    myRCTD <- spacexr::create.RCTD(
      spatialRNA = puck,
      reference = reference,
      max_cores = 1,
      CELL_MIN = 0
    )
    
    
    
    myRCTD <- spacexr::run.RCTD(RCTD = myRCTD,
                                doublet_mode = "multi")
    #?spacexr::run.RCTD()
    
    
    ### Assess results
    results <- myRCTD@results
    
    
    
    #norm_weights <- sweep(results$weights, 1, rowSums(results$weights), '/')
    #rctd_deconv <- myRCTD@results$weights
    rctd_deconv1 = matrix(0, nrow = length(results), ncol = length(unique(cell_type)))
    
    colnames(rctd_deconv1) = unique(cell_type)
    rctd_deconv = rctd_deconv1
    for (k in 1:nrow(rctd_deconv1)) {
      rctd_deconv1[k, names(results[[k]]$sub_weights)] = results[[k]]$sub_weights
      if (sum(rctd_deconv1[k, ]) == 0) {
        rctd_deconv1[k, names(results[[k]]$all_weights)] = results[[k]]$all_weights
        
        for (z in 1:ncol(rctd_deconv1)) {
          rctd_deconv[k, z] = rctd_deconv1[k, z] / sum(rctd_deconv1[k, ])
        }
        for (z in 1:ncol(rctd_deconv1)) {
          if (rctd_deconv[k, z] < 0.09) {
            rctd_deconv1[k, z] = 0
          }
        }
        for (z in 1:ncol(rctd_deconv1)) {
          rctd_deconv[k, z] = rctd_deconv1[k, z] / sum(rctd_deconv1[k, ])
        }
      } else{
        rctd_deconv[k, ] = rctd_deconv1[k, ]
      }
    }
    
    
    
    ### Save results
    
    setwd(
      "D:\\study\\undergraduate\\spolight\\code\\fda_ST\\test_common_synthetic\\test2\\result\\decon\\spacexr"
    )
    saveRDS(object = rctd_deconv,
            file = sprintf("rctd_deconv_%s_%s.RDS", i, j))
    



 colnames(rctd_deconv) <-
   gsub(
     pattern = "[[:punct:]]|[[:blank:]]",
     ".",
     x = colnames(rctd_deconv),
     perl = TRUE
   )
# 
# spot_perform=SPOTlight::test_synthetic_performance(
#   test_spots_metadata_mtrx = rctd_deconv[, colnames(synthetic_comp)],
#   spot_composition_mtrx = synthetic_comp)

synthetic_comp=as.matrix(synthetic_mixtures[[2]])
colnames(synthetic_comp) <-
  gsub(
    pattern = "[[:punct:]]|[[:blank:]]",
    ".",
    x = colnames(synthetic_comp),
    perl = TRUE
  )
colnames(rctd_deconv) <-
  gsub(
    pattern = "[[:punct:]]|[[:blank:]]",
    ".",
    x = colnames(rctd_deconv),
    perl = TRUE
  )


spot_perform=SPOTlight::test_synthetic_performance(
  test_spots_metadata_mtrx = rctd_deconv[, colnames(synthetic_comp)],
  spot_composition_mtrx = synthetic_comp)

setwd("D:\\research\\undergraduate\\fda_ST\\STRIDE_data_set\\data2\\RCTD")
saveRDS(myRCTD,"myRCTD.rds")
saveRDS(rctd_deconv,"rctd_deconv.rds")