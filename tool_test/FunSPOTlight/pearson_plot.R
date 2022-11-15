library(ggplot2)
name_list=c("pearson_FUNSPOTlight",
            "pearson_cell2location",
            "pearson_DSTG",
            "pearson_RCTD",
            "pearson_STRIDE",
            "pearson_SPOTlight")
pearson_FUNSPOTlight=as.matrix(pearson_FUNSPOTlight)
pearson_cell2location=as.matrix(pearson_cell2location)
pearson_DSTG=as.matrix(pearson_DSTG)
pearson_RCTD=as.matrix(pearson_RCTD)   
pearson_STRIDE=as.matrix(pearson_STRIDE)
pearson_SPOTlight=as.matrix(pearson_SPOTlight)
pearson_Data_Frame=rbind(pearson_FUNSPOTlight,pearson_cell2location)
pearson_Data_Frame=rbind(pearson_Data_Frame,pearson_DSTG)
pearson_Data_Frame=rbind(pearson_Data_Frame,pearson_RCTD)
pearson_Data_Frame=rbind(pearson_Data_Frame,pearson_STRIDE)
pearson_Data_Frame=rbind(pearson_Data_Frame,pearson_SPOTlight)

dim(pearson_Data_Frame)
name_matrix=matrix(c(rep(name_list[1],1000),
              rep(name_list[2],1000),
              rep(name_list[3],1000),
              rep(name_list[4],1000),
              rep(name_list[5],1000),
              rep(name_list[6],1000)),
              nrow=6000)
pearson_DataFrame=cbind(pearson_Data_Frame,name_matrix)
pearson_DataFrame=as.data.frame(pearson_DataFrame)
colnames(pearson_DataFrame)=c("pearson","method")

p=ggplot(data = pearson_DataFrame, aes(x=method,y=pearson,group=method,fill=method))+
  stat_boxplot()+
  geom_boxplot()

p1=ggplot(pearson_DataFrame)+geom_boxplot(aes(x=method,y=pearson))

pearson_benchmark <- pearson_DataFrame %>%
  dplyr::mutate(
    pearson = as.numeric(pearson),
    method = factor(x = method,
                  levels = name_list)) %>%
  ggplot(aes(x=method,y=pearson,group=method,fill=method)) +
  stat_boxplot(geom="errorbar",wideth=0.15)+
  geom_boxplot()+
  # geom_jitter(alpha = 0.5) +
  #geom_violin(alpha = 0.8) +
  #geom_boxplot()+
  theme_classic() +
  #scale_fill_brewer(palette = "Set2") +
  #scale_color_brewer(palette = "Set2") +
  #scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set3"))(8))+
  theme(
    #plot.title = element_text(hjust = 0.5, size = 30),
    #strip.text = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 5),
    axis.title = element_text(size = 10),
    #legend.title = element_text(size=10),
    #legend.text = element_text(size=5)
    ) +
  ylim(0,1)
  guides(color=guide_legend(override.aes = list(size=5)))

