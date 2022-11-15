library(ggplot2)
library(ggsci)


table(data_performance_pearson$cl_n)
#table(data_performance_pearson$value)
table(data_performance_pearson$num_repeat)
data_performance_pearson=data_performance_pearson[data_performance_pearson$cl_n!=0,]
num_repeat=5
num_cluster=length(table(data_performance_pearson$cl_n))
df_pearson=data.frame(matrix(rep(0,num_cluster*num_repeat*3),ncol = 3))
colnames(df_pearson)=c("cl_n","quantile","value")
k=1
for (i in names(table(data_performance_pearson$cl_n))) {
  df_pearson[(1+(k-1)*num_repeat*3):(k*num_repeat*3),1]=i
  for (j in 1:num_repeat) {
    df_pearson1=data_performance_pearson[data_performance_pearson$cl_n==i,]
    df_pearson1=df_pearson1[df_pearson1$num_repeat==j,]
    df_pearson[((k-1)*num_repeat*3+1+(j-1)*3):((k-1)*num_repeat*3+j*3),2]=c('0.25','0.5','0.75')
    df_pearson[((k-1)*num_repeat*3+1+(j-1)*3):((k-1)*num_repeat*3+j*3),3]=c(quantile(df_pearson1$value,0.25),
                                                                            quantile(df_pearson1$value,0.5),
                                                                            quantile(df_pearson1$value,0.75))
  }
  k=k+1
}

df_pearson$cl_n=factor(df_pearson$cl_n,levels = unique(df_pearson$cl_n))
p1=ggplot()+
  geom_point(data = df_pearson,
             aes(x=cl_n,y=value,color=quantile,shape=quantile))+
  ylim(0.65,1)+
  theme_classic()+
  geom_jitter(width = 5)+
  labs(x="knot number", y="Pearson's correlation")

for (i in 1:num_cl_n) {
  #cl_n = i * 10
  #data_performance_pearson[(1 + (i - 1) * 1000 * num_repeat):(1000*i*num_repeat), 1] = rep(cl_n, 1000*num_repeat)
  #data_performance_RMSE[(1 + 2*(i - 1) * 1000 * num_repeat):(2*1000*i*num_repeat),1]=rep(cl_n, 2000*num_repeat)
  
  for (j in 1:num_repeat) {
    data_performance_pearson[(1 + (i - 1) * 1000 * num_repeat+(j-1)*1000):((i - 1) * 1000 * num_repeat+j*1000),3]=j
    data_performance_RMSE[(1 + 2*(i - 1) * 1000 * num_repeat+(j-1)*2000):(2*(i - 1) * 1000 * num_repeat+j*2000),3]=j
  }
}
#num_cluster=2
num_repeat=5
table(data_performance_RMSE$cl_n)
data_performance_RMSE=data_performance_RMSE[1:100000,]
data_performance_RMSE=data_performance_RMSE[data_performance_RMSE$cl_n!=0,]
data_performance_RMSE[which(is.na(data_performance_RMSE$value)),]$value=0

df_RMSE=data.frame(matrix(rep(0,num_cluster*num_repeat*2*4),ncol = 4))
colnames(df_RMSE)=c("cl_n","quantile","value","cluster")
k=1
for (i in names(table(data_performance_pearson$cl_n))) {
  df_RMSE[(1+(k-1)*num_repeat*3*2):(k*num_repeat*3*2),1]=i
  for (j in 1:num_repeat) {
    df_RMSE1=data_performance_RMSE[data_performance_RMSE$cl_n==i,]
    df_RMSE1=df_RMSE1[df_RMSE1$num_repeat==j,]
    df_RMSE2=df_RMSE1[df_RMSE1$metric=='absent',]
    df_RMSE3=df_RMSE1[df_RMSE1$metric=='present',]
    df_RMSE[((k-1)*num_repeat*3*2+1+(j-1)*3*2):((k-1)*num_repeat*3*2+j*3*2),2]=c('0.25','0.5','0.75',
                                                                                 '0.25','0.5','0.75')
    df_RMSE[((k-1)*num_repeat*3*2+1+(j-1)*6):((k-1)*num_repeat*3*2+j*6-3),3]=c(quantile(df_RMSE2$value,0.25),
                                                                               quantile(df_RMSE2$value,0.5),
                                                                               quantile(df_RMSE2$value,0.75))
    df_RMSE[((k-1)*num_repeat*3*2+1+(j-1)*6):((k-1)*num_repeat*3*2+j*6-3),4]='absent'
    df_RMSE[((k-1)*num_repeat*3*2+(j-1)*6+4):((k-1)*num_repeat*3*2+j*6),3]=c(quantile(df_RMSE3$value,0.25),
                                                                             quantile(df_RMSE3$value,0.5),
                                                                             quantile(df_RMSE3$value,0.75)) 
    df_RMSE[((k-1)*num_repeat*3*2+(j-1)*6+4):((k-1)*num_repeat*3*2+j*6),4]='present'
  }
  k=k+1
}

#df_RMSE$cl_n=factor(df_RMSE$cl_n,levels = (c((1:20)*10)))
df_RMSE1=df_RMSE[df_RMSE$cluster=='present',]

library(tibble)
seg <- tibble(x = unique(df_RMSE$quantile)),
              y = c(0,0.05,0.10,0.15,0.2)

df_RMSE$cl_n=as.numeric(df_RMSE$cl_n)
p2=ggplot() +
  geom_tile(data = df_RMSE, 
            aes(x = cl_n, y =value, width =0, 
                height = 0,fill=quantile))+
  geom_point(data = df_RMSE[df_RMSE$cluster=="present",],
             aes(x=cl_n-3,y=value,shape=cluster,color=quantile))+
  geom_point(data = df_RMSE[df_RMSE$cluster=="absent",],
             aes(x=cl_n+3,y=value,shape=cluster,color=quantile))+
  labs(x="knot number", y="RMSE")+
  theme_classic()
df_RMSE$cl_n=as.numeric(df_RMSE$cl_n)

p2=ggplot() +
  geom_tile(data = df_RMSE, 
            aes(x = 0, y =0, width =0, 
                height = 0,fill=cluster))+
  geom_point(data = df_RMSE1,
             aes(x=cl_n,y=value,shape=cluster))+
  geom_point(data = filter(df_RMSE, cluster == "absent"),
             aes(x=cl_n+0.2,y=value,shape=quantile))+
  theme_classic()

df_RMSE1=df_RMSE[df_RMSE$cluster=='absent',]
p1=ggplot()+
  geom_point(data = df_RMSE1,
             aes(x=cl_n,y=value,color=quantile,shape=quantile))+
  theme_classic()+
  labs(x="knot number", y="RMSE")+
  geom_jitter(width = 5)

