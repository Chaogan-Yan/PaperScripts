# hierachical clustering
rm(list=ls())
require(R.matlab)
require(dplyr)
require(reshape2)
require(factoextra)

###################
##### funcs #######
###################
load_data <- function(datapath,methods,metric,ext){
  path <- datapath
  l <- list()
  for (i in 1:length(methods)){
    d <- paste0(path,'/',metric,'_',methods[i],ext[i])
    if (ext[i]=='.mat'){
      data <- readMat(d)
      if (metric=='FC_D142'){
        data <- data.frame(matrix(unlist(data),ncol=10011,byrow=F))
      }else{
        data <- data.frame(matrix(unlist(data),ncol=38810,byrow=F))}
    }else{
      data <-read.csv(d)
      if (nrow(data)!=123){
        data <-read.csv(d,header = FALSE)
      }
    }
    l <- append(l,list(data))
  }
  return(l)
}
#header
path<-'./info/'
filename <- file.path(path,'sitesheader.csv')
testdata <- read.csv(filename) # including "subid" which shall be used in the function below

CA4C <-function(data,distance.method,clustering.method){ ## you can add param-"treenum" if you have different number of traveling subjects
  dist.r <- get_dist(data,method = distance.method)
  hc.r <- hclust(dist.r,method=clustering.method) # "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
  clusterGroups <- cutree(hc.r,k=41) #treenum k= num of traveling-subjects
  acc.df <- data.frame(cbind(testdata$subid,clusterGroups))

  colnames(acc.df) <- c('subid','groupNum')
  acc.mat <- dcast(acc.df,subid~groupNum)
  acc <-length(which(acc.mat==3))/41
  return(acc)
}
#######################
#TSP3 - identifiability 
#######################

datapath <-'/Users/dianewang/Desktop/RfMRIHarmonization/scanner_harmonization/data_results/TSP3/processing/ResultsS'
methods <- c('raw','reg','adj','lmm','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','adj_covbat','unadj_covbat','SMA','ICVAE')
ext <- c('.mat','.mat','.mat','.mat','.mat','.mat','.mat','.mat','.mat','.mat','.mat','.mat')
metrics <- c('ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142')

cluster_accS <- c()
cluster_accC <- c()
cluster_accA <- c()
cluster_accW <- c()
cluster_accW2 <- c()

cluster_accSE <- c()
cluster_accCE <- c()
cluster_accAE <- c()
cluster_accWE <- c()
cluster_accWE2 <- c()
for (m in metrics){
  l <- load_data(datapath,methods,m,ext)
  cluster_accA[[m]]<-lapply(l, CA4C,distance.method='pearson', clustering.method='average')
  cluster_accC[[m]]<-lapply(l, CA4C,distance.method='pearson',clustering.method='complete')
  cluster_accS[[m]]<-lapply(l, CA4C,distance.method='pearson',clustering.method='single')
  cluster_accW[[m]]<-lapply(l, CA4C,distance.method='pearson',clustering.method='ward.D')
  cluster_accW2[[m]]<-lapply(l, CA4C,distance.method='pearson',clustering.method='ward.D2')
  
  cluster_accAE[[m]]<-lapply(l, CA4C,distance.method='euclidean', clustering.method='average')
  cluster_accCE[[m]]<-lapply(l, CA4C,distance.method='euclidean',clustering.method='complete')
  cluster_accSE[[m]]<-lapply(l, CA4C,distance.method='euclidean',clustering.method='single')
  cluster_accWE[[m]]<-lapply(l, CA4C,distance.method='euclidean',clustering.method='ward.D')
  cluster_accWE2[[m]]<-lapply(l, CA4C,distance.method='euclidean',clustering.method='ward.D2')
}

acc<-list(cluster_accA,cluster_accAE,cluster_accC,cluster_accCE,cluster_accS,cluster_accSE,cluster_accW,cluster_accWE,cluster_accW2,cluster_accWE2)
dis_link <- c('Pearson_average','Pearson_complete','Pearson_single','Pearson_wardD','Pearson_wardD2','Euclidean_average','Euclidean_complete','Euclidean_single','Euclidean_wardD','Euclidean_wardD2')
acc_long <- melt(acc)
colnames(acc_long) <- c('Accuracy','Method','Metric','DistanceLinkage')

acc_long$Method <- factor(acc_long$Method,levels=c('raw','reg','adj','lmm','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','adj_covbat','unadj_covbat','SMA','ICVAE'),labels=methods)
acc_long$Metric <- as.factor(acc_long$Metric)
acc_long$DistanceLinkage<- factor(acc_long$DistanceLinkage,levels=c(1:10),labels=dis_link)
acc_long_new <- acc_long %>%
  mutate(Metric_dis_link = paste0(DistanceLinkage,'_',Metric))%>%
  select(Accuracy,Method,Metric_dis_link,Metric)

acc_wide <- reshape(acc_long_new,
                    idvar = "Method",
                    timevar = "Metric_dis_link",
                    direction = "wide")
write.csv(acc_wide,file='/Users/dianewang/Desktop/RfMRIHarmonization/revision_analysis/Results/linkage-influence(1)/acc.csv')

############
## here insert the  "exchange target site experiment"
############
datapath <-'./TSP3-harmonized'
methods <- c('raw','reg','adj','lmm','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','adj_covbat','unadj_covbat','SMA','ICVAE','SMA_ipcas_ge')
ext <- c('.mat','.mat','.mat','.mat','.mat','.mat','.mat','.mat','.mat','.mat','.mat','.mat','.mat')

cluster_accS <- c()
cluster_accC <- c()
cluster_accA <- c()
cluster_accW <- c()
cluster_accW2 <- c()

cluster_accSE <- c()
cluster_accCE <- c()
cluster_accAE <- c()
cluster_accWE <- c()
cluster_accWE2 <- c()


for (m in metrics){
  l <- load_data(datapath,methods,m,ext)
  cluster_accA[[m]]<-lapply(l, CA4C,distance.method='pearson', clustering.method='average')
  cluster_accC[[m]]<-lapply(l, CA4C,distance.method='pearson',clustering.method='complete')
  cluster_accS[[m]]<-lapply(l, CA4C,distance.method='pearson',clustering.method='single')
  cluster_accW[[m]]<-lapply(l, CA4C,distance.method='pearson',clustering.method='ward.D')
  cluster_accW2[[m]]<-lapply(l, CA4C,distance.method='pearson',clustering.method='ward.D2')
  
  cluster_accAE[[m]]<-lapply(l, CA4C,distance.method='euclidean', clustering.method='average')
  cluster_accCE[[m]]<-lapply(l, CA4C,distance.method='euclidean',clustering.method='complete')
  cluster_accSE[[m]]<-lapply(l, CA4C,distance.method='euclidean',clustering.method='single')
  cluster_accWE[[m]]<-lapply(l, CA4C,distance.method='euclidean',clustering.method='ward.D')
  cluster_accWE2[[m]]<-lapply(l, CA4C,distance.method='euclidean',clustering.method='ward.D2')
  
  ## organize  to the wide format
  acc<-list(cluster_accA,cluster_accC,cluster_accS,cluster_accW,cluster_accW2,cluster_accAE,cluster_accCE,cluster_accSE,cluster_accWE,cluster_accWE2)
  dis_link <- c('Pearson_average','Pearson_complete','Pearson_single','Pearson_wardD','Pearson_wardD2','Euclidean_average','Euclidean_complete','Euclidean_single','Euclidean_wardD','Euclidean_wardD2')
  acc_long <- melt(acc)
  colnames(acc_long) <- c('Accuracy','Method','Metric','DistanceLinkage')
  
  acc_long$Method <- methods
  acc_long$Metric <- as.factor(acc_long$Metric)
  acc_long$DistanceLinkage<- factor(acc_long$DistanceLinkage,levels=c(1:10),labels=dis_link)
  acc_long_new <- acc_long %>%
    mutate(Metric_dis_link = paste0(DistanceLinkage,'_',Metric))%>%
    select(Accuracy,Method,Metric_dis_link) 
  
  acc_long_new$Accuracy <- scales::percent(round(as.numeric(acc_long_new$Accuracy),2))
  acc_wide <- reshape(acc_long_new,
                      idvar = "Method",
                      timevar = "Metric_dis_link",
                      direction = "wide")
}

################
#### Rank ######
################
data <- as.matrix(acc_wide[,c(2:51)])
rank.acc <- matrix(nrow=12,ncol=50)
for (i in c(1:50)){
  rank.acc[,i]<- rank(data[,i])
}  
rank.acc <- cbind(methods,data.frame(rank.acc))
colnames(rank.acc) = colnames(acc_wide)

library(reshape2)
acc_rank_long <- melt(rank.acc)
acc_rank_long_new <-data.frame(cbind(acc_rank_long$value,acc_rank_long$Method,as.character(acc_long$Metric),as.character(acc_long$DistanceLinkage)))
colnames(acc_rank_long_new) <- c('AccuracyRank','Method','Metric','DistanceLinkage')

acc_rank_long_new$AccuracyRank <- as.numeric(acc_rank_long_new$AccuracyRank)
acc_rank_long_new$Method<-factor(acc_rank_long_new$Method,levels=c('raw','reg','adj','lmm','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','adj_covbat','unadj_covbat','SMA','ICVAE'),labels=methods)
acc_rank_long_new$Metric<-as.factor(acc_rank_long_new$Metric)
acc_rank_long_new$DistanceLinkage<-as.factor(acc_rank_long_new$DistanceLinkage)


###################################
#### for each metric, we apply permutation test between combat-series/covbat and SMA
###################################
library(coin)
library(dplyr)
combat.series <- c('para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','adj_covbat','unadj_covbat')
pval = matrix(nrow=length(metrics),ncol=length(combat.series))
 
rank.pval = matrix(nrow=length(metrics),ncol=length(combat.series)) 
r=1
for (m in metrics){ 
  c=1
  for (cb in combat.series){
    
    test.data <- acc_rank_long_new %>%
      filter(Method=='SMA' | Method== cb) %>%
      filter(Metric==m)
    
    res <- oneway_test(AccuracyRank~Method|Metric,data=test.data,
                       distribution=approximate(nresample = 50000),
                       ties.method = "average-scores",  #same results by "median-ranks"
                       alternative="larger",
                       conf.level=0.99)
    temp <- data.frame(pvalue(res))
    
    rank.pval[r,c] <- round(temp$pvalue.res.,2)
    c=c+1
  }
  r=r+1
}


