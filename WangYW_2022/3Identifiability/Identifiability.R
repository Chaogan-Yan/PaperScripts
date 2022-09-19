# hierachical clustering
rm(list=ls())
require(R.matlab)
require(dplyr)
require(reshape2)
require(factoextra)

load_data <- function(datapath,methods,metric,ext){
  path <- datapath
  l <- list()
  for (i in 1:10){
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
CA4CP <-function(data){
  dist.r <- get_dist(data,method ="pearson")
  hc.r <- hclust(dist.r,method='ward.D2') # "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
  clusterGroups <- cutree(hc.r,k=41)
  acc.df <- data.frame(cbind(testdata$subid,clusterGroups))
  colnames(acc.df) <- c('subid','groupNum')
  acc.mat <- dcast(acc.df,subid~groupNum)
  acc <-length(which(acc.mat==3))/41
  return(acc)
}
CA4CE <-function(data){
  dist.r <- dist(data,method='euclidean')
  hc.r <- hclust(dist.r,method='ward.D2') # "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
  #res_coph <- cophenetic(hc.r)
  #cor(dist.r, res_coph)
  clusterGroups <- cutree(hc.r,k=41)
  acc.df <- data.frame(cbind(testdata$subid,clusterGroups))
  colnames(acc.df) <- c('subid','groupNum')
  acc.mat <- dcast(acc.df,subid~groupNum)
  acc <-length(which(acc.mat==3))/41
  return(acc)
}

#header
path<-'/Users/dianewang/Desktop/RfMRIHarmonization/scanner_harmonization/data_results/TSP3/infos/info/'
filename <- file.path(path,'sitesheader.csv')
testdata <- read.csv(filename)


#######################
#TSP3 - identifiability 
#######################
datapath <-'/Users/dianewang/Desktop/RfMRIHarmonization/scanner_harmonization/1Restart/targetsitechoice/TSP3'
methods <- c('raw','reg','adj','lmm','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','SMA_pku_ge','vae')
ext <- c('.mat','.csv','.csv','.csv','.mat','.mat','.mat','.mat','.mat','.csv')
metrics <- c('ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142')

cluster_accP <- c()
cluster_accE <- c()
for (m in metrics){
  l <- load_data(datapath,methods,m,ext)
  cluster_accP[[m]]<-lapply(l, CA4CP)
  cluster_accE[[m]]<-lapply(l, CA4CE)
}
clusteracc <- cbind(clusteracce,clusteraccp)


