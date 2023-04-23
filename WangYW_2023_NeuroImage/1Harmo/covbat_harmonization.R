# CovBat- harmonize 
rm(list=ls())
library(R.matlab)
library(matrixStats)
require(reshape2)
require(factoextra)
source('~/R/covbat.R')
source('./R/utils.R')
source('./R/combat.R')

################################
#TSP3 - identifiability -covbat
################################
path<-'./TSP3/infos/info/'
filename <- file.path(path,'sitesheader.csv')
testdata <- read.csv(filename)
library(stats)
mod <- testdata[,c(2:3)] 
mod$gender <- as.factor(mod$gender)
mod$age <-as.numeric(mod$age)
mod<- matrix(data=unlist(mod),nrow=123,ncol=2,byrow = FALSE)

metrics <- c('ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142')
batch = as.factor(rep(c(1,2,3),each=41))
#setwd()
for (m in metrics){
  rawdata = readMat(paste0('./',m,'_raw.mat'))
  rawdata = rawdata$alldata
  covbat_harmonized <- covbat(t(rawdata),batch)
  writeMat(paste0('./',m,'_unadj_covbat.mat'),covbat = t(covbat_harmonized$dat.covbat))
  
  covbat_harmonized <- covbat(t(rawdata),batch,mod=mod)
  writeMat(paste0('./',m,'_adj_covbat.mat'),covbat = t(covbat_harmonized$dat.covbat))
  
}

################################
# CoRR - covbat
################################
corr_info <- readMat('./SubInfo_420.mat')
corr_info <-data.frame(matrix(unlist(corr_info),ncol=6,byrow=F))
colnames(corr_info) <- c('motion1','motion2','Subid','age','sex','site')
corr_batch <- corr_info$site
corr_mod <- corr_info[,c(4:5)]
corr_mod$sex <-as.numeric(corr_mod$sex)
corr_mod$sex[corr_mod$sex== 1]=2 #M
corr_mod$sex[corr_mod$sex==-1]=1 #F
corr_mod$sex <-as.factor(corr_mod$sex)
corr_mod$age <- as.numeric(corr_mod$age)
corrmod <-  matrix(data=unlist(corr_mod),nrow=420,ncol=2,byrow = FALSE)
sessdirs <- c('1','2')
for (ses in sessdirs){
for (m in metrics){
  rawdata = readMat(paste0('./',ses,'/',m,'_raw.mat'))
  rawdata = rawdata$raw
  covbat_harmonized <- covbat(t(rawdata),corr_batch)
  writeMat(paste0('./',ses,'/',m,'_unadj_covbat.mat'),covbat = t(covbat_harmonized$dat.covbat))
  
  covbat_harmonized <- covbat(t(rawdata),corr_batch,mod = corrmod)
  writeMat(paste0('./',ses,'/',m,'_adj_covbat.mat'),covbat = t(covbat_harmonized$dat.covbat))
}}

################################
# FCP - covbat with CoRR
################################
corrpath <- './'
FCP_info<-readMat('./600_subinfo.mat')
FCP_names <- names(FCP_info)
FCP_info <-data.frame(matrix(unlist(FCP_info),ncol=6,byrow=F))
colnames(FCP_info) <- FCP_names
FCP_batch <- FCP_info$SiteID
FCP_mod <- data.frame(as.numeric(FCP_info$Age),as.factor(FCP_info$Sex))
FCP_mod <- matrix(unlist(FCP_mod),nrow=600,ncol=2,byrow=FALSE)
mod <- rbind(corrmod,FCP_mod)
batch <- as.factor(c(corr_batch,FCP_info$SiteName))
sessdirs <- c('1','2')
for (ses in sessdirs){
  for (m in metrics){
    corr_rawdata = readMat(paste0(corrpath,ses,'/',m,'_raw.mat'))
    corr_rawdata = corr_rawdata$raw
    fcp_rawdata = readMat(paste0('./',m,'_raw.mat'))
    fcp_rawdata = fcp_rawdata$raw
    rawdata <- rbind(corr_rawdata,fcp_rawdata)
    
    covbat_harmonized <- covbat(t(rawdata),batch)
    writeMat(paste0('./',ses,'/',m,'_unadj_covbat.mat'),covbat = t(covbat_harmonized$dat.covbat))
    
    covbat_harmonized <- covbat(t(rawdata),batch,mod)
    writeMat(paste0('./',ses,'/',m,'_adj_covbat.mat'),covbat = t(covbat_harmonized$dat.covbat))
  }}

################################
# CovBat - SWU-beijing-cambridge
################################
sessdirs <- c('1','2')
for (ses in sessdirs){
  
    rawdata = readMat(paste0('./SWU_beijing_cambridge_',ses,'.mat'))
    data = rawdata$data
    batch = rawdata$batch
    mod = rawdata$mod
    covbat_harmonized <- covbat(t(rawdata),corr_batch)
    writeMat(paste0('./',ses,'/',m,'_unadj_covbat.mat'),covbat = t(covbat_harmonized$dat.covbat))
    
    covbat_harmonized <- covbat(t(data),batch,mod = mod)
    writeMat(paste0('./',ses,'/',m,'_adj_covbat.mat'),covbat = t(covbat_harmonized$dat.covbat))
  }
