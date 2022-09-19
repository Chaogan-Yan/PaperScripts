## regression function
## arguments
  # data ——  file path where the regressed data, n sub* m features,
  #           each column is a feature
  # mod -- unadjusted or adjusted by age, gender and so on
  # formula -- the linear regression model,feature ~ 1+mod+batch
  # batch -- batch list,e.p.[1,1,1,2,2,2], where 1 refer reference batch
  # @return the regressed data 
## wyw 2020/12/16
rm(list=ls())
library(lme4) 
library(R.matlab)

args = commandArgs(trailingOnly=TRUE)
InputName <- args[1]
feature_name <- args[2]
feature_num <- args[3]
OutputPath <- args[4]


lmm_byfeatures <- function(features, mod, formula, batch, batchvar){ 
  # number of features
  V <- colnames(features)
  # total number of observations
  n <- nrow(features)
  l <- ncol(features)
  #one_hot_batch = as.matrix(data.frame(batch$scanner.A,batch$scanner.B,batch$scanner.C))
  b <- batch[1]
  x<-1
  batch_effect = data.frame(matrix(0,n,l))
  for (i in V) {
    feature_df <- cbind(features[i],mod,b)
    colnames(feature_df)<- c('feature',colnames(mod),'scanner')
    # unadjusted
    #
    form <- as.formula(paste0('feature', '~', formula,'+',batchvar))    
    fit <- lme4::lmer(form, data=feature_df, REML=TRUE, control=lme4::lmerControl(optimizer='bobyqa'))
    
    # get the batch coef and remove it from data
    # coefbatch <- as.matrix(data.frame(ranef(fit))$condval)
    # batcheff <- one_hot_batch%*%coefbatch
    
    # batch_effect[x] = batcheff
    
    batch_effect[x] <- predict(fit,re.form=~0)+residuals(fit)
    x<-x+1
  }
  # reg_data <- features - batch_effect 
  return(batch_effect)
  #return(reg_data)
}
regressite_byfeatures <- function(features, mod, formula, batch, batchvar){ 
  #debugonce(regressite_byfeatures)
  # number of features
  V <- colnames(features)
  # total number of observations
  n <- nrow(features)
  l <-ncol(features)
  one_hot_batch = as.matrix(batch[2:4])
  b <- batch[1]
  x<-1
  batch_effect = data.frame(matrix(0,n,l))
  for (i in V) {
    feature_df <- cbind(features[i],mod,b)
    colnames(feature_df)<- c('feature',colnames(mod),batchvar)
    # unadjusted
    #
    form <- as.formula(paste0('feature', '~', formula,'+' , batchvar))    
    fit <- lm(form, data <- feature_df, x = TRUE)
    
    # get the batch coef and remove it from data
    coefbatch <- matrix(c(0,coef(fit)["scannerPKU_GE"],coef(fit)["scannerPKU_SIEMENS"]))
    batcheff <- one_hot_batch%*%coefbatch
    
    batch_effect[x] = batcheff
    x<-x+1
    }
  
  reg_data <- features - batch_effect 
  return(reg_data)
}
library(openxlsx)
path<-'/mnt/Data3/RfMRILab/Wangyw/harmonization_project/TST/TSTafterDpabi/info'
filename <- file.path(path,'subinfo.xlsx')
testdata <- read.xlsx(filename)
header <- as.data.frame(testdata[,c(1,4,6)])
colnames(header)<-c('subid','gender','age')
mod <- header[c(1:41),c(2,3)]
mod$gender <- as.factor(mod$gender)
scanner <- data.frame(rep(c('IPCAS_GE','PKU_GE','PKU_SIEMENS'),c(41,41,41)))
colnames(scanner)<-c('scanner')
for(unique_value in unique(scanner$scanner)){
  scanner[unique_value] <- ifelse(scanner$scanner == unique_value, 1, 0)
}
OutputPath <- '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/TSP3/ResultsS'
# feature_name <- 'FC_D142'
#  raw_FC <-readMat()
#  data<- data.frame(matrix(unlist(raw_FC),ncol=10011,byrow=F))
names = c('ALFF_FunImgARCW','fALFF_FunImgARCW','ReHo_FunImgARCWF','DegreeCentrality_FunImgARCWF','FC_D142')
# #data <- read.csv('/Users/dianewang/Desktop/project/scanner_harmonization/model_code/stability/test.csv',header = FALSE)
# setwd('/Users/dianewang/Desktop/project/scanner_harmonization/data_results/3SIte41Sub/temp')
 require(doParallel)
 cl1 <- makeCluster(5)
registerDoParallel(cl1)
foreach (name = names,.packages=c("R.matlab","lmerTest")) %dopar% {
  InputName <- paste0(OutputPath,'/',name,'_raw.mat')
  raw <- readMat(InputName)
  data <- raw[['alldata']]

  lmm_reg <- lmm_byfeatures(features <- data,
                            mod <-mod,
                            formula <- '1+age+gender',
                            batch <- scanner,
                            batchvar <- '(1|scanner)')
  writeMat(paste0(OutputPath,'/',name,'_lmm.mat'),lmm=lmm_reg)
  #rm(lmm_reg)
  adj <- regressite_byfeatures(features=data,
                             mod=mod,
                             formula = '1+gender+age',
                             batch = scanner,
                             batchvar = 'scanner')
  writeMat(paste0(OutputPath,'/',name,'_adj.mat'),adj=adj)
  #rm(adj)
  reg <- regressite_byfeatures(features=data,
                             mod=mod,
                             formula = '1',
                             batch= scanner,
                             batchvar = 'scanner')
  writeMat(paste0(OutputPath,'/',name,'_reg.mat'),reg=reg)
  #rm(reg)
}
clusterExport(cl, 'l')