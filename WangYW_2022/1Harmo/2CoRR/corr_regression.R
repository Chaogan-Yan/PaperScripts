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
library(data.table)
library(R.matlab)

lmm_byfeatures <- function(features, mod, formula, batch, batchvar,rf){ 
  # number of features
  V <- colnames(features)
  # total number of observations
  n <- nrow(features)
  l <-ncol(features)
  #one_hot_batch = as.matrix(data.frame(batch$scanner.A,batch$scanner.B,batch$scanner.C))
  b <- batch[1]
  x<-1
  batch_effect = data.frame(matrix(0,n,l))
  for (i in V) {
    feature_df <- cbind(features[i],mod,b)
    colnames(feature_df)<- c('feature',colnames(mod),batchvar)
    # unadjusted
    #
    form <- as.formula(paste0('feature', '~', formula,'+', rf))
    fit <- lmerTest::lmer(form, data=feature_df, REML=TRUE, control=lme4::lmerControl(optimizer='bobyqa'))
    
    # get the batch coef and remove it from data
    # coefbatch <- as.matrix(data.frame(ranef(fit))$condval)
    # batcheff <- one_hot_batch%*%coefbatch
    
    # batch_effect[x] = batcheff
    
    batch_effect[x] <- predict(fit,re.form=~0)+residuals(fit)
    x<-x+1
    print(i)
  }
  # reg_data <- features - batch_effect 
  return(batch_effect)
  #return(reg_data)
}
regressite_byfeatures <- function(features, mod, formula, batch, batchvar){ 
  
  # number of features
  V <- colnames(features)
  # total number of observations
  n <- nrow(features)
  l <-ncol(features)
  batch_effect = data.frame(matrix(0,n,l))
  # make coef basket
  one_hot_batch <- as.matrix(batch[2:length(batch)])
  n_level <- ncol(one_hot_batch)
  nm <- colnames(batch[2:length(batch)])
  coefbatch <- matrix(nrow=n_level,ncol=1,dimnames = list(nm))
  #prepare 
  b <- batch[1]
  x<-1
  for (i in V) {
    feature_df <- cbind(features[i],mod,b)
    colnames(feature_df)<- c('feature',colnames(mod),batchvar)
   
    form <- as.formula(paste0('feature', '~', formula,'+' , batchvar))    
    fit <- lm(form, data <- feature_df, x = TRUE)
    
    # get the batch coef and remove it from data
    coefbatch <- matrix(c(0,coef(fit)[nm[2:length(nm)]]))
    batcheff <- one_hot_batch%*%coefbatch
    
    batch_effect[x] = batcheff
    x<-x+1
  }
  
  reg_data <- features - batch_effect 
  return(reg_data)
}

## read head
head <- readMat('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/SubInfo/SubInfo_420.mat')
Motion <- data.frame(head["Motion"])
subid <- data.frame(unlist(head["SubID"]))
gender <-  data.frame(head["Sex"])
age <- data.frame(head["Age"])
site <- data.frame(head["Site"])

header <- cbind(subid,age,gender,Motion)
colnames(header) <- c('SubID','Age','Gender','Motion1','Motion2')


## do regresson
# prepare site dummyvar
for(unique_value in unique(site$Site)){
  site[paste("Site", unique_value, sep = "")] <- ifelse(site$Site == unique_value, 1, 0)
}
site$Site <- as.factor(site$Site)

# read data
dir_nm <- c('Results','S2_Results')
#index <- c('VMHC_FunImgARCWFsymS')
index <- c('ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142')
for (idir in dir_nm) {
  for (ndex in index){ 
  setwd(paste0('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',idir))
  raw <- readMat(paste0('./',ndex,'_raw.mat'))
  data <- data.frame(raw[['raw']])  
  if (idir == 'Results') 
    formu <- '1++Age+Gender+Motion1' else
    formu <- '1+Age+Gender+Motion2'   
  cat('begin to reck! /n')
  reg_corr <- regressite_byfeatures(features=data,
                               mod=header,
                               formula = 1,
                               batch= site,
                               batchvar = 'Site')
  writeMat(paste0('./',ndex,'_reg.mat'),reg=reg_corr)
  adj_corr <- regressite_byfeatures(features=data,
                               mod=header,
                               formula = formu,
                               batch= site,
                               batchvar = 'Site')
  writeMat(paste0('./',ndex,'_adj.mat'),adj=adj_corr)
  lmm_corr <- lmm_byfeatures(features=data,
                                         mod=header,
                                         formula = formu,
                                         batch= site,
                                         batchvar = 'Site',
                                          rf = '(1|Site)')
  writeMat(paste0('./',ndex,'_lmm.mat'),lmm=lmm_corr)
  #  write.csv(reg_corr,paste0('./',ndex,'_reg.csv'),row.names = FALSE)
  # write.csv(adj_corr,paste0('./',ndex,'_adj.csv'),row.names = FALSE)
  # write.csv(lmm_corr,paste0('./',ndex,'_lmm.csv'),row.names = FALSE)
}}

