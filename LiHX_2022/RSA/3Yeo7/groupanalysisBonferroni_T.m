clear;clc;
path='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/convergentTwoEMeanBold/1Yeo7';

maskid=dir([path,'/Yeo*']);
for imask=1:size(maskid,1)
    imaskpath=fullfile(path,maskid(imask).name);
    
    E1=load([imaskpath,'/E1',maskid(imask).name,'_SubResultsrp.mat']);
    E2=load([imaskpath,'/E2',maskid(imask).name,'_SubResultsrp.mat']);
    
    %one sample-t test 
    convergent2E=[E1.SubResults(:,3);E2.SubResults(:,3)];
    [h,p,ci,stats]=ttest(convergent2E);
    Ttestresult{imask,1}=maskid(imask).name;
    Ttestresult{imask,2}=stats.tstat;
    Ttestresult{imask,3}=p;
     
    %%%partial spearman correlation 
    convergent2E_partial=[E1.SubResults(:,6);E2.SubResults(:,6)];
    [h_partial,p_partial,ci_partial,stats_partial]=ttest(convergent2E_partial);
    Ttestresult_partial{imask,1}=maskid(imask).name;
    Ttestresult_partial{imask,2}=stats_partial.tstat;
    Ttestresult_partial{imask,3}=p_partial; 
end
save([path,'/yeo7result.mat'],'Ttestresult')
save([path,'/partial_yeo7result.mat'],'Ttestresult_partial')  

