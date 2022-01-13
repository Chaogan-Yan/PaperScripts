clear;clc;
path='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/convergentTwoexperiment/RSAROIanalysis/schaefer400PatrialSpearman/1paper86partialSpearman/';
outpath=[path,'statisticanalysis/'];
mkdir(outpath)
E1=load([path,'subresults/E1PartialSpearman_schaefer400_SubResultsrp.mat']);
E2=load([path,'subresults/E2PartialSpearman_schaefer400_SubResultsrp.mat']);

%one sample-t test
convergent2E=[E1.ROI_MultipleLabel(:,:,3);E2.ROI_MultipleLabel(:,:,3)];
[h,p,ci,stats]=ttest(convergent2E);
significantnum=length(find(h==1));
nonsignificantnum=length(find(h==0));

% correct: Bonferroni method    0.05/400
correctsig_Bonferroninum=length(find(p<0.05/400));
ROI_index=find(p<0.05/400);
save([outpath,'ttestresultBonferroni_New.mat'],'h','p','ci','stats','significantnum','nonsignificantnum','correctsig_Bonferroninum','ROI_index')

Schaeferpath='/mnt/Data3/RfMRILab/Lihuixian/DPABI_V5.1_201230/Templates/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii';
[MaskData,MaskVox,MaskHead]=y_ReadRPI(Schaeferpath);
MaskROI=reshape(MaskData,1,[]);
Element = unique(MaskROI);
Element(find(isnan(Element))) = [];
Element(find(Element==0)) = [];
nonsigROI=setdiff(Element,ROI_index);
for i=1:size(nonsigROI,2)
    
   MaskData(find(MaskData==nonsigROI(i)))=0;
end
Tvalue=stats.tstat(ROI_index);
for j=1:size(ROI_index,2)
    
   MaskData(find(MaskData==ROI_index(j)))=Tvalue(j);
end
mint=min(Tvalue)
maxt=max(Tvalue)
y_Write(MaskData,MaskHead,[outpath,'sigpartialspearmanROIresults_Tmap'])
