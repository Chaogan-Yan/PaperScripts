%Perform paired t test on the ROI level (within DMN)
%Controlling head motion info
%20190620 ChenXiao

clear;clc;

%initialization
Workdir='/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing';
Resultdir='/MD3860F/RfMRILab/ChenX/Rumination_project/Analysis/Analysis4Publish';
SubList = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/Analysis/IPCAS_Sublist.txt');

SiteSet = {'IPCAS','PKUGE','PKUSIMENS'};

%Extract head motion
%1,Rest; 2, Happy; 3, Sad; 4, Rum; 5, Dis;
HeadMotion = {};
for iSite = 1:length(SiteSet)
    for i=1:length(SubList)
        Temp=load(['/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/',SiteSet{iSite},'_rest/RealignParameter/',SubList{i},'/FD_Jenkinson_',SubList{i},'.txt']);
        HeadMotion{iSite}(i,1)=mean(Temp);
        Temp=load(['/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/',SiteSet{iSite},'_task/RealignParameter/',SubList{i},'/FD_Jenkinson_',SubList{i},'.txt']);
        HeadMotion{iSite}(i,2)=mean(Temp);
        Temp=load(['/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/',SiteSet{iSite},'_task/RealignParameter/',SubList{i},'/S2_FD_Jenkinson_',SubList{i},'.txt']);
        HeadMotion{iSite}(i,3)=mean(Temp);
        Temp=load(['/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/',SiteSet{iSite},'_task/RealignParameter/',SubList{i},'/S3_FD_Jenkinson_',SubList{i},'.txt']);
        HeadMotion{iSite}(i,4)=mean(Temp);
        Temp=load(['/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/',SiteSet{iSite},'_task/RealignParameter/',SubList{i},'/S4_FD_Jenkinson_',SubList{i},'.txt']);
        HeadMotion{iSite}(i,5)=mean(Temp);
    end
end

for iSite = 1:length(SiteSet)

    RestFCDir=['/',SiteSet{iSite},'_rest/Results/ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
    HappyFCDir=['/',SiteSet{iSite},'_task/Results/ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
    SadFCDir=['/',SiteSet{iSite},'_task/S2_Results/S2_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
    RumFCDir=['/',SiteSet{iSite},'_task/S3_Results/S3_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
    DisFCDir=['/',SiteSet{iSite},'_task/S4_Results/S4_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];

    FCDirset = {RestFCDir, HappyFCDir,SadFCDir,RumFCDir,DisFCDir};



    AllFCMatrix = zeros(41,5);
    rum_rest_tMap = [];
    rum_rest_pMap = [];
    

    rum_dis_tMap = [];
    rum_dis_pMap = [];
    
    % Using y_regress_ss to conduct t test with GLM
    nSub = length(SubList);
    Regressors = [ones(nSub,1);-1*ones(nSub,1)];
    for i=1:nSub
        SubjectRegressors(:,i) = zeros(nSub*2,1);
        SubjectRegressors(i:nSub:nSub*2,i) = 1;
    end
    OtherCovariatesMatrix = [];
    
    for iROI = 1:23
       for jROI = iROI+1:24
            for iFC = 1:length(FCDirset)
                for i = 1:length(SubList)
                   load([Workdir,FCDirset{iFC},SubList{i},'.mat']);
                   AllFCMatrix(i,iFC) = ROICorrelation_FisherZ(iROI,jROI);
                end 
            end

            y = [AllFCMatrix(:,4);AllFCMatrix(:,5)];
            OtherCovariatesMatrix = [HeadMotion{iSite}(:,4);HeadMotion{iSite}(:,5)];
            FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
            Contrast = zeros(1,size(FullRegressors,2));
            Contrast(1) = 1;
            TF_Flag = 'T';
            [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
            rum_dis_tMap(iROI,jROI) = TF_ForContrast;
            rum_dis_pMap(iROI,jROI) = 2*tcdf(-abs(TF_ForContrast),nSub-1);

            y = [AllFCMatrix(:,4);AllFCMatrix(:,1)];
            OtherCovariatesMatrix = [HeadMotion{iSite}(:,4);HeadMotion{iSite}(:,1)];
            FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
            Contrast = zeros(1,size(FullRegressors,2));
            Contrast(1) = 1;
            TF_Flag = 'T';
            [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
            rum_rest_tMap(iROI,jROI) = TF_ForContrast;
            rum_rest_pMap(iROI,jROI) = 2*tcdf(-abs(TF_ForContrast),nSub-1);
            
            AllFCMatrix = zeros(41,5);
       end
    end
    %fdr correction
    icount = 0;
    for i = 1:23
        for j = i+1:24
            icount = icount +1;
            rum_dis_pMap4fdr(icount,1) = rum_dis_pMap(i,j);
            rum_rest_pMap4fdr(icount,1) = rum_rest_pMap(i,j);
        end
    end

    rum_dis_pMapAfterfdr = mafdr(rum_dis_pMap4fdr,'BHFDR',true);
    rum_rest_pMapAfterfdr = mafdr(rum_rest_pMap4fdr,'BHFDR',true);
        
    icount = 0;
    for i = 1:23
        for j = i+1:24
            icount = icount +1;
            rum_dis_fdrcorrected_pMap(i,j) = rum_dis_pMapAfterfdr(icount);
            rum_rest_fdrcorrected_pMap(i,j) = rum_rest_pMapAfterfdr(icount);
        end
    end
    save([Resultdir,'/',SiteSet{iSite},'_ROI_Level_Stats_Smooth.mat'],'*Map');  
end








