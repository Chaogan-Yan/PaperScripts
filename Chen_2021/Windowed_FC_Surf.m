% To output the dynamic ROI-level FC maps (n x n x timepoints)
% n is the number of ROIs
% Originally written by ChenXiao on 20181130
% Modified to deal with surf data on 20210125

clear all;clc;

DataDir = '/mnt/Data3/RfMRILab/ChenX/Rumination_surf/';
ResultDir = '/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/ROI_FCMatrix_Dynamic/';
Sub_List = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/Sublist_NoLeftHanded.txt');
Condition_Set = {'rum','dis'};
Site_Set = {'IPCAS','PKUGE','PKUSIEMENS'};
ROI_Path = '/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/PairedT_Standardized_NoLeftHanded/Sig_ROIs';

ROI_Cell_lh = {'lh_DLPFC','lh_hippocampus'};
ROI_Cell_rh = {'rh_hippocampus', 'rh_MCC', 'rh_MPFC', 'rh_Occipital', 'rh_PostCentral', ...
                'rh_TPJdor', 'rh_TPJros'};
            
ROI_Cell_full = {'lh_DLPFC','lh_hippocampus','rh_hippocampus', 'rh_MCC', 'rh_MPFC', 'rh_Occipital', 'rh_PostCentral', ...
                'rh_TPJdor', 'rh_TPJros'};

WindowSize = 32;
WindowStep = 2;
WindowType = 'rectwin';

%%%%%%%%%%%%%%%%%%% Extract roi signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dealing with the left hemisphere
parpool(20);
SeedSeries_rum = cell(3,length(Sub_List));
SeedSeries_dis = cell(3,length(Sub_List));

%IPCAS
for iROI = 1:length(ROI_Cell_lh)
    ROIDef{1,1} = [ROI_Path,'/PKUSIEMENS_',ROI_Cell_lh{iROI},'_Mask.gii'];
    parfor iSub = 1:length(Sub_List)
        Img = [DataDir,'IPCAS_task/S3_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-3_task-rest_space-fsaverage5_hemi-L.func.gii'];
        SeedSeries_rum{1,iSub}(:,iROI) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
        Img = [DataDir,'IPCAS_task/S4_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-4_task-rest_space-fsaverage5_hemi-L.func.gii'];
        SeedSeries_dis{1,iSub}(:,iROI) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
    end
end
%PKUGE
for iROI = 1:length(ROI_Cell_lh)
    ROIDef{1,1} = [ROI_Path,'/PKUSIEMENS_',ROI_Cell_lh{iROI},'_Mask.gii'];
    parfor iSub = 1:length(Sub_List)
        Img = [DataDir,'PKUGE_task/S3_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-3_task-rest_space-fsaverage5_hemi-L.func.gii'];
        SeedSeries_rum{2,iSub}(:,iROI) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
        Img = [DataDir,'PKUGE_task/S4_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-4_task-rest_space-fsaverage5_hemi-L.func.gii'];
        SeedSeries_dis{2,iSub}(:,iROI) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
    end
end
%PKUSIEMENS
for iROI = 1:length(ROI_Cell_lh)
    ROIDef{1,1} = [ROI_Path,'/PKUSIEMENS_',ROI_Cell_lh{iROI},'_Mask.gii'];
    parfor iSub = 1:length(Sub_List)
        Img = [DataDir,'PKUSIEMENS_task/S3_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-3_task-rest_space-fsaverage5_hemi-L.func.gii'];
        SeedSeries_rum{3,iSub}(:,iROI) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
        Img = [DataDir,'PKUSIEMENS_task/S4_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-4_task-rest_space-fsaverage5_hemi-L.func.gii'];
        SeedSeries_dis{3,iSub}(:,iROI) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
    end
end

%dealing with the right hemisphere
%IPCAS
for iROI = 1:length(ROI_Cell_rh)
    ROIDef{1,1} = [ROI_Path,'/PKUSIEMENS_',ROI_Cell_rh{iROI},'_Mask.gii'];
    parfor iSub = 1:length(Sub_List)
        Img = [DataDir,'IPCAS_task/S3_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-3_task-rest_space-fsaverage5_hemi-R.func.gii'];
        SeedSeries_rum{1,iSub}(:,iROI+length(ROI_Cell_lh)) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
        Img = [DataDir,'IPCAS_task/S4_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-4_task-rest_space-fsaverage5_hemi-R.func.gii'];
        SeedSeries_dis{1,iSub}(:,iROI+length(ROI_Cell_lh)) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
    end
end
%PKUGE
for iROI = 1:length(ROI_Cell_rh)
    ROIDef{1,1} = [ROI_Path,'/PKUSIEMENS_',ROI_Cell_rh{iROI},'_Mask.gii'];
    parfor iSub = 1:length(Sub_List)
        Img = [DataDir,'PKUGE_task/S3_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-3_task-rest_space-fsaverage5_hemi-R.func.gii'];
        SeedSeries_rum{2,iSub}(:,iROI+length(ROI_Cell_lh)) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
        Img = [DataDir,'PKUGE_task/S4_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-4_task-rest_space-fsaverage5_hemi-R.func.gii'];
        SeedSeries_dis{2,iSub}(:,iROI+length(ROI_Cell_lh)) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
    end
end
%PKUSIEMENS
for iROI = 1:length(ROI_Cell_rh)
    ROIDef{1,1} = [ROI_Path,'/PKUSIEMENS_',ROI_Cell_rh{iROI},'_Mask.gii'];
    parfor iSub = 1:length(Sub_List)
        Img = [DataDir,'PKUSIEMENS_task/S3_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-3_task-rest_space-fsaverage5_hemi-R.func.gii'];
        SeedSeries_rum{3,iSub}(:,iROI+length(ROI_Cell_lh)) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
        Img = [DataDir,'PKUSIEMENS_task/S4_FunSurfWCFS/sub-',Sub_List{iSub},'/ssub-',Sub_List{iSub}, ...
            '_ses-4_task-rest_space-fsaverage5_hemi-R.func.gii'];
        SeedSeries_dis{3,iSub}(:,iROI+length(ROI_Cell_lh)) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);
    end
end

nDimTimePoints = 235;
nROI = 9;
ROICorrelationDynamic = [];
ROICorrelationDynamic_std_rum = {};
for iSite = 1:3
    for iSub = 1:length(Sub_List)
        %rum
        SeedSeries = SeedSeries_rum{iSite,iSub};

        nWindow = fix((nDimTimePoints - WindowSize)/WindowStep) + 1;
        nDimTimePoints_WithinWindow = WindowSize;

        if ischar(WindowType)
            eval(['WindowType = ',WindowType,'(WindowSize);'])
        end
        WindowMultiplier = repmat(WindowType(:),1,nROI);
        for iWindow = 1:nWindow
            fprintf('\nProcessing window %g of total %g windows\n', iWindow, nWindow);
            SeedSeriesWindow = SeedSeries((iWindow-1)*WindowStep+1:(iWindow-1)*WindowStep+WindowSize,:);
            SeedSeriesWindow = SeedSeriesWindow.*repmat(WindowType(:),1,size(SeedSeriesWindow,2));
            SeedSeriesWindow = SeedSeriesWindow-repmat(mean(SeedSeriesWindow),size(SeedSeriesWindow,1),1);

            ROICorrelationWindow = corrcoef(SeedSeriesWindow);
            ROICorrelation_FisherZ_Window = 0.5 * log((1 + ROICorrelationWindow)./(1- ROICorrelationWindow));
            ROICorrelation_FisherZ_Window(find(isinf(ROICorrelation_FisherZ_Window))) = 0;
            ROICorrelationDynamic(:,:,iWindow) = ROICorrelation_FisherZ_Window;
        end
        
        dFCSerries_rh_MPFC_rh_hippocampus_rum{iSite}(:,iSub) = ROICorrelationDynamic(3,5,:);
        dFCSerries_lh_hippocampus_rh_TPJros_rum{iSite}(:,iSub) = ROICorrelationDynamic(2,9,:);
        
        ROICorrelationDynamic_std = std(ROICorrelationDynamic,0,3);
        ROICorrelationDynamic_std_rum{iSite,iSub} = ROICorrelationDynamic_std;
        
        %dis
        SeedSeries = SeedSeries_dis{iSite,iSub};

        nWindow = fix((nDimTimePoints - WindowSize)/WindowStep) + 1;
        nDimTimePoints_WithinWindow = WindowSize;

        if ischar(WindowType)
            eval(['WindowType = ',WindowType,'(WindowSize);'])
        end
        WindowMultiplier = repmat(WindowType(:),1,nROI);
        for iWindow = 1:nWindow
            fprintf('\nProcessing window %g of total %g windows\n', iWindow, nWindow);
            SeedSeriesWindow = SeedSeries((iWindow-1)*WindowStep+1:(iWindow-1)*WindowStep+WindowSize,:);
            SeedSeriesWindow = SeedSeriesWindow.*repmat(WindowType(:),1,size(SeedSeriesWindow,2));
            SeedSeriesWindow = SeedSeriesWindow-repmat(mean(SeedSeriesWindow),size(SeedSeriesWindow,1),1);

            ROICorrelationWindow = corrcoef(SeedSeriesWindow);
            ROICorrelation_FisherZ_Window = 0.5 * log((1 + ROICorrelationWindow)./(1- ROICorrelationWindow));
            ROICorrelation_FisherZ_Window(find(isinf(ROICorrelation_FisherZ_Window))) = 0;
            ROICorrelationDynamic(:,:,iWindow) = ROICorrelation_FisherZ_Window;
        end
        
        dFCSerries_rh_MPFC_rh_hippocampus_dis{iSite}(:,iSub) = ROICorrelationDynamic(3,5,:);
        dFCSerries_lh_hippocampus_rh_TPJros_dis{iSite}(:,iSub) = ROICorrelationDynamic(2,9,:);
        
        ROICorrelationDynamic_std = std(ROICorrelationDynamic,0,3);
        ROICorrelationDynamic_std_dis{iSite,iSub} = ROICorrelationDynamic_std;
    end
    
    dlmwrite([ResultDir,Site_Set{iSite},'_rh_MPFC-rh_hippocampus_rum_dFCSerries.csv'], ...
                    dFCSerries_rh_MPFC_rh_hippocampus_rum{iSite},'delimiter',',');
    dlmwrite([ResultDir,Site_Set{iSite},'_rh_MPFC-rh_hippocampus_dis_dFCSerries.csv'], ...
                    dFCSerries_rh_MPFC_rh_hippocampus_dis{iSite},'delimiter',',');
                
    dlmwrite([ResultDir,Site_Set{iSite},'_lh_hippocampus_rh_TPJros_rum_dFCSerries.csv'], ...
                    dFCSerries_lh_hippocampus_rh_TPJros_rum{iSite},'delimiter',',');
    dlmwrite([ResultDir,Site_Set{iSite},'_lh_hippocampus_rh_TPJros_dis_dFCSerries.csv'], ...
                    dFCSerries_lh_hippocampus_rh_TPJros_dis{iSite},'delimiter',',');
end



% read in head motion
% 1: rum 2: dis
Head_Motion = {};
for iSite = 1:length(Site_Set)
    for i=1:length(Sub_List)
        Temp=load(['/mnt/Data3/RfMRILab/ChenX/Rumination_surf/',Site_Set{iSite},'_task/RealignParameter/sub-',Sub_List{i},'/S3_FD_Jenkinson_sub-',Sub_List{i},'.txt']);
        Head_Motion{iSite}(i,1)=mean(Temp);
        Temp=load(['/mnt/Data3/RfMRILab/ChenX/Rumination_surf/',Site_Set{iSite},'_task/RealignParameter/sub-',Sub_List{i},'/S4_FD_Jenkinson_sub-',Sub_List{i},'.txt']);
        Head_Motion{iSite}(i,2)=mean(Temp);
    end
end

%generating subject regressors 
nSub = length(Sub_List);
Regressors = [ones(nSub,1);-1*ones(nSub,1)];
for i=1:nSub
    SubjectRegressors(:,i) = zeros(nSub*2,1);
    SubjectRegressors(i:nSub:nSub*2,i) = 1;
end
OtherCovariatesMatrix = [];

% do statistical testing
Stats = {};
Stats{1,2} = 'IPCAS';
Stats{1,5} = 'PKUGE';
Stats{1,8} = 'PKUSIEMENS';
for iSite = 1:3
    iCount = 0;
    for iROI = 1:8
        for jROI = iROI+1:9
            iCount = iCount + 1;
            Stats{iCount+1,1} = [ROI_Cell_full{iROI},'-',ROI_Cell_full{jROI}];
            
            y_rum = [];y_dis = [];
            for iSub = 1:length(Sub_List)
                y_rum(iSub,1) = ROICorrelationDynamic_std_rum{iSite, iSub}(iROI,jROI);
                y_dis(iSub,1) = ROICorrelationDynamic_std_dis{iSite, iSub}(iROI,jROI);
            end
            y = [y_rum;y_dis];
            OtherCovariatesMatrix = [];
            OtherCovariatesMatrix = [Head_Motion{iSite}(:,1);Head_Motion{iSite}(:,2)];
            FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
            Contrast = zeros(1,size(FullRegressors,2));
            Contrast(1) = 1;
            TF_Flag = 'T';
            [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
            Stats{iCount+1,iSite*3-1} = TF_ForContrast;
            if 2*tcdf(-abs(TF_ForContrast),nSub-2) < 0.05
                Stats{iCount+1,iSite*3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
            end
            p4fdr(iCount,1) = 2*tcdf(-abs(TF_ForContrast),nSub-2);
            
            if (iROI == 2) && (jROI == 9)
                %write out csv files for violin plots
                adjusted = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
                adjusted(1:nSub,2) = 1;
                adjusted(nSub+1:2*nSub,2) = 2;
                dlmwrite([ResultDir,Site_Set{iSite},'_',ROI_Cell_full{iROI},'-',ROI_Cell_full{jROI},'_rum_dis_adjustedROIstd.csv'], ...
                    adjusted,'delimiter',',');
            end
            
            if (iROI == 3) && (jROI == 5)
                %write out csv files for violin plots
                adjusted = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
                adjusted(1:nSub,2) = 1;
                adjusted(nSub+1:2*nSub,2) = 2;
                dlmwrite([ResultDir,Site_Set{iSite},'_',ROI_Cell_full{iROI},'-',ROI_Cell_full{jROI},'_rum_dis_adjustedROIstd.csv'], ...
                    adjusted,'delimiter',',');
            end
        end
    end
    pfdr = mafdr(p4fdr,'BHFDR',true);
    pfdr(find(pfdr>0.05)) = 0;
    iCount = 0;
    for iROI = 1:8
        for jROI = iROI+1:9
            iCount = iCount + 1;
            Stats{iCount+1,iSite*3+1} = pfdr(iCount,1);
        end
    end
end

save([ResultDir,'std_Stats.mat'],'Stats');
