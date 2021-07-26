% Paried t tests comparing rum vs. dis and rum vs. rest on stability

clear;clc;
% Sub_List = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/Analysis/IPCAS_Sublist.txt');
% no Left handed
Sub_List = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/Sublist_NoLeftHanded.txt');
Data_Dir = '/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Stability_Surf/smoothed_WL_64';
Metric_Dir = '/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Stability_Surf/Organized_Metrics_NoLeftHanded_WL_64';
Metric_Standardized_Dir = '/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Stability_Surf/Organized_Metrics_NoLeftHanded_WL_64_Standardized';
Result_Dir = '/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/PairedT_Standardized_NoLeftHanded_WL_64';
Temp_Dir = '/mnt/Data/RfMRILab/ChenX/CX_software/DPABI_V5.1_201230/DPABISurf/SurfTemplates';
Hemisphere_Set = {'LH','RH'};
HemisphereName_Set = {'lh','rh'};
Site_Set  = {'IPCAS','PKUGE','PKUSIEMENS'};
Condition_Set = {'rum','dis','rest'};

% organize data
for iSite = 1:3
    for iCondition = 1:3
        targetdir_lh = [Metric_Dir,'/',Site_Set{iSite},'/',Condition_Set{iCondition},'/lh'];
        if ~exist(targetdir_lh); mkdir(targetdir_lh); end
        targetdir_rh = [Metric_Dir,'/',Site_Set{iSite},'/',Condition_Set{iCondition},'/rh'];
        if ~exist(targetdir_rh); mkdir(targetdir_rh); end
        for iSubject = 1:length(Sub_List)
            originfile_lh = [Data_Dir,'/',Site_Set{iSite},'/',Condition_Set{iCondition},'/',Sub_List{iSubject},'-LH.gii'];
            copyfile(originfile_lh,targetdir_lh);
            originfile_rh = [Data_Dir,'/',Site_Set{iSite},'/',Condition_Set{iCondition},'/',Sub_List{iSubject},'-RH.gii'];
            copyfile(originfile_rh,targetdir_rh);
        end
    end
end

% standardization
MaskData_LH = [Temp_Dir,'/fsaverage5_lh_cortex.label.gii'];
MaskData_RH = [Temp_Dir,'/fsaverage5_rh_cortex.label.gii'];
MethodType = 4; %Z - Standardization
for iSite = 1:3
    for iCondition = 1:3
        ImgCells_LH = {};ImgCells_RH = {};
        for iSub = 1:length(Sub_List)
            ImgCells_LH{iSub,1} = [Metric_Dir,'/',Site_Set{iSite},'/',Condition_Set{iCondition},'/lh/',Sub_List{iSub},'-LH.gii'];
            ImgCells_RH{iSub,1} = [Metric_Dir,'/',Site_Set{iSite},'/',Condition_Set{iCondition},'/rh/',Sub_List{iSub},'-RH.gii'];
        end
        Output_Dir = [Metric_Standardized_Dir,'/',Site_Set{iSite},'/',Condition_Set{iCondition}];
        if ~exist(Metric_Standardized_Dir); mkdir(Metric_Standardized_Dir); end
        y_Standardization_Surf(ImgCells_LH, ImgCells_RH, MaskData_LH, MaskData_RH, MethodType, Output_Dir, '');
    end
end

% read in head motion
% 1: rum 2: dis 3: rest
Head_Motion = {};
for iSite = 1:length(Site_Set)
    for i=1:length(Sub_List)
        Temp=load(['/mnt/Data3/RfMRILab/ChenX/Rumination_surf/',Site_Set{iSite},'_task/RealignParameter/sub-',Sub_List{i},'/S3_FD_Jenkinson_sub-',Sub_List{i},'.txt']);
        Head_Motion{iSite}(i,1)=mean(Temp);
        Temp=load(['/mnt/Data3/RfMRILab/ChenX/Rumination_surf/',Site_Set{iSite},'_task/RealignParameter/sub-',Sub_List{i},'/S4_FD_Jenkinson_sub-',Sub_List{i},'.txt']);
        Head_Motion{iSite}(i,2)=mean(Temp);
        Temp=load(['/mnt/Data3/RfMRILab/ChenX/Rumination_surf/',Site_Set{iSite},'_rest/RealignParameter/sub-',Sub_List{i},'/FD_Jenkinson_sub-',Sub_List{i},'.txt']);
        Head_Motion{iSite}(i,3)=mean(Temp);
    end
end

PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=0;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'

% % read scale and demographical data
% load '/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Demographic_data/DemographicInfo_V3.mat';
% ScaleMat = [Rumination,Brooding,Reflection];
% DemographicMat = [Age,Sex];
% 
% %exclude left handed sub010
% wantedSubMatrix = ones(length(Sublist),1);
% %the variable Sublist here is stored in the "DemographicInfo_V2.mat", which contains the Sub010, not
% %the Sub_List above!
% wantedSubMatrix(10,1) = 0;
% WantedSubIndex = find(wantedSubMatrix);
% Sublist = Sublist(WantedSubIndex);
% %now the Sublist is equal to Sub_List
% ScaleMat = ScaleMat(WantedSubIndex,:);
% DemographicMat = DemographicMat(WantedSubIndex,:);
% avg_Demographic = mean(DemographicMat);
% maxage = max(DemographicMat(:,1));
% minage = min(DemographicMat(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%One sample t test%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Base = 0;
CovariateDirs = [];
for iSite = 1:3
    TopOutput_Dir = ['/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/OneSampleT_Standardized_NoLeftHanded/',Site_Set{iSite}];
    if ~exist(TopOutput_Dir) mkdir(TopOutput_Dir); end
    for iCondition = 1:3
        for iHem = 1:2
            PALMSettings.SurfFile = [Temp_Dir,'/fsaverage5_',HemisphereName_Set{iHem},'_white.surf.gii'];
            PALMSettings.SurfAreaFile = [Temp_Dir,'/fsaverage5_',HemisphereName_Set{iHem},'_white_avg.area.gii'];
            MaskFile =  [Temp_Dir,'/fsaverage5_',HemisphereName_Set{iHem},'_cortex.label.gii'];
            DependentDirs{1,1} = [Metric_Standardized_Dir,'/',Site_Set{iSite},'/',Condition_Set{iCondition},'/',Hemisphere_Set{iHem},'/',HemisphereName_Set{iHem}];
            OutputName = [TopOutput_Dir,'/',Condition_Set{iCondition},'_',Hemisphere_Set{iHem}];
            OtherCovariates{1,1} = [DemographicMat,Head_Motion{iSite}(:,iCondition)];
            [TTest1_T,Header] = y_TTest1_Image(DependentDirs,OutputName,MaskFile,CovariateDirs,OtherCovariates,Base,PALMSettings);
        end
    end
end
% write out corrected data
for iSite = 1:3
    TopOutput_Dir = ['/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/OneSampleT_Standardized_NoLeftHanded/',Site_Set{iSite}];
    for iCondition = 1:3
        for iHem = 1:2
            [Data T T Header]=y_ReadAll([TopOutput_Dir,'/',Condition_Set{iCondition},'_',Hemisphere_Set{iHem},'_dpv_tstat.gii']);
            P=y_ReadAll([TopOutput_Dir,'/',Condition_Set{iCondition},'_',Hemisphere_Set{iHem},'_tfce_tstat_fwep.gii']);
            NewData=Data.*(P<=0.025);
            y_Write(NewData,Header,[TopOutput_Dir,'/',Condition_Set{iCondition},'_',Hemisphere_Set{iHem},'_dpv_tstat_TFCECorrected025.gii']);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Paired t test%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%rum vs dis
TopOutput_Dir = [Result_Dir,'/rumvsdis/'];
if ~exist(TopOutput_Dir); mkdir(TopOutput_Dir); end
for iSite = 1:3
    OtherCovariates = {};
    OtherCovariates{1,1} = Head_Motion{iSite}(:,1);%rum
    OtherCovariates{2,1} = Head_Motion{iSite}(:,2);%dis
     for iHem = 1:2
        PALMSettings.SurfFile = [Temp_Dir,'/fsaverage5_',HemisphereName_Set{iHem},'_white.surf.gii'];
        PALMSettings.SurfAreaFile = [Temp_Dir,'/fsaverage5_',HemisphereName_Set{iHem},'_white_avg.area.gii'];
        MaskFile =  [Temp_Dir,'/fsaverage5_',HemisphereName_Set{iHem},'_cortex.label.gii'];
        DependentDirs{1,1} = [Metric_Standardized_Dir,'/',Site_Set{iSite},'/rum/',Hemisphere_Set{iHem},'/',HemisphereName_Set{iHem}];
        DependentDirs{2,1} = [Metric_Standardized_Dir,'/',Site_Set{iSite},'/dis/',Hemisphere_Set{iHem},'/',HemisphereName_Set{iHem}];
        OutputName = [TopOutput_Dir,'/',Site_Set{iSite},'_',Hemisphere_Set{iHem}];
        [TTestPaired_T,Header] = y_TTestPaired_Image(DependentDirs,OutputName,MaskFile,[],OtherCovariates,PALMSettings);
     end
end

%rum vs rest
TopOutputDir = [Result_Dir,'/rumvsrest/'];
if ~exist(TopOutputDir); mkdir(TopOutputDir); end
for iSite = 1:3
    OtherCovariates = {};
    OtherCovariates{1,1} = Head_Motion{iSite}(:,1);%rum
    OtherCovariates{2,1} = Head_Motion{iSite}(:,3);%rest
     for iHem = 1:2
        PALMSettings.SurfFile = [Temp_Dir,'/fsaverage5_',HemisphereName_Set{iHem},'_white.surf.gii'];
        PALMSettings.SurfAreaFile = [Temp_Dir,'/fsaverage5_',HemisphereName_Set{iHem},'_white_avg.area.gii'];
        MaskFile =  [Temp_Dir,'/fsaverage5_',HemisphereName_Set{iHem},'_cortex.label.gii'];
        DependentDirs{1,1} = [Metric_Standardized_Dir,'/',Site_Set{iSite},'/rum/',Hemisphere_Set{iHem},'/',HemisphereName_Set{iHem}];
        DependentDirs{2,1} = [Metric_Standardized_Dir,'/',Site_Set{iSite},'/rest/',Hemisphere_Set{iHem},'/',HemisphereName_Set{iHem}];
        OutputName = [TopOutputDir,'/',Site_Set{iSite},'_',Hemisphere_Set{iHem}];
        [TTestPaired_T,Header] = y_TTestPaired_Image(DependentDirs,OutputName,MaskFile,[],OtherCovariates,PALMSettings);
     end
end