% Organize left and right RMPFC as well as left and right MTL ROI maps
% and use paired t test to compare FC patterns between rumination and
% distraction

%organize data
clear;clc;
DataDir = '/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing';
ResultDir = '/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis4Publish/ROI_Analysis/DataMaps';

SubList = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/Analysis/IPCAS_Sublist.txt');
SiteSet = {'IPCAS','PKUGE','PKUSIMENS'};
ROINameSet = {'left_RMPFC','right_RMPFC','left_DMPFC','right_DMPFC','left_MTL','right_MTL'};
ROISet = {'ROI4','ROI9','ROI12','ROI17','ROI21','ROI24'};
ConditionSet = {'S3_','S4_'};
ConditionNameSet = {'Rumination','Distraction'};

for iSite = 1:length(SiteSet)
    for iCondition = 1:length(ConditionSet)
       for iROI = 1: length(ROINameSet)
            for iSubject = 1:length(SubList)
                CurrentFilePath = [DataDir,'/',SiteSet{iSite},'_task/',ConditionSet{iCondition},'Results/FC_FunImgARCWFS/',ROISet{iROI},'FCMap_',SubList{iSubject},'.nii'];
                TargetFilePath = [ResultDir,'/',SiteSet{iSite},'/',ConditionNameSet{iCondition},'/',ROINameSet{iROI}];
                if ~exist(TargetFilePath,'dir'); mkdir(TargetFilePath); end
                copyfile(CurrentFilePath,TargetFilePath);
            end
       end
    end
end

MaskFile = '/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Regional_Metrics/Allsite_Greymatter_GroupMask.nii';
OutputDir = '/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis4Publish/ROI_Analysis';

PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=1;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'

PALMSettings31=PALMSettings;
PALMSettings31.ClusterFormingThreshold=3.1;
PALMSettings258=PALMSettings;
PALMSettings258.ClusterFormingThreshold=2.58;
PALMSettings329=PALMSettings;
PALMSettings329.ClusterFormingThreshold=3.29;

for iSite = 1:length(SiteSet)
    for iROI = 1: length(ROINameSet)
        Dependentdirs{1} = [ResultDir,'/',SiteSet{iSite},'/Rumination/',ROINameSet{iROI}];
        Dependentdirs{2} = [ResultDir,'/',SiteSet{iSite},'/Distraction/',ROINameSet{iROI}];
        Outputname = [OutputDir,'/',SiteSet{iSite},'_',ROINameSet{iROI}];
        [TTestPaired_T,Header] = y_TTestPaired_Image(Dependentdirs,Outputname,MaskFile,[],[],PALMSettings);
    end
end