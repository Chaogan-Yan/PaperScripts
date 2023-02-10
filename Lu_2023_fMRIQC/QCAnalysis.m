%% Download data

%% Unzip data

%% REST - organize data to DPABI input format
clc;clear;
InDir = '/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/Data/fmri-open-qc-rest';
OutDir = '/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/ForDPABI/Rest_Surf_New';
mkdir(OutDir);

SubList = dir([InDir,filesep,'sub*']);
for iSub = 1:length(SubList)
    T1Dir = [OutDir,filesep,'T1Img',filesep,SubList(iSub).name];
    mkdir(T1Dir);
    FunDir = [OutDir,filesep,'FunImg',filesep,SubList(iSub).name]; 
    mkdir(FunDir);

%     copyfile([InDir,filesep,SubList(iSub).name,filesep,'ses-01',filesep,'anat',filesep,SubList(iSub).name,'*run-01*'],T1Dir);
    try
        copyfile([InDir,filesep,SubList(iSub).name,filesep,'ses-01',filesep,'func',filesep,SubList(iSub).name,'*run-02*'],FunDir);% sub-608 has two run in func ses-01
    catch
        copyfile([InDir,filesep,SubList(iSub).name,filesep,'ses-01',filesep,'func',filesep,SubList(iSub).name,'*run-01*'],FunDir);% sub-608 has two run in func ses-01
    end
    disp(['Sub ',num2str(iSub)])
end


%% REST_Volume - organize data to DPABI input format
clc;clear;
InDir = '/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/Data/fmri-open-qc-rest';
OutDir = '/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/ForDPABI/Rest_Volume';
mkdir(OutDir);

SubList = dir([InDir,filesep,'sub*']);
for iSub = 1:length(SubList)
    T1Dir = [OutDir,filesep,'T1Img',filesep,SubList(iSub).name];
    mkdir(T1Dir);
    FunDir = [OutDir,filesep,'FunImg',filesep,SubList(iSub).name]; 
    mkdir(FunDir);

%     copyfile([InDir,filesep,SubList(iSub).name,filesep,'ses-01',filesep,'anat',filesep,SubList(iSub).name,'*run-01*'],T1Dir);
    try
        copyfile([InDir,filesep,SubList(iSub).name,filesep,'ses-01',filesep,'func',filesep,SubList(iSub).name,'*run-02*'],FunDir);% sub-608 has two run in func ses-01
    catch
        copyfile([InDir,filesep,SubList(iSub).name,filesep,'ses-01',filesep,'func',filesep,SubList(iSub).name,'*run-01*'],FunDir);% sub-608 has two run in func ses-01
    end
    disp(['Sub ',num2str(iSub)])
end


%% DPABISurf Preprocessing 

% Error1, sub-608 has two run in func ses-01

% Error2, time points were not consistent

% Error3, Site 2 and 5 failed to run fmriprep
% failed BIDS validator

% report from  http://bids-standard.github.io/bids-validator/
% "SliceTiming" value/s contains invalid value as it is greater than RepetitionTime. SliceTiming values should be in seconds not milliseconds (common mistake).

% After re-run fmriprep, these subjects have not run fmriprep yet:
%     {'sub-501'}
%     {'sub-502'}
%     {'sub-503'}
%     {'sub-504'}
%     {'sub-509'}

% Error4, Some sites don't have slice-timing information and some sites have WRONG slice-timing information.
% cancel slice-timing correction


%% Compute sex-difference in the functional metrics - with QC
clc;clear;
load('/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/ForDPABI/Rest_Volume/Pred_Sex.mat');
load('/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/ForDPABI/Rest_Surf_New/Motion_WithID.mat');
InDir = '/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/ForDPABI/Rest_Surf_New/ResultsS/';
OutDir = '/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/Statistic/Rest_Surf_QC/';

[C,IA,IB] = intersect(ID,Motion_ID);
SubID = C;
Sex = round(Pred_Sex(IA));
Sex(find(Sex==0)) = -1;
Motion = Motion_FD(IB);
Site = cellfun(@(ID) str2num(ID(5)),SubID);


Metrics  = {'ALFF_FunSurfWC','fALFF_FunSurfWC','ReHo_FunSurfWCF','DegreeCentrality_FunSurfWCF'};
Prefixs = {'szALFF_','szfALFF_','szReHo_','szDegreeCentrality_Bilateral_PositiveWeightedSumBrain_'};
Hemispheres = {'FunSurfLH','FunSurfRH','FunVolu'};
Masks = {'/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V5.1_201201/DPABISurf/SurfTemplates/fsaverage5_lh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V5.1_201201/DPABISurf/SurfTemplates/fsaverage5_rh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/Reslice_freesurfer_subcortical_mask_1mm.nii'};


SiteCov=[];
SiteIndex = unique(Site);
for i=1:length(SiteIndex)-1
    SiteCov(:,i) = Site==SiteIndex(i);
end

BadImageList = {
    'sub-102',...
    'sub-114',...
    'sub-115',...
    'sub-203',...
    'sub-305',...
    'sub-307',...
    'sub-312',...
    'sub-315',...
    'sub-501',...
    'sub-502',...
    'sub-503',...
    'sub-504',...
    'sub-508',...
    'sub-509',...
    'sub-518',...
    'sub-519',...
    'sub-620',...
    'sub-706',...
    'sub-714',...
    'sub-719',...
    };
QCFlag = ones(length(SubID),1);
[C,IA,IB] = intersect(SubID,BadImageList);
QCFlag(IA) = 0;
SubIndex = find(QCFlag);

for iHemi = 1:length(Hemispheres)
    for iMetric = 1:length(Metrics)
        if iHemi < 3
            DependentVolume = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric},filesep,...
                Prefixs{iMetric},ID,'.func.gii'],SubID(SubIndex),'UniformOutput',false);
        else
            DependentVolume = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,regexprep(Metrics{iMetric},'Surf','Volu'),...
                filesep,regexprep(Prefixs{iMetric},'_Bilateral_','_'),ID,'.nii'],SubID(SubIndex),'UniformOutput',false);
        end
        Predictor = [ones(length(SubIndex),1),Sex(SubIndex),Motion(SubIndex),SiteCov(SubIndex,:)];
        OutputDir = [OutDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric}];
        mkdir(OutputDir)
        Contrast = zeros(1,size(Predictor,2));
        Contrast(2) = 1;
        y_GroupAnalysis_Image(DependentVolume,Predictor,[OutputDir,filesep,'T2'],Masks{iHemi},[],Contrast,'T');
        disp(['Have down ',Hemispheres{iHemi},' ',Metrics{iMetric}])
    end
end



%% Compute sex-difference in the functional metrics - without QC
clc;clear;
load('/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/ForDPABI/Rest_Volume/Pred_Sex.mat');
load('/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/ForDPABI/Rest_Surf_New/Motion_WithID.mat');
InDir = '/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/ForDPABI/Rest_Surf_New/ResultsS/';
OutDir = '/mnt/Data4/RfMRILab/Lubin/Project/QC_Frontiers/Statistic/Rest_Surf_NoQC/';

[C,IA,IB] = intersect(ID,Motion_ID);
SubID = C;
Sex = round(Pred_Sex(IA));
Sex(find(Sex==0)) = -1;
Motion = Motion_FD(IB);
Site = cellfun(@(ID) str2num(ID(5)),SubID);


Metrics  = {'ALFF_FunSurfWC','fALFF_FunSurfWC','ReHo_FunSurfWCF','DegreeCentrality_FunSurfWCF'};
Prefixs = {'szALFF_','szfALFF_','szReHo_','szDegreeCentrality_Bilateral_PositiveWeightedSumBrain_'};
Hemispheres = {'FunSurfLH','FunSurfRH','FunVolu'};
Masks = {'/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V5.1_201201/DPABISurf/SurfTemplates/fsaverage5_lh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V5.1_201201/DPABISurf/SurfTemplates/fsaverage5_rh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/Reslice_freesurfer_subcortical_mask_1mm.nii'};


SiteCov=[];
SiteIndex = unique(Site);
for i=1:length(SiteIndex)-1
    SiteCov(:,i) = Site==SiteIndex(i);
end

BadImageList = {
    'none'
    };
QCFlag = ones(length(SubID),1);
[C,IA,IB] = intersect(SubID,BadImageList);
QCFlag(IA) = 0;

for iHemi = 1:length(Hemispheres)
    for iMetric = 1:length(Metrics)
        SubIndex = find(QCFlag);
        if iHemi < 3
            DependentVolume = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric},filesep,...
                Prefixs{iMetric},ID,'.func.gii'],SubID(SubIndex),'UniformOutput',false);
        else
            DependentVolume = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,regexprep(Metrics{iMetric},'Surf','Volu'),...
                filesep,regexprep(Prefixs{iMetric},'_Bilateral_','_'),ID,'.nii'],SubID(SubIndex),'UniformOutput',false);
        end
        Predictor = [ones(length(SubIndex),1),Sex(SubIndex),Motion(SubIndex),SiteCov(SubIndex,:)];
        OutputDir = [OutDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric}];
        mkdir(OutputDir)
        Contrast = zeros(1,size(Predictor,2));
        Contrast(2) = 1;
        y_GroupAnalysis_Image(DependentVolume,Predictor,[OutputDir,filesep,'T2'],Masks{iHemi},[],Contrast,'T');
        disp(['Have down ',Hemispheres{iHemi},' ',Metrics{iMetric}])
    end
end




