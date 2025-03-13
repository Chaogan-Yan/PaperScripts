% To Perform stats on the sgACC and DLPFC's FC maps in the DIRECT data
% Xiao Chen 220809
% chenxiaophd@gmail.com

% to select participants 
% modified 221228 by Xiao Chen

%% initilization
clear;clc;

addpath('/mnt/Data3/RfMRILab/Wangyw/codes/matlab/HarmonizationFunctions');
addpath(genpath('/mnt/Data5/RfMRILab/Chenxiao/CX_Software/ComBatHarmonization-master'));

load /mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/AllPhenotype.mat
load /mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Data_MDD_NC.mat

WorkDir='/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Data';
DataDir = '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/group_diff_target/FC_maps/ResultsS';
OutputDir = '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/group_diff_target';
if ~exist(OutputDir, "dir"); mkdir(OutputDir); end

MeasureStringSuffix='/fsaverage5';
MeasureStringSet = {'FC_SeedSurfLHSurfRHVolu_FunSurfWCF', ...
                    'FC_SeedSurfLHSurfRHVolu_FunSurfWglobalCF'};
MeasureStringVoluSet = {'FC_SeedSurfLHSurfRHVolu_FunVoluWCF', ...
                        'FC_SeedSurfLHSurfRHVolu_FunVoluWglobalCF'};

[DPABIPath, fileN, extn] = fileparts(which('DPABI.m'));

%%%%%% change the following codes to modify analyses %%%%%
% All sgACC and DLPFC masks
% MaskSet = {'Brodmann_Mai_Matajnik_BA25_fsaverage5_lh',  ...
%            'glasser_BA25_lh', 'Schaefer_sgACC_lh',...
%            '5cmMethod_20mm_Grey_Matter_fsaverage5_masked_lh',  ...
%            'BA_9_20mm_Grey_Matter_fsaverage5_masked_lh', ...
%            'BA_46_20mm_Grey_Matter_fsaverage5_masked_lh', ...
%            'F3Beam_20mm_Grey_Matter_fsaverage5_masked_lh', ...
%            'SGC_group_target_Grey_Matter_fsaverage5_masked_lh', ...
%            'Brodmann_Mai_Matajnik_BA25_fsaverage5_rh', ...
%            'glasser_BA25_rh', 'Schaefer_sgACC_rh', ...
%            'Fox_2012_Grey_Matter_fsaverage5_masked_rh', ...
%            'Liston_2014_Grey_Matter_fsaverage5_masked_rh', ...
%            'Brodmann_Mai_Matajnik_BA25_lh', 'Brodmann_Mai_Matajnik_BA25_rh',...
%            'Fox_2012_Grey_Matter', 'glasser_BA25_MNI_lh', ...
%            'glasser_BA25_MNI_rh', ...
%            'Liston_2014_Grey_Matter', ...
%            '5cmMethod_20mm_Grey_Matter', ...
%            'BA_9_20mm_Grey_Matter', ...
%            'BA_46_20mm_Grey_Matter', 'F3Beam_20mm_Grey_Matter', ...
%            'SGC_group_target_Grey_Matter'
%                       };
                  
MaskSet = {'group_diff_target_V2_C1_Grey_Matter_fsaverage5_masked_lh', ...
           'group_diff_target_V2_C2_Grey_Matter_fsaverage5_masked_lh', ...
           'group_diff_target_V2_C1_Grey_Matter', 'group_diff_target_V2_C2_Grey_Matter'
           };
                        
%% select subjects and consutruct matrix
SubID=[MDD.ID;NC.ID];
Dx=[ones(length(MDD.ID),1);-1*ones(length(NC.ID),1)];
Age=[MDD.Age;NC.Age];
Sex=[MDD.Sex;NC.Sex];
Edu=[MDD.Education;NC.Education];

Site=zeros(length(SubID),1);
for i=1:length(SubID)
    k = strfind(SubID{i},'-');
    Site(i)=str2num(SubID{i}(3:k-1));
end

% deal with bad coverage, 15 participants were discarded
[Mask_Thrd, ~, ~, Header] = y_ReadAll( ...
    [WorkDir, '/Masks/GroupMask/GroupMask90Percent.nii']);

S = dir([WorkDir,'/Masks/AutoMasks/']);
C = struct2cell(S);
MaskCell = C(1,3:end)';
WantedMatrix = ones(size(MaskCell,1),1);
% some subjects were preprocessesed, but are not included in the
% demographic info table
for i = 1:numel(MaskCell)
    current_subj = MaskCell{i}(1:12);
    if ~ismember(current_subj, SubID)
        WantedMatrix(i) = 0;
    end
end
WantedIndex = find(WantedMatrix);
MaskCell = MaskCell(WantedIndex);
% sort MaskCell according to SubID
MaskCell_reorder = [];
for i = 1:numel(SubID)
   for j = 1:numel(MaskCell)
      if strcmp(SubID{i},MaskCell{j}(1:12))
          MaskCell_reorder{i,1} = MaskCell{j};
          break
      end
   end
end

MaskAll=MaskCell_reorder;
[MaskAll, Vox, FileCell, Header] =y_ReadAll(MaskAll);
GMask=logical(Mask_Thrd);
CoverageVector=zeros(size(MaskCell_reorder,1),1);
for i=1:numel(MaskCell_reorder)
    Mask=logical(MaskAll(:,:,:,i));
    Inter=(Mask & GMask);   
    CoverageVector(i, 1)=sum(Inter(:))/sum(GMask(:)); 
end
CoverageVectorLess90 = CoverageVector<0.9;
WantedSubMatrix=ones(length(SubID),1);
WantedSubMatrix(find(CoverageVectorLess90))=0;

% dealing with head motion, 82 participants were discarded
% get head motion
HeadMotion = [];
for i=1:length(SubID)
    Temp=load(['/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Data/RealignParameter/', ...
                SubID{i},'/FD_Jenkinson_',SubID{i},'.txt']);
    HeadMotion(i,1)=mean(Temp);
end

MotionGreater02 = HeadMotion > 0.2;
WantedSubMatrix(find(MotionGreater02))=0;

% rule out NaN, 20 subjects were discarded
matrix4NaN = [Age,Sex,Edu,HeadMotion];
HasNaN = isnan(sum(matrix4NaN,2));
WantedSubMatrix(HasNaN) = 0;

%select participants
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
HeadMotion=HeadMotion(WantedSubIndex,:);

save('/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/Sub_Info_1574_1308.mat', ...
    'SubID', 'Dx', 'Age', 'Sex', 'Edu', 'Site', 'HeadMotion', '-mat');

%%%%%%% constructing regression matrix %%%%%%%
% use dummy code to regress out site effect
SiteCov=[];
SiteIndex = unique(Site);
for i=1:length(SiteIndex)-1
    SiteCov(:,i) = Site==SiteIndex(i);
end
AllCov = [ones(length(SubID),1),Dx,Age,Sex,Edu,HeadMotion, SiteCov];


%Centering: Let the first column (constant) have the mean effect.
AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
Contrast=zeros(1,size(AllCov,2));
Contrast(2)=1;
DOF = size(AllCov,1) - size(AllCov,2);


% %% deal with the EMBARC dataset
% sub_info = readtable(['/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses/', ...
%                       'FCMaps_EMBARC/EMBARC_demographic_info.csv']);
% mask1 = sub_info.dx == 1;
% mask2 = sub_info.session == 1;
% mask = mask1.*mask2;
% MDDID = sub_info.SubjectID_BIDS(find(mask));  
% MDDSite = sub_info.site(find(mask));
% 
% mask1 = sub_info.dx == 2;
% mask2 = sub_info.session == 1;
% mask = mask1.*mask2;
% HCID = sub_info.SubjectID_BIDS(find(mask)); 
% HCSite = sub_info.site(find(mask)); 
% 
% SubID=[MDDID;HCID];
% Dx=[ones(length(MDDID),1);-1*ones(length(HCID),1)];
% Site = [MDDSite; HCSite];
% 
% % use dummy code to regress out site effect
% SiteCov=[];
% SiteIndex = unique(Site);
% for i=1:length(SiteIndex)-1
%     SiteCov(:,i) = Site==SiteIndex(i);
% end
% AllCov = [ones(length(SubID),1), Dx,SiteCov];
% 
% %Centering: Let the first column (constant) have the mean effect.
% AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
% Contrast=zeros(1,size(AllCov,2));
% Contrast(2)=1;
% DOF = size(AllCov,1) - size(AllCov,2);
% Table.SubjID=SubID;

%% open par pool
poolobj = gcp('nocreate');
delete(poolobj);
parpool(12);

%% do stats on the surface
for iGSR = 1:length(MeasureStringSet)
    MeasureString = MeasureStringSet{iGSR};
    parfor iMask = 1:length(MaskSet)
        ROI = MaskSet{iMask};

        FileListLH=[];
        FileListRH=[];
        for iSub=1:length(Table.SubjID)
            FileListLH{iSub,1}=sprintf('%s/FunSurfLH/%s_%s/szFC_%s.func.gii', ...
                DataDir,MeasureString,ROI,Table.SubjID{iSub});
            FileListRH{iSub,1}=sprintf('%s/FunSurfRH/%s_%s/szFC_%s.func.gii', ...
                DataDir,MeasureString,ROI,Table.SubjID{iSub});
        end


        MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage5_lh_cortex.label.gii');
        mkdir([OutputDir,'/FuncSurfLH/',ROI])
        OutputName=[OutputDir,'/FuncSurfLH/',ROI,'/Dx_', MeasureString,'.gii'];
        y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);

        MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage5_rh_cortex.label.gii');
        mkdir([OutputDir,'/FuncSurfRH/',ROI])
        OutputName=[OutputDir,'/FuncSurfRH/',ROI,'/Dx_', MeasureString,'.gii'];
        y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    end
end


%% do stats in volume space
MaskFile='/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Data/AllResampled_BrainMask_05_91x109x91.nii';
for iGSR = 1:length(MeasureStringVoluSet)
    MeasureString = MeasureStringVoluSet{iGSR};
    parfor iMask = 1:length(MaskSet)
        ROI = MaskSet{iMask};

        FileList=[];

        for iSub=1:length(Table.SubjID)
            FileList{iSub,1}=sprintf('%s/FunVolu/%s_%s/szFC_%s.nii',DataDir, ...
                MeasureString,ROI,Table.SubjID{iSub}); 
        end

        mkdir([OutputDir,'/FunVolu/',ROI])
        OutputName=[OutputDir,'/FunVolu/',ROI,'/Dx_', MeasureString,'.nii'];
        y_GroupAnalysis_Image(FileList,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% comBat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for comBat
load '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/Sub_Info_1574_1308.mat'

SiteIndex = unique(Site);
Batch=[];
for i=1:length(SiteIndex)
    Batch(Site==SiteIndex(i))=i;
end

AllCov_comBat = [Dx,Age,Sex,Edu];

ComBatParam.batch = Batch;
ComBatParam.mod = AllCov_comBat;
ComBatParam.isparametric = 1;

% %% deal with the EMBARC dataset
% sub_info = readtable(['/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses/', ...
%                       'FCMaps_EMBARC/EMBARC_demographic_info.csv']);
% mask1 = sub_info.dx == 1;
% mask2 = sub_info.session == 1;
% mask = mask1.*mask2;
% MDDID = sub_info.SubjectID_BIDS(find(mask));  
% MDDSite = sub_info.site(find(mask));
% 
% mask1 = sub_info.dx == 2;
% mask2 = sub_info.session == 1;
% mask = mask1.*mask2;
% HCID = sub_info.SubjectID_BIDS(find(mask)); 
% HCSite = sub_info.site(find(mask)); 
% 
% SubID=[MDDID;HCID];
% Dx=[ones(length(MDDID),1);-1*ones(length(HCID),1)];
% Site = [MDDSite; HCSite];
% 
% SiteIndex = unique(Site);
% Batch=[];
% for i=1:length(SiteIndex)
%     Batch(Site==SiteIndex(i))=i;
% end
% 
% AllCov_comBat = Dx;
% 
% ComBatParam.batch = Batch;
% ComBatParam.mod = AllCov_comBat;
% ComBatParam.isparametric = 1;

%% open par pool
poolobj = gcp('nocreate');
delete(poolobj);
parpool(4);
%% do comBat on surface
for iGSR = 1:length(MeasureStringSet)
    MeasureString = MeasureStringSet{iGSR};
    parfor iMask = 1:length(MaskSet)
        ROI = MaskSet{iMask};

        FileListLH=[];
        FileListRH=[];
        for iSub=1:length(SubID)
            FileListLH{iSub,1}=sprintf('%s/FunSurfLH/%s_%s/szFC_%s.func.gii', ...
                DataDir,MeasureString,ROI,SubID{iSub});
            FileListRH{iSub,1}=sprintf('%s/FunSurfRH/%s_%s/szFC_%s.func.gii', ...
                DataDir,MeasureString,ROI,SubID{iSub});
        end

        MaskFile_LH =fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage5_lh_cortex.label.gii');
        MaskFile_RH =fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage5_rh_cortex.label.gii');
        OutputName = [OutputDir, '/comBat'];
        if ~exist(OutputName, "dir"); mkdir(OutputName); end
        yw_Harmonization_Surf(FileListLH,FileListRH,MaskFile_LH,MaskFile_RH,1,ComBatParam,OutputName,[]);
    end
end

%% do combat in the volume space
OutputName = [OutputDir,'/comBat/Volu/'];
if ~exist(OutputName, "dir"); mkdir(OutputName); end
MaskFile='/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Data/AllResampled_BrainMask_05_91x109x91.nii';
MaskData = y_Read(MaskFile);
for iGSR = 1:length(MeasureStringVoluSet)
    MeasureString = MeasureStringVoluSet{iGSR};
    parfor iMask = 1:length(MaskSet)
        ROI = MaskSet{iMask};

        FileList=[];

        for iSub=1:length(SubID)
            FileList{iSub,1}=sprintf('%s/FunVolu/%s_%s/szFC_%s.nii',DataDir, ...
                MeasureString,ROI,SubID{iSub}); 
        end

       yw_Harmonization(FileList,MaskData,1,ComBatParam,OutputName,[]);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%% do stats on images after comBat %%%%%%%%%%%%%%%%%
%% select subjects and consutruct matrix for two-sample t test
AllCov = [];
load '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/Sub_Info_1574_1308.mat';

AllCov = [ones(length(SubID),1),Dx,Age,Sex,Edu,HeadMotion];%site was not included after comBat

%Centering: Let the first column (constant) have the mean effect.
AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
Contrast=zeros(1,size(AllCov,2));
Contrast(2)=1;
DOF = size(AllCov,1) - size(AllCov,2);

%% open par pool
poolobj = gcp('nocreate');
delete(poolobj);
parpool(4);
%% do stats on the surface
for iGSR = 1:length(MeasureStringSet)
    MeasureString = MeasureStringSet{iGSR};
    parfor iMask = 1:length(MaskSet)
        ROI = MaskSet{iMask};


        FileListLH=[];
        FileListRH=[];
        for iSub=1:length(SubID)
            FileListLH{iSub,1}=sprintf('%s/LH/%s_%s/szFC_%s.func.gii', ...
                [OutputDir,'/comBat'],MeasureString,ROI,SubID{iSub});
            FileListRH{iSub,1}=sprintf('%s/RH/%s_%s/szFC_%s.func.gii', ...
                [OutputDir,'/comBat'],MeasureString,ROI,SubID{iSub});
        end


        MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage5_lh_cortex.label.gii');
        mkdir([OutputDir,'/stats_comBat/FuncSurfLH/',ROI])
        OutputName=[OutputDir,'/stats_comBat/FuncSurfLH/',ROI,'/Dx_', MeasureString,'.gii'];
        y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);

        MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage5_rh_cortex.label.gii');
        mkdir([OutputDir,'/stats_comBat/FuncSurfRH/',ROI])
        OutputName=[OutputDir,'/stats_comBat/FuncSurfRH/',ROI,'/Dx_', MeasureString,'.gii'];
        y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    end
end

%% do stats in volume space
MaskFile='/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Data/AllResampled_BrainMask_05_91x109x91.nii';
for iGSR = 1:length(MeasureStringVoluSet)
    MeasureString = MeasureStringVoluSet{iGSR};
    parfor iMask = 1:length(MaskSet)
        ROI = MaskSet{iMask};

        FileList=[];

        for iSub=1:length(SubID)
            FileList{iSub,1}=sprintf('%s/Volu/%s_%s/szFC_%s.nii',[OutputDir,'/comBat'], ...
                MeasureString,ROI,SubID{iSub}); 
        end

        mkdir([OutputDir,'/stats_comBat/FunVolu/',ROI])
        OutputName=[OutputDir,'/stats_comBat/FunVolu/',ROI,'/Dx_', MeasureString,'.nii'];
        y_GroupAnalysis_Image(FileList,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    end
end

%% select subjects and consutruct matrix for one-sample t test
AllCov = [];
load /mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses/Fox_group_target_redo/AllPhenotype.mat
load /mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses/Fox_group_target_redo/Data_MDD_NC.mat

% % MDD
% SubID = MDD.ID;
% Age = MDD.Age;
% Sex = MDD.Sex;
% Edu = MDD.Education;

%NC
SubID = NC.ID;
Age = NC.Age;
Sex = NC.Sex;
Edu = NC.Education;


HeadMotion = [];
for i=1:length(SubID)
    Temp=load(['/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Data/RealignParameter/', ...
                SubID{i},'/FD_Jenkinson_',SubID{i},'.txt']);
    HeadMotion(i,1)=mean(Temp);
end

% codes for site, unnecessary
Site=zeros(length(SubID),1);
for i=1:length(SubID)
    k = strfind(SubID{i},'-');
    Site(i)=str2num(SubID{i}(3:k-1));
end

% use dummy code to regress out site effect
SiteCov=[];
SiteIndex = unique(Site);
for i=1:length(SiteIndex)-1
    SiteCov(:,i) = Site==SiteIndex(i);
end

AllCov = [ones(length(SubID),1),Age,Sex,Edu,HeadMotion];%site was not included after comBat

%Deal with NaN
HasNaN = isnan(sum(AllCov,2));
AllCov(HasNaN,:)=[];
SubID(HasNaN,:)=[];

% IS017-1-0047 bad coverage
tmep_idx = find(strcmp(SubID,'IS017-1-0047'));
AllCov(tmep_idx,:)=[];
SubID(tmep_idx,:)=[];

%Centering: Let the first column (constant) have the mean effect.
AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
Contrast=zeros(1,size(AllCov,2));
Contrast(1)=1;
DOF = size(AllCov,1) - size(AllCov,2);

%% do stats on the surface
iGSR = 2;
MeasureString = MeasureStringSet{iGSR};
for iMask = 1:length(MaskSet)
    ROI = MaskSet{iMask};

    FileListLH=[];
    FileListRH=[];
    for iSub=1:length(SubID)
        FileListLH{iSub,1}=sprintf('%s/LH/%s_%s/szFC_%s.func.gii', ...
            [OutputDir,'/comBat'],MeasureString,ROI,SubID{iSub});
        FileListRH{iSub,1}=sprintf('%s/RH/%s_%s/szFC_%s.func.gii', ...
            [OutputDir,'/comBat'],MeasureString,ROI,SubID{iSub});
    end

    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage5_lh_cortex.label.gii');
    mkdir([OutputDir,'/stats_comBat_NC/FuncSurfLH/',ROI])
    OutputName=[OutputDir,'/stats_comBat_NC/FuncSurfLH/',ROI,'/Dx_', MeasureString,'.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);

    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage5_rh_cortex.label.gii');
    mkdir([OutputDir,'/stats_comBat_NC/FuncSurfRH/',ROI])
    OutputName=[OutputDir,'/stats_comBat_NC/FuncSurfRH/',ROI,'/Dx_', MeasureString,'.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end

% average method from Weigand et al., 2018
[data,~,~,header] = y_ReadAll(FileListLH);
data_average = mean(data,2);
y_Write(data_average, header, [OutputDir,'/group_mean_sgACC_FCMaps_MDD_lh.gii']);

%% do stats in volume space
MaskFile='/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Data/AllResampled_BrainMask_05_91x109x91.nii';
iGSR = 2;
MeasureString = MeasureStringVoluSet{iGSR};
for iMask = 1:length(MaskSet)
    ROI = MaskSet{iMask};

    FileList=[];

    for iSub=1:length(SubID)
        FileList{iSub,1}=sprintf('%s/Volu/%s_%s/szFC_%s.nii',[OutputDir,'/comBat'], ...
            MeasureString,ROI,SubID{iSub}); 
    end

    mkdir([OutputDir,'/stats_comBat_NC/FunVolu/',ROI])
    OutputName=[OutputDir,'/stats_comBat_NC/FunVolu/',ROI,'/Dx_', MeasureString,'.nii'];
    y_GroupAnalysis_Image(FileList,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end

% 
[data,~,~,header] = y_ReadAll(FileList);
data_average = squeeze(mean(data,4));
y_Write(data_average, header, [OutputDir,'/group_mean_sgACC_FCMaps_HC.nii']);
