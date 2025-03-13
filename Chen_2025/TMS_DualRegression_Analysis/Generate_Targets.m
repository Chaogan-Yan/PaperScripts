%% 1. [Dual Regression] Generate DR-based TMS targets
%  Generate TMS targets using multiple parameter combinations suggested by reviewers
%  - DLPFC seeds: 
%       sgACC-FC group difference map-based cluster (FDR p=0.1) at MFG
%       Mean MDD sgACC-FC map
%       Mean HC sgACC-FC map
%  - Masks:
%       BA46
%       MFG from AAL3
%       Intersection between MFG and BA46+BA9 (i.e. 'MFG-n-BA46+BA9') 

clc;clear;
%%%% Set parameters %%%%
Approach = 'DR';
SaveImage = 1; % Do not write brain images to save time
nNeighbor = 18; % Criteria for defining a cluster
Verbose = 0; % Mute
PeakMode = 'Size'; % Largest cluster
ParameterDir = '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/';


% Dataset info
Datasets = {'SID','CUD'};
InDirs = {'/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/TMS_WangHuaNing/ForDPABI_Baseline_New';...
    '/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/CUD_TMS/ForDPABI/'};
OutDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Generated_Targets/SID',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Generated_Targets/CUD'};
Phenotypes = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Phenotype/Phenotype_RevisedTargets.mat',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Phenotype/Phenotype_CUD.mat'};
HeadMotinoSubjects = {{'Sub_022'},{'sub-025'}};

% Seedmaps
DLPFCSeeds = {'T2-FDR10','Mean-MDD','Mean-HC'};
SeedDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/T2_sgACC-FC_Fox2012/Reslice_3mm_Target1_FDR10.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_MDD.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_HC.nii'};
    
% DLPFC masks for defining searching scope
DLPFCMasks = {'BA46','MFG','MFG-n-BA46+BA9'};
MaskDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46/Reslice_3mm_BA46_mask.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/AAL3_MFG/Reslice_MFG_mask.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46+9_n_MFG/BA46+9_n_MFG.nii'};

%%%% Generate targets %%%%
nData = length(Datasets);
nSeed = length(DLPFCSeeds);
nMask = length(DLPFCMasks);
for iData = 1:length(Datasets)
    % Load phenotype
    load(Phenotypes{iData});
    
    % Load head motion
    MotionInfo = readtable([InDirs{iData},filesep,'RealignParameter/HeadMotion.tsv'],'FileType','text');
    ExcludeSub = HeadMotinoSubjects{iData}; % max headmotion exclude 3mm or 3degree
    MotionFD = zeros(length(SubID),1);
    for iSub = 1:length(SubID)
        if any(contains(ExcludeSub,SubID{iSub}))
            MotionFD(iSub,1) = 3; % indicate 3mm or 3degree
        else
            Index = find(contains(MotionInfo{:,1},SubID{iSub}));
            MotionFD(iSub,1) = MotionInfo{Index,21}; % line 21 is mean FD-Jenkinson
        end
    end
    
    % Exclude subjects with extremly large head motion
    Index = find(MotionFD~=3);
    SubList_TMS = SubID(Index);
    HAMD_TMS = HAMD(Index,:);
    ActualTarget_TMS = ActualTarget(Index,:);
    Age_TMS = Age(Index);
    Sex_TMS = Sex(Index);
    Motion_TMS = MotionFD(Index);
    
    HAMD_Reducing_Rate = zeros(length(SubList_TMS),1);
    for iSub = 1:length(SubList_TMS)
        HAMD_Reducing_Rate(iSub,1) = (HAMD_TMS(iSub,1) -HAMD_TMS(iSub,2))/HAMD_TMS(iSub,1);
    end
    
    % Load images
    AllData = cell(length(SubList_TMS),1);
    for iSub = 1:length(SubList_TMS)
        [Data,Header] = y_Read([InDirs{iData},filesep,'FunImgARCWF',filesep,SubList_TMS{iSub},filesep,'Filtered_4DVolume.nii']);
        AllData{iSub} = Data;
        disp(['Loading subject ',num2str(iSub)])
    end
    
    for iSeed = 1:nSeed
        [SeedData,~] = y_Read(SeedDirs{iSeed});
        for iMask = 1:nMask
            [MaskData,~] = y_Read(MaskDirs{iMask});
            OutputDir = [OutDirs{iData},filesep,'DLPFCSEED_',DLPFCSeeds{iSeed},'_MASK_',DLPFCMasks{iMask}];
            mkdir(OutputDir);
    
            Optimized_Targets = zeros(length(SubList_TMS),3);
            Offset_Distances = zeros(length(SubList_TMS),1);
            
            % Retrieve the best-performing cluster threshold 
            Parameters = load([ParameterDir,filesep,Datasets{iData},filesep,...
                'GSR_0',filesep,'PEAKMODE_Size',filesep,'MASK_',DLPFCMasks{iMask},filesep,'DLPFCSEED_',DLPFCSeeds{iSeed},filesep,'Summary.mat']);
            [~, Threshold] = min(Parameters.ClinicalSignificance_r);
            
            % Compute targets
            for iSub = 1:length(SubList_TMS)
                [~, ~, Optimized_Targets(iSub,:)] = TMS_Targeting_DR_V4(squeeze(AllData{iSub}), SeedData, [OutputDir,filesep,SubList_TMS{iSub}],MaskData,PeakMode,Threshold,nNeighbor,SaveImage,Header,Verbose);
                Offset_Distances(iSub,1) = pdist([squeeze(Optimized_Targets(iSub,:));ActualTarget_TMS(iSub,:)]);
            end
            
            Regressors = [Offset_Distances,Age_TMS,Sex_TMS,Motion_TMS,ones(length(Age_TMS),1)];
            Contrast = zeros(1,size(Regressors,2));
            Contrast(1) = 1;
            [b,r,SSE,SSR, t, TF_ForContrast, Cohen_f2] = y_regress_ss(HAMD_Reducing_Rate,Regressors,Contrast,'T');
            Df_E = size(Regressors,1) - size(Contrast,2);
            ClinicalSignificance_r = TF_ForContrast./(sqrt(Df_E+TF_ForContrast.*TF_ForContrast));
            ClinicalSignificance_p = 1-tcdf(abs(TF_ForContrast),Df_E);
            
            % Regress out covariates for illustration
            [b,residual,SSE,SSR] = y_regress_ss(HAMD_Reducing_Rate,Regressors(:,2:end));
            HAMD_Reducing_Rate_Regressed=b(end).*ones(size(HAMD_Reducing_Rate))+residual;
            
            save([OutputDir,filesep,'Targets_and_Clinical_Significance.mat'],...
                'SubList_TMS','ClinicalSignificance_r','ClinicalSignificance_p','Optimized_Targets','Offset_Distances',...
                'HAMD_Reducing_Rate','HAMD_Reducing_Rate_Regressed','Threshold');
            
%             % Plot figure
%             figure;
%             set(gcf, "position", [150, 80, 800, 600]);
%             plot(Offset_Distances, HAMD_Reducing_Rate_Regressed, "o");
%             lsline;
%             title({['DLPFCSEED_',DLPFCSeeds{iSeed},'_MASK_',DLPFCMasks{iMask}],...
%                 ['Correlation is ', num2str(ClinicalSignificance_r),' p is ', num2str(ClinicalSignificance_p)]});
%             ylabel("HAMD (%)");
%             xlabel("Euclidean Distance");
%             set(gca, "fontsize", 18);
        end
    end
end


%% 2. [Dual Regression] Generate DR-based TMS targets - Using Alternative Group Difference Map-based Clusters
%  Generate TMS targets using multiple parameter combinations suggested by reviewers
%  - DLPFC seeds: 
%       sgACC-FC group difference map-based cluster (FDR p=0.1) at SFG
%       sgACC-FC group difference map-based cluster (FDR p=0.1) at SFG and MFG
%  - Masks:
%       MFG from AAL3

clc;clear;
%%%% Set parameters %%%%
Approach = 'DR';
SaveImage = 1; % Do not write brain images to save time
nNeighbor = 18; % Criteria for defining a cluster
Verbose = 0; % Mute
PeakMode = 'Size'; % Largest cluster
ParameterDir = '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/';


% Dataset info
Datasets = {'SID','CUD'};
InDirs = {'/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/TMS_WangHuaNing/ForDPABI_Baseline_New';...
    '/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/CUD_TMS/ForDPABI/'};
OutDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Generated_Targets/SID',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Generated_Targets/CUD'};
Phenotypes = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Phenotype/Phenotype_RevisedTargets.mat',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Phenotype/Phenotype_CUD.mat'};
HeadMotinoSubjects = {{'Sub_022'},{'sub-025'}};

% Seedmaps
DLPFCSeeds = {'T2-SFG-FDR10', 'T2-SFG+MFG-FDR10'};
SeedDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/T2_sgACC-FC_Fox2012/Reslice_3mm_Target_SFG_FDR10.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/T2_sgACC-FC_Fox2012/Reslice_3mm_Target_SFG+MFG_FDR10.nii'};
    
% DLPFC masks for defining searching scope
DLPFCMasks = {'MFG'};
MaskDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/AAL3_MFG/Reslice_MFG_mask.nii'};

%%%% Generate targets %%%%
nData = length(Datasets);
nSeed = length(DLPFCSeeds);
nMask = length(DLPFCMasks);
for iData = 1:length(Datasets)
    % Load phenotype
    load(Phenotypes{iData});
    
    % Load head motion
    MotionInfo = readtable([InDirs{iData},filesep,'RealignParameter/HeadMotion.tsv'],'FileType','text');
    ExcludeSub = HeadMotinoSubjects{iData}; % max headmotion exclude 3mm or 3degree
    MotionFD = zeros(length(SubID),1);
    for iSub = 1:length(SubID)
        if any(contains(ExcludeSub,SubID{iSub}))
            MotionFD(iSub,1) = 3; % indicate 3mm or 3degree
        else
            Index = find(contains(MotionInfo{:,1},SubID{iSub}));
            MotionFD(iSub,1) = MotionInfo{Index,21}; % line 21 is mean FD-Jenkinson
        end
    end
    
    % Exclude subjects with extremly large head motion
    Index = find(MotionFD~=3);
    SubList_TMS = SubID(Index);
    HAMD_TMS = HAMD(Index,:);
    ActualTarget_TMS = ActualTarget(Index,:);
    Age_TMS = Age(Index);
    Sex_TMS = Sex(Index);
    Motion_TMS = MotionFD(Index);
    
    HAMD_Reducing_Rate = zeros(length(SubList_TMS),1);
    for iSub = 1:length(SubList_TMS)
        HAMD_Reducing_Rate(iSub,1) = (HAMD_TMS(iSub,1) -HAMD_TMS(iSub,2))/HAMD_TMS(iSub,1);
    end
    
    % Load images
    AllData = cell(length(SubList_TMS),1);
    for iSub = 1:length(SubList_TMS)
        [Data,Header] = y_Read([InDirs{iData},filesep,'FunImgARCWF',filesep,SubList_TMS{iSub},filesep,'Filtered_4DVolume.nii']);
        AllData{iSub} = Data;
        disp(['Loading subject ',num2str(iSub)])
    end
    
    for iSeed = 1:nSeed
        [SeedData,~] = y_Read(SeedDirs{iSeed});
        for iMask = 1:nMask
            [MaskData,~] = y_Read(MaskDirs{iMask});
            OutputDir = [OutDirs{iData},filesep,'DLPFCSEED_',DLPFCSeeds{iSeed},'_MASK_',DLPFCMasks{iMask}];
            mkdir(OutputDir);
    
            Optimized_Targets = zeros(length(SubList_TMS),3);
            Offset_Distances = zeros(length(SubList_TMS),1);
            
            % Retrieve the best-performing cluster threshold 
            Parameters = load([ParameterDir,filesep,Datasets{iData},filesep,...
                'GSR_0',filesep,'PEAKMODE_Size',filesep,'MASK_',DLPFCMasks{iMask},filesep,'DLPFCSEED_',DLPFCSeeds{iSeed},filesep,'Summary.mat']);
            [~, Threshold] = min(Parameters.ClinicalSignificance_r);
            
            % Compute targets
            for iSub = 1:length(SubList_TMS)
                [~, ~, Optimized_Targets(iSub,:)] = TMS_Targeting_DR_V4(squeeze(AllData{iSub}), SeedData, [OutputDir,filesep,SubList_TMS{iSub}],MaskData,PeakMode,Threshold,nNeighbor,SaveImage,Header,Verbose);
                Offset_Distances(iSub,1) = pdist([squeeze(Optimized_Targets(iSub,:));ActualTarget_TMS(iSub,:)]);
            end
            
            Regressors = [Offset_Distances,Age_TMS,Sex_TMS,Motion_TMS,ones(length(Age_TMS),1)];
            Contrast = zeros(1,size(Regressors,2));
            Contrast(1) = 1;
            [b,r,SSE,SSR, t, TF_ForContrast, Cohen_f2] = y_regress_ss(HAMD_Reducing_Rate,Regressors,Contrast,'T');
            Df_E = size(Regressors,1) - size(Contrast,2);
            ClinicalSignificance_r = TF_ForContrast./(sqrt(Df_E+TF_ForContrast.*TF_ForContrast));
            ClinicalSignificance_p = 1-tcdf(abs(TF_ForContrast),Df_E);
            
            % Regress out covariates for illustration
            [b,residual,SSE,SSR] = y_regress_ss(HAMD_Reducing_Rate,Regressors(:,2:end));
            HAMD_Reducing_Rate_Regressed=b(end).*ones(size(HAMD_Reducing_Rate))+residual;
            
            save([OutputDir,filesep,'Targets_and_Clinical_Significance.mat'],...
                'SubList_TMS','ClinicalSignificance_r','ClinicalSignificance_p','Optimized_Targets','Offset_Distances',...
                'HAMD_Reducing_Rate','HAMD_Reducing_Rate_Regressed','Threshold');
            
%             % Plot figure
%             figure;
%             set(gcf, "position", [150, 80, 800, 600]);
%             plot(Offset_Distances, HAMD_Reducing_Rate_Regressed, "o");
%             lsline;
%             title({['DLPFCSEED_',DLPFCSeeds{iSeed},'_MASK_',DLPFCMasks{iMask}],...
%                 ['Correlation is ', num2str(ClinicalSignificance_r),' p is ', num2str(ClinicalSignificance_p)]});
%             ylabel("HAMD (%)");
%             xlabel("Euclidean Distance");
%             set(gca, "fontsize", 18);
        end
    end
end



%% 3. [Group Target] Calculated Clinical Significance of Group Targets Derived from Group-level Statistical Maps

clc;clear;
%%%% Set parameters %%%%
% Dataset info
Datasets = {'SID','CUD'};
InDirs = {'/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/TMS_WangHuaNing/ForDPABI_Baseline_New';...
    '/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/CUD_TMS/ForDPABI/'};
OutDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Generated_Targets/SID',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Generated_Targets/CUD'};
Phenotypes = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Phenotype/Phenotype_RevisedTargets.mat',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Phenotype/Phenotype_CUD.mat'};
HeadMotinoSubjects = {{'Sub_022'},{'sub-025'}};

% Seedmaps
DLPFCSeeds = {'T2-FDR10','Mean-MDD','Mean-HC'};
SeedDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/T2_sgACC-FC_Fox2012/Reslice_3mm_Target1_FDR10.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_MDD.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_HC.nii'};

MaskDir = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46+9_n_MFG/BA46+9_n_MFG.nii'};

% Group targets
GroupTargets = {[-34,36,40],[-40,50,22],[-42,38,32]};

%%%% Calculate clinical significance %%%%
nData = length(Datasets);
nSeed = length(DLPFCSeeds);
for iData = 1:length(Datasets)
    % Load phenotype
    load(Phenotypes{iData});
    
    % Load head motion
    MotionInfo = readtable([InDirs{iData},filesep,'RealignParameter/HeadMotion.tsv'],'FileType','text');
    ExcludeSub = HeadMotinoSubjects{iData}; % max headmotion exclude 3mm or 3degree
    MotionFD = zeros(length(SubID),1);
    for iSub = 1:length(SubID)
        if any(contains(ExcludeSub,SubID{iSub}))
            MotionFD(iSub,1) = 3; % indicate 3mm or 3degree
        else
            Index = find(contains(MotionInfo{:,1},SubID{iSub}));
            MotionFD(iSub,1) = MotionInfo{Index,21}; % line 21 is mean FD-Jenkinson
        end
    end
    
    % Exclude subjects with extremly large head motion
    Index = find(MotionFD~=3);
    SubList_TMS = SubID(Index);
    HAMD_TMS = HAMD(Index,:);
    ActualTarget_TMS = ActualTarget(Index,:);
    Age_TMS = Age(Index);
    Sex_TMS = Sex(Index);
    Motion_TMS = MotionFD(Index);
    
    for iSeed = 1:nSeed
        [SeedData,Header] = y_Read(SeedDirs{iSeed});
        
%         % Calculate the targets (same for every subject)
%         if iSeed == 1
%             [I,J,K] = ind2sub(size(SeedData),find(SeedData));
%         else % add a mask to unthresholded group masks
%             [MaskData,MaskHeader] = y_Read(MaskDir);
%             SeedData = SeedData.*MaskData;
%             [I,J,K] = ind2sub(size(Data),find(SeedData==max(reshape(SeedData,1,[]))));
%         end
%         Coordinates = cor2mni([I,J,K], Header.mat);
%         Group_Target = mean(Coordinates,1);
        
        % Directly use the peak points from the group-level statistical anayses
        Group_Target = GroupTargets{iSeed};
        
        Offset_Distances = zeros(length(SubList_TMS),1);
        HAMD_Reducing_Rate = zeros(length(SubList_TMS),1);
        for iSub = 1:length(SubList_TMS)
            HAMD_Reducing_Rate(iSub,1) = (HAMD_TMS(iSub,1) -HAMD_TMS(iSub,2))/HAMD_TMS(iSub,1);
            Offset_Distances(iSub,1) = pdist([Group_Target;ActualTarget_TMS(iSub,:)]);
        end
         
        Regressors = [Offset_Distances,Age_TMS,Sex_TMS,Motion_TMS,ones(length(Age_TMS),1)];
        Contrast = zeros(1,size(Regressors,2));
        Contrast(1) = 1;
        [b,r,SSE,SSR, t, TF_ForContrast, Cohen_f2] = y_regress_ss(HAMD_Reducing_Rate,Regressors,Contrast,'T');
        Df_E = size(Regressors,1) - size(Contrast,2);
        ClinicalSignificance_r = TF_ForContrast./(sqrt(Df_E+TF_ForContrast.*TF_ForContrast));
        ClinicalSignificance_p = 1-tcdf(abs(TF_ForContrast),Df_E);
        
        % Regress out covariates for illustration
        [b,residual,SSE,SSR] = y_regress_ss(HAMD_Reducing_Rate,Regressors(:,2:end));
        HAMD_Reducing_Rate_Regressed=b(end).*ones(size(HAMD_Reducing_Rate))+residual;
        
        OutputDir = [OutDirs{iData},filesep,'GroupTarget_DLPFCSEED_',DLPFCSeeds{iSeed}];
        mkdir(OutputDir);
        
        save([OutputDir,filesep,'Targets_and_Clinical_Significance.mat'],...
            'SubList_TMS','ClinicalSignificance_r','ClinicalSignificance_p','Group_Target','Offset_Distances',...
            'HAMD_Reducing_Rate','HAMD_Reducing_Rate_Regressed');
        
        % Plot figure
        figure;
        set(gcf, "position", [150, 80, 800, 600]);
        plot(Offset_Distances, HAMD_Reducing_Rate_Regressed, "o");
        lsline;
        title({['DATASET ',Datasets{iData},' DLPFCSEED_',DLPFCSeeds{iSeed}],...
            ['Correlation is ', num2str(ClinicalSignificance_r),' p is ', num2str(ClinicalSignificance_p)]});
        ylabel("HAMD (%)");
        xlabel("Euclidean Distance");
        set(gca, "fontsize", 18);
    end
end


