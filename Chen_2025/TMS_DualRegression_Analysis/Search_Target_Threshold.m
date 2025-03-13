%% [Searching Parameters] Search best thresholds (number of voxels to keep) for dual-regression-generated maps to determine final targets
clc;clear;

%%%% Set parameters %%%%
CluSizeRange = [5,100];
SaveImage = 0; % Do not write brain images to save time
nNeighbor = 18; % Criteria for defining a cluster
Verbose = 0; % Mute

% Datasets
Datasets = {'SID','CUD'};
InDirs = {'/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/TMS_WangHuaNing/ForDPABI_Baseline_New';...
    '/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/CUD_TMS/ForDPABI/'};
OutDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/SID',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/CUD'};
Phenotypes = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Phenotype/Phenotype_RevisedTargets.mat',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Phenotype/Phenotype_CUD.mat'};
HeadMotinoSubjects = {{'Sub_022'},{'sub-025'}};

% Seedmaps
DLPFCSeeds = {'T2-FDR10','Mean-MDD','Mean-HC'};
SeedDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/T2_sgACC-FC_Fox2012/Reslice_3mm_Target1_FDR10.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_MDD.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_HC.nii'};
    
% DLPFC masks for defining searching scope
DLPFCMasks = {'BA46','MFG','MFG+BA46','MFG-n-BA46+BA9'};
MaskDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46/Reslice_3mm_BA46_mask.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/AAL3_MFG/Reslice_MFG_mask.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46+MFG/BA46+MFG.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46+9_n_MFG/BA46+9_n_MFG.nii'};

% Approaches for defining peak points as targets
PeakModes = {'Size','WeightedSize'}; 

% Do GSR or not 
Suffixes = {'FunImgARCWF','FunImgARglobalCWF'};


%%%% Grid search parameters%%%%
nData = length(Datasets);
nGSR = length(Suffixes);
nPeak = length(PeakModes);
nSeed = length(DLPFCSeeds);
nMask = length(DLPFCMasks);
nStep = CluSizeRange(2)-CluSizeRange(1)+1;
TotalSteps = nData*nSeed*nMask*nPeak*nGSR*nStep;
WaitBar = waitbar(0,'Searching proper cluster size');

for iData = 1:nData
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
    
    for iGSR = 1:nGSR
        % Load images
        AllData = cell(length(SubList_TMS),1);
        for iSub = 1:length(SubList_TMS)
            [Data,Header] = y_Read([InDirs{iData},filesep,Suffixes{iGSR},filesep,SubList_TMS{iSub},filesep,'Filtered_4DVolume.nii']);
            AllData{iSub} = Data;
            disp(['Loading subject ',num2str(iSub)])
        end
        for iPeak = 1:nPeak
            for iSeed = 1:nSeed
                [SeedData,~] = y_Read(SeedDirs{iSeed});
                for iMask = 1:nMask
                    [MaskData,~] = y_Read(MaskDirs{iMask});
                    OutputDir = [OutDirs{iData},filesep,'GSR_',num2str(iGSR-1),filesep,'PEAKMODE_',PeakModes{iPeak},filesep,'MASK_',DLPFCMasks{iMask},filesep,'DLPFCSEED_',DLPFCSeeds{iSeed}];
                    mkdir(OutputDir);
                    Optimized_Targets = zeros(length(SubList_TMS),CluSizeRange(2),3);
                    Offset_Distances = zeros(length(SubList_TMS),CluSizeRange(2));
                    ClinicalSignificance_r = zeros(CluSizeRange(2),1);
                    ClinicalSignificance_p = zeros(CluSizeRange(2),1);
                    for iStep = CluSizeRange(1):CluSizeRange(2)
                        CurrentStep = ...
                            (iData - 1) * (nGSR * nPeak * nSeed * nMask * nStep) + ...
                            (iGSR  - 1) * (nPeak * nSeed * nMask * nStep) + ...
                            (iPeak - 1) * (nSeed * nMask * nStep)          + ...
                            (iSeed - 1) * (nMask * nStep)                  + ...
                            (iMask - 1) * nStep                            + ...
                            (iStep - CluSizeRange(1) +1);
                        Progress = CurrentStep / TotalSteps;
                        waitbar(Progress, WaitBar,  ...
                            sprintf('Data=%d/%d, GSR=%d/%d, Peak=%d/%d, Seed=%d/%d, Mask=%d/%d, Steps=%d/%d', ...
                            iData,nData,...
                            iGSR, nGSR, ...
                            iPeak, nPeak, ...
                            iSeed, nSeed, ...
                            iMask, nMask, ...
                            iStep - 4, nStep));
                        
                        for iSub = 1:length(SubList_TMS)
                            [~, ~, Optimized_Targets(iSub,iStep,:)] = TMS_Targeting_DR_V4(squeeze(AllData{iSub}), SeedData, [OutputDir,filesep,SubList_TMS{iSub}], MaskData,PeakModes{iPeak},iStep,nNeighbor,SaveImage,Header,Verbose);
                            Offset_Distances(iSub,iStep) = pdist([squeeze(Optimized_Targets(iSub,iStep,:))';ActualTarget_TMS(iSub,:)]);
                        end
                        
                        Regressors = [Offset_Distances(:,iStep),Age_TMS,Sex_TMS,Motion_TMS,ones(length(Age_TMS),1)];
                        Contrast = zeros(1,size(Regressors,2));
                        Contrast(1) = 1;
                        [b,r,SSE,SSR, t, TF_ForContrast, Cohen_f2] = y_regress_ss(HAMD_Reducing_Rate,Regressors,Contrast,'T');
                        Df_E = size(Regressors,1) - size(Contrast,2);
                        ClinicalSignificance_r(iStep,1) = TF_ForContrast./(sqrt(Df_E+TF_ForContrast.*TF_ForContrast));
                        ClinicalSignificance_p(iStep,1) = 1-tcdf(abs(TF_ForContrast),Df_E);
                    end
                    save([OutputDir,filesep,'Summary.mat'],...
                        'SubList_TMS','ClinicalSignificance_r','ClinicalSignificance_p','Optimized_Targets','Offset_Distances','HAMD_Reducing_Rate');
                end
            end
        end
    end
end


%% [Searching Parameters] Test Alternative Seeds derived from T2 Group Difference Map (i.e. SFG cluster, SFG+MFG cluster) 
% Suggested by a reviewer
clc;clear;

%%%% Setting parameters %%%%
CluSizeRange = [5,100];
SaveImage = 0; % Do not write brain images to save time
nNeighbor = 18; % Criteria for defining a cluster
Verbose = 0; % Mute

% Datasets
Datasets = {'SID','CUD'};
InDirs = {'/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/TMS_WangHuaNing/ForDPABI_Baseline_New';...
    '/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/CUD_TMS/ForDPABI/'};
OutDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/SID',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/CUD'};
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

% Approaches for defining peak points as targets
PeakModes = {'Size'}; 

% Do GSR or not 
Suffixes = {'FunImgARCWF'};


%%%% Grid search parameters%%%%
nData = length(Datasets);
nGSR = length(Suffixes);
nPeak = length(PeakModes);
nSeed = length(DLPFCSeeds);
nMask = length(DLPFCMasks);
nStep = CluSizeRange(2)-CluSizeRange(1)+1;
TotalSteps = nData*nSeed*nMask*nPeak*nGSR*nStep;
WaitBar = waitbar(0,'Searching proper cluster size');

for iData = 1:nData
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
    
    for iGSR = 1:nGSR
        % Load images
        AllData = cell(length(SubList_TMS),1);
        for iSub = 1:length(SubList_TMS)
            [Data,Header] = y_Read([InDirs{iData},filesep,Suffixes{iGSR},filesep,SubList_TMS{iSub},filesep,'Filtered_4DVolume.nii']);
            AllData{iSub} = Data;
            disp(['Loading subject ',num2str(iSub)])
        end
        for iPeak = 1:nPeak
            for iSeed = 1:nSeed
                [SeedData,~] = y_Read(SeedDirs{iSeed});
                for iMask = 1:nMask
                    [MaskData,~] = y_Read(MaskDirs{iMask});
                    OutputDir = [OutDirs{iData},filesep,'GSR_',num2str(iGSR-1),filesep,'PEAKMODE_',PeakModes{iPeak},filesep,'MASK_',DLPFCMasks{iMask},filesep,'DLPFCSEED_',DLPFCSeeds{iSeed}];
                    mkdir(OutputDir);
                    Optimized_Targets = zeros(length(SubList_TMS),CluSizeRange(2),3);
                    Offset_Distances = zeros(length(SubList_TMS),CluSizeRange(2));
                    ClinicalSignificance_r = zeros(CluSizeRange(2),1);
                    ClinicalSignificance_p = zeros(CluSizeRange(2),1);
                    for iStep = CluSizeRange(1):CluSizeRange(2)
                        CurrentStep = ...
                            (iData - 1) * (nGSR * nPeak * nSeed * nMask * nStep) + ...
                            (iGSR  - 1) * (nPeak * nSeed * nMask * nStep) + ...
                            (iPeak - 1) * (nSeed * nMask * nStep)          + ...
                            (iSeed - 1) * (nMask * nStep)                  + ...
                            (iMask - 1) * nStep                            + ...
                            (iStep - CluSizeRange(1) +1);
                        Progress = CurrentStep / TotalSteps;
                        waitbar(Progress, WaitBar,  ...
                            sprintf('Data=%d/%d, GSR=%d/%d, Peak=%d/%d, Seed=%d/%d, Mask=%d/%d, Steps=%d/%d', ...
                            iData,nData,...
                            iGSR, nGSR, ...
                            iPeak, nPeak, ...
                            iSeed, nSeed, ...
                            iMask, nMask, ...
                            iStep - 4, nStep));
                        
                        for iSub = 1:length(SubList_TMS)
                            [~, ~, Optimized_Targets(iSub,iStep,:)] = TMS_Targeting_DR_V4(squeeze(AllData{iSub}), SeedData, [OutputDir,filesep,SubList_TMS{iSub}], MaskData,PeakModes{iPeak},iStep,nNeighbor,SaveImage,Header,Verbose);
                            Offset_Distances(iSub,iStep) = pdist([squeeze(Optimized_Targets(iSub,iStep,:))';ActualTarget_TMS(iSub,:)]);
                        end
                        
                        Regressors = [Offset_Distances(:,iStep),Age_TMS,Sex_TMS,Motion_TMS,ones(length(Age_TMS),1)];
                        Contrast = zeros(1,size(Regressors,2));
                        Contrast(1) = 1;
                        [b,r,SSE,SSR, t, TF_ForContrast, Cohen_f2] = y_regress_ss(HAMD_Reducing_Rate,Regressors,Contrast,'T');
                        Df_E = size(Regressors,1) - size(Contrast,2);
                        ClinicalSignificance_r(iStep,1) = TF_ForContrast./(sqrt(Df_E+TF_ForContrast.*TF_ForContrast));
                        ClinicalSignificance_p(iStep,1) = 1-tcdf(abs(TF_ForContrast),Df_E);
                    end
                    save([OutputDir,filesep,'Summary.mat'],...
                        'SubList_TMS','ClinicalSignificance_r','ClinicalSignificance_p','Optimized_Targets','Offset_Distances','HAMD_Reducing_Rate');
                end
            end
        end
    end
end



%% [Compare Paremeters] 1. Compare PeakMode - Peak Mode "Size" is better
clc;clear;close all;

OutDir = '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Code/ForManuscript_TMS_SciBul/Figure/Search_Proper_Target_Threshold/1_Compare_PeakMode';
mkdir(OutDir);

% Datasets
Datasets = {'SID','CUD'};
InDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/SID',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/CUD'};

% Seedmaps
DLPFCSeeds = {'T2-FDR10','Mean-MDD','Mean-HC'};
SeedDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/T2_sgACC-FC_Fox2012/Reslice_3mm_Target1_FDR10.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_MDD.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_HC.nii'};
    
% DLPFC masks for defining searching scope
DLPFCMasks = {'BA46','MFG','MFG+BA46','MFG-n-BA46+BA9'};
MaskDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46/Reslice_3mm_BA46_mask.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/AAL3_MFG/Reslice_MFG_mask.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46+MFG/BA46+MFG.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46+9_n_MFG/BA46+9_n_MFG.nii'};

% Approaches for defining peak points as targets
PeakModes = {'Size','WeightedSize'}; 

% Do GSR or not 
Suffixes = {'FunImgARCWF','FunImgARglobalCWF'};

CluSizeRange = [5,100];

nData = length(Datasets);
nGSR = length(Suffixes);
nPeak = length(PeakModes);
nSeed = length(DLPFCSeeds);
nMask = length(DLPFCMasks);
nStep = CluSizeRange(2)-CluSizeRange(1)+1;

for iData = 1:nData
    for iSeed = 1:nSeed
        figure('Units','normalized','Position',[0.1,0.1,0.8,0.6],'Color','w');
        tiledlayout(nGSR, nMask, 'TileSpacing', 'compact', 'Padding', 'compact');
        sgtitle(['DATASET-',Datasets{iData},'  DLPFCSEED-',DLPFCSeeds{iSeed}], 'FontWeight', 'bold');
        for iGSR = 1:nGSR
            for iMask = 1:nMask
                nexttile
                for iPeak = 1:nPeak
                    load([InDirs{iData},filesep,'GSR_',num2str(iGSR-1),filesep,'PEAKMODE_',PeakModes{iPeak},filesep,'MASK_',DLPFCMasks{iMask},filesep,'DLPFCSEED_',DLPFCSeeds{iSeed},filesep,'Summary.mat']);
                    
                    plot(1:CluSizeRange(2), ClinicalSignificance_r, 'LineWidth', 1.5);
                    [~,MinIndex] = min(ClinicalSignificance_r);
                    xline(MinIndex,'--','LineWidth', 1.5, 'Color', [0.5,0,0]);
                    grid on;
                    box on;
                    hold on;
                    
                    xlim([1,CluSizeRange(2)]);
                    ylim([-0.6,0]);
                    ylabel({'Clinical Significance', '(Pearson''s Correlation between Targeting', 'Offset and Clinical Improvement)'}, 'FontSize', 9);
                    xlabel({'Threshold for Individualized Targeting Map', '(Number of Voxel)'}, 'FontSize', 9);
                    title(regexprep(['GSR_',num2str(iGSR-1),filesep,'MASK_',DLPFCMasks{iMask}],'_','-'), 'FontSize', 10, 'FontWeight', 'normal');
                    legend('Size','WeightedSize','Location','best');
                end
            end
        end
        
        fig = gcf;
        exportgraphics(fig,[OutDir,filesep,'DATASET-',Datasets{iData},'  DLPFCSEED-',DLPFCSeeds{iSeed},'.jpg'],'Resolution',300);
    end
end


%% [Compare Paremeters] 2. Compare GSR - No GSR is slightly better
% Fixed PeakMode-"Size"
clc;clear;close all;

OutDir = '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Code/ForManuscript_TMS_SciBul/Figure/Search_Proper_Target_Threshold/2_Compare_PeakGSR';
mkdir(OutDir);

% Datasets
Datasets = {'SID','CUD'};
InDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/SID',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/CUD'};

% Seedmaps
DLPFCSeeds = {'T2-FDR10','Mean-MDD','Mean-HC'};
SeedDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/T2_sgACC-FC_Fox2012/Reslice_3mm_Target1_FDR10.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_MDD.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_HC.nii'};
    
% DLPFC masks for defining searching scope
DLPFCMasks = {'BA46','MFG','MFG+BA46','MFG-n-BA46+BA9'};
MaskDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46/Reslice_3mm_BA46_mask.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/AAL3_MFG/Reslice_MFG_mask.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46+MFG/BA46+MFG.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46+9_n_MFG/BA46+9_n_MFG.nii'};

% Approaches for defining peak points as targets
PeakModes = {'Size'}; 

% Do GSR or not 
Suffixes = {'FunImgARCWF','FunImgARglobalCWF'};

CluSizeRange = [5,100];

nData = length(Datasets);
nGSR = length(Suffixes);
nPeak = length(PeakModes);
nSeed = length(DLPFCSeeds);
nMask = length(DLPFCMasks);
nStep = CluSizeRange(2)-CluSizeRange(1)+1;

for iData = 1:nData
    for iPeak = 1:nPeak
        figure('Units','normalized','Position',[0.1,0.1,0.9,0.9],'Color','w');
        tiledlayout(nSeed, nMask, 'TileSpacing', 'compact', 'Padding', 'compact');
        sgtitle(['DATASET-',Datasets{iData}], 'FontWeight', 'bold');
        for iSeed = 1:nSeed
            for iMask = 1:nMask
                nexttile
                for iGSR = 1:nGSR
                    load([InDirs{iData},filesep,'GSR_',num2str(iGSR-1),filesep,'PEAKMODE_',PeakModes{iPeak},filesep,'MASK_',DLPFCMasks{iMask},filesep,'DLPFCSEED_',DLPFCSeeds{iSeed},filesep,'Summary.mat']);
                    
                    plot(1:CluSizeRange(2), ClinicalSignificance_r, 'LineWidth', 1.5);
                    [~,MinIndex] = min(ClinicalSignificance_r);
                    xline(MinIndex,'--','LineWidth', 1.5, 'Color', [0.5,0,0]);
                    grid on;
                    box on;
                    hold on;
                    
                    xlim([1,CluSizeRange(2)]);
                    ylim([-0.6,0]);
                    ylabel({'Clinical Significance', '(Pearson''s Correlation between Targeting', 'Offset and Clinical Improvement)'}, 'FontSize', 9);
                    xlabel({'Threshold for Individualized Targeting Map', '(Number of Voxel)'}, 'FontSize', 9);
                    title(regexprep(['DLPFCSEED_',DLPFCSeeds{iSeed},filesep,'MASK_',DLPFCMasks{iMask}],'_','-'), 'FontSize', 10, 'FontWeight', 'normal');
                    legend('No GSR','GSR','Location','best');
                end
            end
        end
        
        fig = gcf;
        exportgraphics(fig,[OutDir,filesep,'DATASET-',Datasets{iData},'.jpg'],'Resolution',300);
    end
end


%% [Compare Paremeters] 3. Compare Masks - BA46 fits T2 seeds better, MFG-n-BA46+BA9 fits mean-FC maps better
% Fixed PeakMode-"Size" and GSR-"No GSR"
clc;clear;close all;

OutDir = '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Code/ForManuscript_TMS_SciBul/Figure/Search_Proper_Target_Threshold/3_Compare_Masks';
mkdir(OutDir);

% Datasets
Datasets = {'SID','CUD'};
InDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/SID',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/CUD'};

% Seedmaps
DLPFCSeeds = {'T2-FDR10','Mean-MDD','Mean-HC'};
SeedDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/T2_sgACC-FC_Fox2012/Reslice_3mm_Target1_FDR10.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_MDD.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/Mean_sgACC-FC_Fox2012/Fliped_Reslice_3mm_group_mean_sgACC_FCMaps_HC.nii'};
    
% DLPFC masks for defining searching scope
DLPFCMasks = {'BA46','MFG','MFG-n-BA46+BA9'};  % 'MFG+BA46' is not required by the reviewer and is not promising
% DLPFCMasks = {'BA46','MFG','MFG+BA46','MFG-n-BA46+BA9'};
MaskDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46/Reslice_3mm_BA46_mask.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/AAL3_MFG/Reslice_MFG_mask.nii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/BA46+9_n_MFG/BA46+9_n_MFG.nii'};

% Approaches for defining peak points as targets
PeakModes = {'Size'}; 

% Do GSR or not 
Suffixes = {'FunImgARCWF'};

CluSizeRange = [5,100];

nData = length(Datasets);
nGSR = length(Suffixes);
nPeak = length(PeakModes);
nSeed = length(DLPFCSeeds);
nMask = length(DLPFCMasks);
nStep = CluSizeRange(2)-CluSizeRange(1)+1;

for iPeak = 1:nPeak
    for iGSR = 1:nGSR
        for iSeed = 1:nSeed
            figure('Units','normalized','Position',[0.1,0.1,0.5,0.8],'Color','w');
            tiledlayout(nMask, nData, 'TileSpacing', 'normal', 'Padding', 'compact');
            %     sgtitle(['DLPFCSEED ',DLPFCSeeds{iSeed}], 'FontWeight', 'bold');
            for iMask = 1:nMask
                for iData = 1:nData
                    nexttile
                    load([InDirs{iData},filesep,'GSR_',num2str(iGSR-1),filesep,'PEAKMODE_',PeakModes{iPeak},filesep,'MASK_',DLPFCMasks{iMask},filesep,'DLPFCSEED_',DLPFCSeeds{iSeed},filesep,'Summary.mat']);
                    
                    plot(1:CluSizeRange(2), ClinicalSignificance_r, 'LineWidth', 1.5);
                    [~,MinIndex] = min(ClinicalSignificance_r);
                    xline(MinIndex,'--','LineWidth', 1.5, 'Color', [0.5,0,0]);
                    grid on;
                    box on;
                    
                    xlim([1,CluSizeRange(2)]);
                    ylim([-0.6,0]);
                    ylabel({'Clinical Significance', '(Pearson''s Correlation between Targeting', 'Offset and Clinical Improvement)'}, 'FontSize', 9);
                    xlabel({'Threshold for Individualized Targeting Map', '(Number of Voxel)'}, 'FontSize', 9);
%                     title(regexprep(['DLPFCSEED_',DLPFCSeeds{iSeed},filesep,'MASK_',DLPFCMasks{iMask}],'_','-'), 'FontSize', 10, 'FontWeight', 'normal');
                end
            end
            fig = gcf;
            exportgraphics(fig,[OutDir,filesep,'DLPFCSEED_',DLPFCSeeds{iSeed},'_MASK_',DLPFCMasks{iMask},'_GSR_',num2str(iGSR-1),'_PEAKMODE_',PeakModes{iPeak},'.jpg'],'Resolution',300);
        end
    end
end


%% [Compare Paremeters] 4. Compare Different Seeds Derived from the Group Difference Map - The original seed is the best, but the SFG seed is surprisingly not bad
% Fix PeakMode-"Size" and GSR-"No GSR", set BA46 mask for the original seed but set MFG mask for the SFG seed (because only the large MFG mask from AALV3 can include the SFG seed).
clc;clear;close all;

OutDir = '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Code/ForManuscript_TMS_SciBul/Figure/Search_Proper_Target_Threshold/4_Compare_T2-based_Seeds';
mkdir(OutDir);

% Datasets
Datasets = {'SID','CUD'};
InDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/SID',...
    '/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Results/ForManuscript_SciBul_Revision/Search_Proper_Target_Threshold/CUD'};

% Seedmaps
DLPFCSeeds = {'T2-FDR10','T2-SFG-FDR10','T2-SFG+MFG-FDR10'};
% SeedNames = {'MFG','SFG','SFG+MFG'};

    
% DLPFC masks for defining searching scope
DLPFCMasks = {'MFG'};
MaskDirs = {'/mnt/Data6/RfMRILab/Lubin/Project/DR_Targeting/Template/AAL3_MFG/Reslice_MFG_mask.nii'};

% Approaches for defining peak points as targets
PeakModes = {'Size'}; 

% Do GSR or not 
Suffixes = {'FunImgARCWF'};

CluSizeRange = [5,100];

nData = length(Datasets);
nGSR = length(Suffixes);
nPeak = length(PeakModes);
nSeed = length(DLPFCSeeds);
nMask = length(DLPFCMasks);
nStep = CluSizeRange(2)-CluSizeRange(1)+1;

for iPeak = 1:nPeak
    for iGSR = 1:nGSR
        for iMask = 1:nMask
            figure('Units','normalized','Position',[0.1,0.1,0.5,0.8],'Color','w');
            tiledlayout(nSeed, nData, 'TileSpacing', 'normal', 'Padding', 'compact');
            for iSeed = 1:nSeed
                for iData = 1:nData
                    nexttile
                    if iSeed == 1 % Compare with the benchmarker using BA46 as the mask 
                        load([InDirs{iData},filesep,'GSR_',num2str(iGSR-1),filesep,'PEAKMODE_',PeakModes{iPeak},filesep,'MASK_BA46',filesep,'DLPFCSEED_',DLPFCSeeds{iSeed},filesep,'Summary.mat']);
                    else
                        load([InDirs{iData},filesep,'GSR_',num2str(iGSR-1),filesep,'PEAKMODE_',PeakModes{iPeak},filesep,'MASK_',DLPFCMasks{iMask},filesep,'DLPFCSEED_',DLPFCSeeds{iSeed},filesep,'Summary.mat']);
                    end
                    
                    plot(1:CluSizeRange(2), ClinicalSignificance_r, 'LineWidth', 1.5);
                    [~,MinIndex] = min(ClinicalSignificance_r);
                    xline(MinIndex,'--','LineWidth', 1.5, 'Color', [0.5,0,0]);
                    grid on;
                    box on;
                    
                    xlim([1,CluSizeRange(2)]);
                    ylim([-0.6,0]);
                    ylabel({'Clinical Significance', '(Pearson''s Correlation between Targeting', 'Offset and Clinical Improvement)'}, 'FontSize', 9);
                    xlabel({'Threshold for Individualized Targeting Map', '(Number of Voxel)'}, 'FontSize', 9);
                    %                     title(regexprep(['DLPFCSEED_',DLPFCSeeds{iSeed},filesep,'MASK_',DLPFCMasks{iMask}],'_','-'), 'FontSize', 10, 'FontWeight', 'normal');
                end
            end
            fig = gcf;
            exportgraphics(fig,[OutDir,filesep,'Compare_Alternative_T2_Seeds.jpg'],'Resolution',300);
        end
    end
end


