%% Prepare targets for illustraction using BrainNet Viewer
clc;clear;

%% Prepare DR targets %%
Datasets = {'SID','CUD'};
InDir = '/Users/lubin/Documents/Projects/TMS_Targeting/Results/ForManuscript_SciBul_Revision/ForManuscript_SciBul_Revision/Generated_Targets';
OutDir = '/Users/lubin/Documents/Projects/TMS_Targeting/Results/ForManuscript_SciBul_Revision/ForManuscript_SciBul_Revision/Generated_Targets/ForBrainNet';
mkdir(OutDir);

for iData = 1:length(Datasets)
    List = dir([InDir,filesep,Datasets{iData},filesep,'*SEED*']);
    for iPara = 1:length(List)
        load([List(iPara).folder,filesep,List(iPara).name,filesep,'Targets_and_Clinical_Significance.mat']);
        
        NodeSize = HAMD_Reducing_Rate*6;
        fid = fopen([OutDir, filesep, Datasets{iData},'_',List(iPara).name, '.node'],'wt');
        for iTarget = 1:size(Optimized_Targets, 1)
            fprintf(fid,'%f\t%f\t%f\t%s\t%s\t%s\n', ...
                Optimized_Targets(iTarget,1), ...
                Optimized_Targets(iTarget,2), ...
                Optimized_Targets(iTarget,3), ...
                '1', '4', '-');
        end
        fclose(fid);
    end
end

%% Prepare seedmaps in different masks %%
clc;clear;
SeedNames = {'T2','MeanMDD','MeanHC'};
SeedDirs = {'/Users/lubin/Documents/Projects/TMS_Targeting/Template/T2_sgACC-FC_Fox2012/Dx_FC_SeedSurfLHSurfRHVolu_FunVoluWglobalCF.nii',...
    '/Users/lubin/Documents/Projects/TMS_Targeting/Template/Mean_sgACC-FC_Fox2012/group_mean_sgACC_FCMaps_MDD.nii',...
    '/Users/lubin/Documents/Projects/TMS_Targeting/Template/Mean_sgACC-FC_Fox2012/group_mean_sgACC_FCMaps_HC.nii'};
MaskNames = {'BA46','MFG','BA46+9_n_MFG'};
MaskDirs = {'/Users/lubin/Documents/Projects/TMS_Targeting/Template/ForBrainNet/ReslicedMask/Reslice_BA46_mask.nii',...
    '/Users/lubin/Documents/Projects/TMS_Targeting/Template/ForBrainNet/ReslicedMask/Reslice_MFG_mask.nii',...
    '/Users/lubin/Documents/Projects/TMS_Targeting/Template/ForBrainNet/ReslicedMask/Reslice_BA46+9_n_MFG_mask.nii'};
OutDir = '/Users/lubin/Documents/Projects/TMS_Targeting/Template/ForBrainNet/MaskedSeed';
mkdir(OutDir);

for iSeed = 1:length(SeedDirs)
    for iMask = 1:length(MaskDirs)
        [Seed,Header] = y_Read(SeedDirs{iSeed});
        [Mask,~] = y_Read(MaskDirs{iMask});
        Seed = Mask.*Seed;
        y_Write(Seed,Header,[OutDir,filesep,SeedNames{iSeed},'_',MaskNames{iMask}]);
    end
end
        

    
    
    
    
    
    
    
    
    