% do Ward's hierarchical clustering with FC maps calculated with different
% sgACC ROIs
%
% Xiao Chen 221019
% chenxiaophd@gmail.com

%% initialization
clear; clc;

InputDir = '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/group_diff_maps';
OutputDir = '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/Hierarchical_Cluter';
if ~exist(OutputDir, "dir"); mkdir(OutputDir); end
                  
MeasureStringSet = {'FC_SeedSurfLHSurfRHVolu_FunSurfWCF', ...
                    'FC_SeedSurfLHSurfRHVolu_FunSurfWglobalCF'};
MeasureStringVoluSet = {'FC_SeedSurfLHSurfRHVolu_FunVoluWCF', ...
                        'FC_SeedSurfLHSurfRHVolu_FunVoluWglobalCF'};

ROI = 'DLPFC';
                    
%% ROI on surface to surface maps

% % All sgACC masks
% MaskSet = {'Brodmann_Mai_Matajnik_BA25_fsaverage5_lh',  ...
%            'glasser_BA25_lh', 'Schaefer_sgACC_lh', ...
%            'Brodmann_Mai_Matajnik_BA25_fsaverage5_rh', ...
%            'glasser_BA25_rh', 'Schaefer_sgACC_rh', ...
%            'Fox_2012_Grey_Matter_fsaverage5_masked_rh', ...
%            'Liston_2014_Grey_Matter_fsaverage5_masked_rh', ...
%            'Brodmann_Mai_Matajnik_BA25_lh', 'Brodmann_Mai_Matajnik_BA25_rh',...
%            'Fox_2012_Grey_Matter', 'glasser_BA25_MNI_lh', ...
%            'glasser_BA25_MNI_rh', ...
%            'Liston_2014_Grey_Matter'
%                       };

% % All DLPFC masks on the surface
% MaskSet = {'5cmMethod_20mm_Grey_Matter_fsaverage5_masked_lh',  ...
%            'BA_9_15mm_Grey_Matter_fsaverage5_masked_lh', ...
%            'BA_9_20mm_Grey_Matter_fsaverage5_masked_lh', ...
%            'BA_46_20mm_Grey_Matter_fsaverage5_masked_lh', ...
%            'F3Beam_20mm_Grey_Matter_fsaverage5_masked_lh', ...
%            'SGC_group_target_Grey_Matter_fsaverage5_masked_lh'
%                       };
                  
% All DLPFC masks in the volume
MaskSet = {'5cmMethod_20mm_Grey_Matter',  ...
           'BA_9_20mm_Grey_Matter', ...
           'BA_46_20mm_Grey_Matter', ...
           'F3Beam_20mm_Grey_Matter', ...
           'SGC_group_target_Grey_Matter', ...
           'group_diff_target_V2_C2_Grey_Matter', ...
           'group_diff_target_V2_C1_Grey_Matter'
                      };

%%                 
MeasureString = MeasureStringSet{2};
for iMask = 1:length(MaskSet)
    [LH, ~, ~, ~] = y_ReadAll([InputDir,'/stats_comBat/FuncSurfLH/', ...
                                MaskSet{iMask},'/Dx_', MeasureString, '.gii']);
    [RH, ~, ~, ~] = y_ReadAll([InputDir,'/stats_comBat/FuncSurfRH/', ...
                            MaskSet{iMask},'/Dx_', MeasureString, '.gii']);
    current_vector = [LH; RH];
    if iMask == 1
        cluster_matrix = current_vector;
    else
        cluster_matrix = [cluster_matrix,current_vector];                   
    end
end

x = cluster_matrix';
D = pdist(x, "correlation");
Z = linkage(D, "Ward");
eva = evalclusters(x, 'Linkage', 'gap', 'KList', [1:6]);
plot(eva)

figure;
set(gcf, "position", [150, 80, 1000, 700]);
dendrogram(Z, "Labels", MaskSet, "Orientation", "Left");
set(gca, "fontsize", 18);
print(gcf, [OutputDir,'/',ROI,'_ROI_SurfaceMaps.jpg'], "-djpeg", "-r300");
close all;

%% ROI on the surface to volume space
MeasureString = MeasureStringVoluSet{2}; 
for iMask = 1:length(MaskSet)
    [Data,~] = y_ReadAll([InputDir,'/stats_comBat/FunVolu/', ...
                        MaskSet{iMask},'/Dx_', MeasureString, '.nii']);
    current_vector = Data(:);
    if iMask == 1
        cluster_matrix = current_vector;
    else
        cluster_matrix = [cluster_matrix,current_vector];                   
    end
end

D = pdist(cluster_matrix', "correlation");
Z = linkage(D, "Ward");

figure;
set(gcf, "position", [150, 80, 1000, 700]);
dendrogram(Z, "Labels", MaskSet, "Orientation", "Left");
set(gca, "fontsize", 18);
print(gcf, [OutputDir,'/',ROI,'_SurfaceROI_VolumeMaps.jpg'], "-djpeg", "-r300");
close all;

%% ROI in the volume to surface maps
% % All sgACC masks
% MaskSet = {'Brodmann_Mai_Matajnik_BA25_lh', 'Brodmann_Mai_Matajnik_BA25_rh',...
%            'Fox_2012_Grey_Matter', 'glasser_BA25_MNI_lh', ...
%            'glasser_BA25_MNI_rh', ...
%            'Liston_2014_Grey_Matter'
%                       };
                  
% % All DLPFC masks
% MaskSet = {'5cmMethod_20mm_Grey_Matter', 'BA_9_9mm_Grey_Matter', ...
%            'BA_9_15mm_Grey_Matter', 'BA_9_20mm_Grey_Matter', ...
%            'BA_46_20mm_Grey_Matter', 'F3Beam_20mm_Grey_Matter', ...
%            'SGC_group_target_Grey_Matter'
%                       };
                  
MeasureString = MeasureStringSet{2};
for iMask = 1:length(MaskSet)
    [LH, ~, ~, ~] = y_ReadAll([InputDir,'/stats_comBat/FuncSurfLH/', ...
                                MaskSet{iMask},'/Dx_', MeasureString, '.gii']);
    [RH, ~, ~, ~] = y_ReadAll([InputDir,'/stats_comBat/FuncSurfRH/', ...
                            MaskSet{iMask},'/Dx_', MeasureString, '.gii']);
    current_vector = [LH; RH];
    if iMask == 1
        cluster_matrix = current_vector;
    else
        cluster_matrix = [cluster_matrix,current_vector];                   
    end
end

D = pdist(cluster_matrix', "correlation");
Z = linkage(D, "Ward");

figure;
set(gcf, "position", [150, 80, 1000, 700]);
dendrogram(Z, "Labels", MaskSet, "Orientation", "Left");
set(gca, "fontsize", 18);
print(gcf, [OutputDir,'/',ROI,'_VolumeROI_SurfaceMaps.jpg'], "-djpeg", "-r300");
close all;

%% ROI in the volume space to Volume maps
MeasureString = MeasureStringVoluSet{2}; 
for iMask = 1:length(MaskSet)
    [Data,~] = y_ReadAll([InputDir,'/stats_comBat/FunVolu/', ...
                        MaskSet{iMask},'/Dx_', MeasureString, '.nii']);
    current_vector = Data(:);
    if iMask == 1
        cluster_matrix = current_vector;
    else
        cluster_matrix = [cluster_matrix,current_vector];                   
    end
end

output_matrix = cluster_matrix;
T = array2table(output_matrix);
T.Properties.VariableNames(1:7) = MaskSet;
writetable(T, [OutputDir,'/cluster_matrix_2targets.csv']);

x = cluster_matrix';
x4zero = sum(x,1);
x(:,find(x4zero == 0)) = [];
D = pdist(x, "correlation");
Z = linkage(D, "Ward");
% use gap statistics to evaluate the results
x = single(x);
eva = evalclusters(x, 'linkage', 'gap', 'klist', [1:6]);
plot(eva)
print(gcf, [OutputDir,'/Gap_statistics.jpg'], "-djpeg", "-r300");

figure;
set(gcf, "position", [150, 80, 1000, 700]);
dendrogram(Z, "Labels", MaskSet, "Orientation", "Left");
set(gca, "fontsize", 18);
print(gcf, [OutputDir,'/',ROI,'_VolumeROI_VolumeMaps.jpg'], "-djpeg", "-r300");
close all;