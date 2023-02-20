% calculate the Euclidean distances between previous spot and its
% correlation to the clinical outcomes
% this script is also used to calculate the correlation between group
% difference t and the hamd reduction rate
%
%
% Xiao Chen 221019
% chenxiaophd@gmail.com

%% initialization
clear; clc;

OutputDir = '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/group_diff_target';
load([OutputDir,'/group_diff_target_Fox2012.mat']);
previous_results = readtable([OutputDir,'/previous_clinical_outcomes_V2.txt'], ...
                             "Delimiter", "\t");

%% calculate the Euclidean distances

prev_coord = [previous_results.x, previous_results.y, previous_results.z];
HAMD = previous_results.HAMD_percent;
study = previous_results.studies;

idx =find(study == 1);
prev_coord = prev_coord(idx,:);
HAMD = HAMD(idx);

for iCoord = 1:size(target_coord,1)
    D = sqrt(sum((prev_coord - target_coord(iCoord,:)).^2, 2));
    [RHO, PVAL] = corr(D, HAMD);
    stats(iCoord,1) = RHO;
    stats(iCoord,2) = PVAL;
end

%%
D = sqrt(sum((prev_coord - Fox2012(1,:)).^2, 2));
figure;
set(gcf, "position", [150, 80, 800, 600]);
plot(D, previous_results.HAMD_percent, "o");
lsline;
ylabel("HAMD(%)");
xlabel("Euclidean Distance");
set(gca, "fontsize", 18);
print(gcf, [OutputDir,'/','Fox_ROI1.jpg'], "-djpeg", "-r300");
close all;

%% group diff t maps and HAMD reduction rate
groupDiffMap = ['/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/group_diff_maps/', ...
                'stats_comBat/FunVolu/Fox_2012_Grey_Matter/Dx_FC_SeedSurfLHSurfRHVolu_FunVoluWglobalCF.nii'];
HAMD = previous_results.HAMD_percent;
study = previous_results.studies;
prev_coord = [previous_results.x, previous_results.y, previous_results.z];

idx =find(study == 1);
prev_coord = prev_coord(idx,:);
HAMD = HAMD(idx);
            
radius = 4* ones(length(prev_coord), 1);
sphere_def = [prev_coord, radius];
ROIDef = [];
for i = 1:length(sphere_def)
    ROIDef{i,1} =  sphere_def(i,:); 
end

[ROISignals] = y_ExtractROISignal( ...
                groupDiffMap, ...
                ROIDef, [], '', 0 ); % no brain mask

[RHO, PVAL] = corr(ROISignals', HAMD);

output_matrix = [ROISignals', HAMD];
T = array2table(output_matrix);
T.Properties.VariableNames(1:2) = {'T_Value', 'HAMD_Rate'};
writetable(T, [OutputDir,'/t_HAMDcorr.csv']);
