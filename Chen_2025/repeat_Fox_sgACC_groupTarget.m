% Fox first used a group of MDD patients' sgACC FC maps to get a group
% level target. see Fox et al., 2012
% Then Weigand et al. upgraded to a group of 1000 NC's sgACC FC maps to get
% the most anti-correlated point in DLPFC

%% initialization
clear;clc;
work_dir = '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/Fox_Replication';
combat_FC_dir = '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/group_diff_maps/comBat';

load /mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/Sub_Info_1574_1308.mat

MaskSet = {'Fox_2012_Grey_Matter' ...           
                      };
MeasureStringVoluSet = {'FC_SeedSurfLHSurfRHVolu_FunVoluWCF', ...
                        'FC_SeedSurfLHSurfRHVolu_FunVoluWglobalCF'};

%% surface
[mask,~,~,header_mask] = y_ReadAll([work_dir,'/Vol2Surf_DLPFC_Mask_lh.gii']);

% MDD
[data,~,~,header] = y_ReadAll([work_dir, '/group_mean_sgACC_FCMaps_MDD_lh.gii']);
data_DLPFC = data.*mask;
% y_Write(data_DLPFC, header,'data_DLPFC_MDD.gii');
data_target = zeros(length(data_DLPFC),1);
data_target(find(data_DLPFC == min(data_DLPFC))) = 1;
y_Write(data_target,header,'sgACC_groupTarget_MDD.gii');

% HC
[data,~,~,header] = y_ReadAll([work_dir, '/group_mean_sgACC_FCMaps_NC_lh.gii']);
data_DLPFC = data.*mask;
% y_Write(data_DLPFC, header,'data_DLPFC_HC.gii');
data_target = zeros(length(data_DLPFC),1);
data_target(find(data_DLPFC == min(data_DLPFC))) = 1;
y_Write(data_target,header,'sgACC_groupTarget_HC.gii');

%% volume
ROI = MaskSet{1};
MeasureString = MeasureStringVoluSet{2};

% calculate MDD/HC mean map
wanted_matrix = zeros(length(SubID),1);
% wanted_matrix(Dx == 1) = 1;
wanted_matrix(Dx == -1) = 1;
wanted_idx = find(wanted_matrix);
SubID_Grouped = SubID(wanted_idx);
FileList = cell(length(SubID_Grouped),1);
for iSub=1:length(SubID_Grouped)
    FileList{iSub,1}=sprintf('%s/Volu/%s_%s/szFC_%s.nii',combat_FC_dir, ...
        MeasureString,ROI,SubID_Grouped{iSub}); 
end

[data,~,~,header] = y_ReadAll(FileList);
data_mean = mean(data,4);
% y_Write(data_mean, header, [work_dir, '/group_mean_sgACC_FCMaps_MDD.nii']);
y_Write(data_mean, header, [work_dir, '/group_mean_sgACC_FCMaps_HC.nii']);


% reslice the mask
InputFile = [work_dir,'/DLPFC_Mask.nii'];
OutputFile = [work_dir,'/DLPFC_Mask_Reslice.nii'];
NewVoxSize = [];
hld = 0;
TargetSpace = [work_dir, '/group_mean_sgACC_FCMaps_MDD.nii'];
y_Reslice(InputFile,OutputFile,NewVoxSize,hld, TargetSpace);

[mask,~,~,header_mask] = y_ReadAll([work_dir,'/DLPFC_Mask_Reslice.nii']);

% MDD
[data,~,~,header] = y_ReadAll([work_dir, '/group_mean_sgACC_FCMaps_MDD.nii']);
data_DLPFC = data.*mask;
data_DLPFC_2d = data_DLPFC(:);

data_target = zeros(length(data_DLPFC_2d),1);
data_target(find(data_DLPFC_2d == min(data_DLPFC_2d))) = 1;

data_dimension = size(data);
data_target_3d = reshape(data_target,data_dimension);

y_Write(data_target_3d,header,'sgACC_groupTarget_MDD.nii');

% HC
[data,~,~,header] = y_ReadAll([work_dir, '/group_mean_sgACC_FCMaps_HC.nii']);
data_DLPFC = data.*mask;
data_DLPFC_2d = data_DLPFC(:);

data_target = zeros(length(data_DLPFC_2d),1);
data_target(find(data_DLPFC_2d == min(data_DLPFC_2d))) = 1;

data_dimension = size(data);
data_target_3d = reshape(data_target,data_dimension);

y_Write(data_target_3d,header,'sgACC_groupTarget_HC.nii');
