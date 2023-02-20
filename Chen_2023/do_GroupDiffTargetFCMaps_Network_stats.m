% To first extract the 7 network of the FC maps from the group difference
% FC target and then conducted glm analyses on them
% Xiao Chen 221103
% chenxiaophd@gmail.com

%% initialization

clear;clc;

combat_FC_dir = '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/group_diff_target/comBat';
                 
output_dir = '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/Network_Analyses';

load '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/Sub_Info_1574_1308.mat'

% MaskSet = {'Brodmann_Mai_Matajnik_BA25_lh', 'Brodmann_Mai_Matajnik_BA25_rh',...
%            'Fox_2012_Grey_Matter', 'glasser_BA25_MNI_lh', ...
%            'glasser_BA25_MNI_rh', ...
%            'Liston_2014_Grey_Matter', 'Schaefer_sgACC_lh', ...
%            'Schaefer_sgACC_rh',...
%            '5cmMethod_20mm_Grey_Matter',... 
%            'BA_9_20mm_Grey_Matter', ...
%            'BA_46_20mm_Grey_Matter', 'F3Beam_20mm_Grey_Matter', ...
%            'SGC_group_target_Grey_Matter'            
%                         };
                    
% % select some masks
% MaskSet = {'Fox_2012_Grey_Matter', ...
%            'Liston_2014_Grey_Matter', ...
%            '5cmMethod_20mm_Grey_Matter',... 
%            'BA_9_20mm_Grey_Matter', ...
%            'BA_46_20mm_Grey_Matter', 'F3Beam_20mm_Grey_Matter', ...
%            'SGC_group_target_Grey_Matter'
%                         };
                    
% group difference target
MaskSet = {'group_diff_target_V2_C1_Grey_Matter', 'group_diff_target_V2_C2_Grey_Matter'};

MeasureStringSet = {'FC_SeedSurfLHSurfRHVolu_FunSurfWCF', ...
                    'FC_SeedSurfLHSurfRHVolu_FunSurfWglobalCF'};
MeasureStringVoluSet = {'FC_SeedSurfLHSurfRHVolu_FunVoluWCF', ...
                        'FC_SeedSurfLHSurfRHVolu_FunVoluWglobalCF'};

%% select subjects and consutruct matrix
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

%% extract network FCs on the surface
TemplateDir = '/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Analysis/Matrix_Extraction';
% TemplateDir = '/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Analysis/Matrix_Extraction_Yeo2011';
dpabiSurfTemplateDir = '/mnt/Data2/RfMRILab/Chenxiao/CX_software/DPABI-master/DPABISurf/SurfTemplates';
HemisphereNameSet = {'lh', 'rh'};
HemisphereSet = {'LH', 'RH'};

TopOutputDir = [output_dir, '/FC_network_extracted_Schaefer400Surface'];
if ~exist(TopOutputDir); mkdir(TopOutputDir); end
% surface
IsMultipleLabel = 1;
GHeader = [];
CUTNUMBER = 1;
ROIDef = {};

for iMeasure  = 1:2
    MeasureString = MeasureStringSet{iMeasure};
    for iHem = 1:2
        for iMask = 1:length(MaskSet)
            FileList = cell(length(SubID),1);
            for iSub = 1:length(SubID)
                FileList{iSub,1}=sprintf('%s/%s/%s_%s/szFC_%s.func.gii',combat_FC_dir, ...
                                         HemisphereSet{iHem},MeasureString,MaskSet{iMask},SubID{iSub});
            end
            [AllVolume,~,~,~] = y_ReadAll(FileList); 
            ROIDef{1} = [TemplateDir, '/fsaverage5_', HemisphereNameSet{iHem}, ...
                                                    '_Schaefer2018_400Parcels_17Networks_order.label.gii'];
%             ROIDef{1} = [TemplateDir, '/fsaverage5_', HemisphereNameSet{iHem}, ...
%                                                     '_Yeo2011_17Networks_order.label.gii'];
            OutputName = [TopOutputDir, '/', HemisphereSet{iHem}, '/', MeasureString, '_', ...
                           MaskSet{iMask} ];
            if ~exist([TopOutputDir, '/', HemisphereSet{iHem}], "dir")
                mkdir([TopOutputDir, '/', HemisphereSet{iHem}])
            end
            AMaskFilename = [dpabiSurfTemplateDir, '/fsaverage5_', ...
                                            HemisphereNameSet{iHem}, '_cortex.label.gii'];
            [ROISignals] = y_ExtractROISignal_Surf(AllVolume, ROIDef, OutputName, AMaskFilename, ...
                                                                        IsMultipleLabel, GHeader, CUTNUMBER);
        end
    end
end

%% get network FC
load([output_dir,'/Schaefer2018_400_Info.mat']);
% load([output_dir,'/DPABISurf_Yeo2011_Info.mat']);

% network_FC: n x 9 matrixs, 6 Yeo network and 3 subsystems of DMN
TopOutputDir = [output_dir, '/FC_network_extracted_Schaefer400Surface'];

MeasureString = MeasureStringSet{1};
network_FC = [];
ROISignals_full = [];
ROISignals_full_int = [];
for iMask = 1:length(MaskSet)
    ROI = ['A_',MaskSet{iMask}];
    load([TopOutputDir, '/LH/ROISignals_', MeasureString, '_', ...
                           MaskSet{iMask}, '.mat']);
    ROISignals_lh = ROISignals;
    load([TopOutputDir, '/RH/ROISignals_', MeasureString, '_', ...
                           MaskSet{iMask}, '.mat']);
    ROISignals_rh = ROISignals;
    ROISignals_full.(ROI) = [ROISignals_lh,ROISignals_rh];
    ROISignals_full_int.(ROI) = ROISignals_full.(ROI)(:, ROIIndex_181_Schafer);
    network_label = DPABISurf_Schaefer2018_400_YeoNetwork';
%     ROISignals_full_int.(ROI) = ROISignals_full.(ROI)(:, ROIIndex_Yeo2011);
%     network_label = DPABISurf_Yeo2011_YeoNetwork';
    for iNetwork = unique(network_label)
        network_FC.(ROI)(:,iNetwork) = mean(ROISignals_full_int.(ROI)(:, network_label == iNetwork), 2);
    end
    network_FC.(ROI)(:,10) = mean(network_FC.(ROI)(:,7:9),2);
end

MeasureString = MeasureStringSet{2};
network_FC_GSR = [];
ROISignals_full = [];
for iMask = 1:length(MaskSet)
    ROI = ['A_',MaskSet{iMask}];
    load([TopOutputDir, '/LH/ROISignals_', MeasureString, '_', ...
                           MaskSet{iMask}, '.mat']);
    ROISignals_lh = ROISignals;
    load([TopOutputDir, '/RH/ROISignals_', MeasureString, '_', ...
                           MaskSet{iMask}, '.mat']);
    ROISignals_rh = ROISignals;
    ROISignals_full.(ROI) = [ROISignals_lh,ROISignals_rh];
    ROISignals_full_int.(ROI) = ROISignals_full.(ROI)(:, ROIIndex_181_Schafer);
    network_label = DPABISurf_Schaefer2018_400_YeoNetwork';
%     ROISignals_full_int.(ROI) = ROISignals_full.(ROI)(:, ROIIndex_Yeo2011);
%     network_label = DPABISurf_Yeo2011_YeoNetwork';
    for iNetwork = unique(network_label)
        network_FC_GSR.(ROI)(:,iNetwork) = mean(ROISignals_full_int.(ROI)(:, network_label == iNetwork), 2);
    end
    network_FC_GSR.(ROI)(:,10) = mean(network_FC_GSR.(ROI)(:,7:9),2);
end

%% open par pool
poolobj = gcp('nocreate');
delete(poolobj);
parpool(2);

%% extract network FC in the volume space
TemplateDir = '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Analyses_V2/Network_Analyses';
dpabiSurfTemplateDir = '/mnt/Data5/RfMRILab/Chenxiao/sgACC_Project/Data/Masks';

TopOutputDir = [output_dir, '/FC_network_extracted_Schaefer400Volume_7Network'];
if ~exist(TopOutputDir); mkdir(TopOutputDir); end

% %Reslice mask
% InputFile = [TemplateDir, '/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii'];
% OutputFile = [TemplateDir, '/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm_Reslice.nii'];
% NewVoxSize = [];
% hld = 0;
% TargetSpace = [dpabiSurfTemplateDir, '/AllResampled_BrainMask_05_91x109x91.nii'];
% y_Reslice(InputFile,OutputFile,NewVoxSize,hld, TargetSpace);


ROIDef = {};
ROIDef{1} = [TemplateDir, '/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm_Reslice.nii'];
for iMeasure  = 1:2
    MeasureString = MeasureStringVoluSet{iMeasure};
    parfor iMask = 1:length(MaskSet)
        FileList = cell(length(SubID),1);
        for iSub = 1:length(SubID)
            FileList{iSub,1}=sprintf('%s/Volu/%s_%s/szFC_%s.nii',combat_FC_dir, ...
                                     MeasureString,MaskSet{iMask},SubID{iSub});
        end
        [AllVolume,~,~,Header] = y_ReadAll(FileList); 
       
        OutputName = [TopOutputDir, '/', MeasureString, '_', MaskSet{iMask} ];
       
        AMaskFilename = [dpabiSurfTemplateDir, '/AllResampled_BrainMask_05_91x109x91.nii'];
        [ROISignals] = y_ExtractROISignal(AllVolume, ROIDef, OutputName, AMaskFilename, 1, ...
                                          0,[],[],[],[],[],Header);
    end
end

%% get network FC
load([output_dir,'/Schaefer2018_400_7Network_Info.mat']);
% load([output_dir,'/Schaefer2018_400_Info.mat']);
% load([output_dir,'/DPABISurf_Yeo2011_Info.mat']);

% network_FC: n x 7 matrixs, 7 Yeo networks
MeasureString = MeasureStringVoluSet{1};
network_FC = [];
ROISignals_full=[];
ROISignals_full_int = [];
for iMask = 1:length(MaskSet)
    ROI = ['A_',MaskSet{iMask}];
    load([TopOutputDir, '/ROISignals_', MeasureString, '_', ...
                           MaskSet{iMask}, '.mat']);
    ROISignals_full.(ROI) = ROISignals;
    ROISignals_full_int.(ROI) = ROISignals_full.(ROI)(:, ROIIndex_181_Schafer);
    network_label = DPABISurf_Schaefer2018_400_YeoNetwork';
%     ROISignals_full_int.(ROI) = ROISignals_full.(ROI)(:, ROIIndex_Yeo2011);
%     network_label = DPABISurf_Yeo2011_YeoNetwork';
    for iNetwork = unique(network_label)
        network_FC.(ROI)(:,iNetwork) = mean(ROISignals_full_int.(ROI)(:, network_label == iNetwork), 2);
    end
end

MeasureString = MeasureStringVoluSet{2};
network_FC_GSR = [];
ROISignals_full = [];
ROISignals_full_int = [];
for iMask = 1:length(MaskSet)
    ROI = ['A_',MaskSet{iMask}];
    load([TopOutputDir, '/ROISignals_', MeasureString, '_', ...
                           MaskSet{iMask}, '.mat']);
    ROISignals_full.(ROI) = ROISignals;
    ROISignals_full_int.(ROI) = ROISignals_full.(ROI)(:, ROIIndex_181_Schafer);
    network_label = DPABISurf_Schaefer2018_400_YeoNetwork';
%     ROISignals_full_int.(ROI) = ROISignals_full.(ROI)(:, ROIIndex_Yeo2011);
%     network_label = DPABISurf_Yeo2011_YeoNetwork';
    for iNetwork = unique(network_label)
        network_FC_GSR.(ROI)(:,iNetwork) = mean(ROISignals_full_int.(ROI)(:, network_label == iNetwork), 2);
    end
end

%% write out .csv files for ploting

for iMask = 1:length(MaskSet)
    ROI = ['A_',MaskSet{iMask}];
    output_matrix = [network_FC.(ROI), Dx];
    T = array2table(output_matrix, ...
                'VariableNames', {'Network1', 'Network2', 'Network3', 'Network4', 'Network5', ...
                                  'Network6', 'Network7', 'Group'});
    T_stacked = stack(T, {'Network1', 'Network2', 'Network3', 'Network4', 'Network5', ...
                                  'Network6', 'Network7'}, ...
                            'NewDataVariableName','FC','IndexVariableName', 'Network');
    writetable(T_stacked, [output_dir,'/',MaskSet{iMask},'_4ViolinPlot.csv']);
    
    ROI = ['A_',MaskSet{iMask}];
    output_matrix = [network_FC_GSR.(ROI), Dx];
    T = array2table(output_matrix, ...
                'VariableNames', {'Network1', 'Network2', 'Network3', 'Network4', 'Network5', ...
                                  'Network6', 'Network7', 'Group'});
    T_stacked = stack(T, {'Network1', 'Network2', 'Network3', 'Network4', 'Network5', ...
                                  'Network6', 'Network7'}, ...
                            'NewDataVariableName','FC','IndexVariableName', 'Network');
    writetable(T_stacked, [output_dir,'/',MaskSet{iMask},'_4ViolinPlot_GSR.csv']);
end

%% do stats on group differences
TF_Flag = 'T';
stats = [];
for iMask = 1:length(MaskSet)
    ROI = ['A_',MaskSet{iMask}];
    for iNetwork = 1:size(network_FC.(ROI),2)
        y = network_FC.(ROI)(:, iNetwork);
        X = AllCov;
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,X,Contrast,TF_Flag);
        stats.(ROI)(iNetwork,1) = TF_ForContrast;
        stats.(ROI)(iNetwork,2) = 2*tcdf(-abs(TF_ForContrast),DOF);
        stats.(ROI)(iNetwork,3) = Cohen_f2;
        stats.(ROI)(iNetwork,4) = 14*tcdf(-abs(TF_ForContrast),DOF);
    end
end

TF_Flag = 'T';
stats_GSR = [];
for iMask = 1:length(MaskSet)
    ROI = ['A_',MaskSet{iMask}];
    for iNetwork = 1:size(network_FC_GSR.(ROI),2)
        y = network_FC_GSR.(ROI)(:, iNetwork);
        X = AllCov;
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,X,Contrast,TF_Flag);
        stats_GSR.(ROI)(iNetwork,1) = TF_ForContrast;
        stats_GSR.(ROI)(iNetwork,2) = 2*tcdf(-abs(TF_ForContrast),DOF);
        stats_GSR.(ROI)(iNetwork,3) = Cohen_f2;
        stats_GSR.(ROI)(iNetwork,4) = 14*tcdf(-abs(TF_ForContrast),DOF);
    end
end

save([output_dir,'/t_stats_volu_7networks_groupDiff_target.mat'],'stats','stats_GSR','-mat');

%% and do correlation analyses between HAMD scores and network FCs
load /mnt/Data7/RfMRILab/Yan/YAN_Work/REST_meta-MDD_Surf/Analysis/SubInfo/AllPhenotype.mat
load /mnt/Data7/RfMRILab/Yan/YAN_Work/REST_meta-MDD_Surf/Analysis/SubCorticalVolume/Data_MDD_NC.mat

SubID_full =[MDD.ID;NC.ID];
HAMD_total = [MDD.HAMDTotal17; zeros(length(NC.ID),1)];
HAMA_total = [MDD.HAMATotal; zeros(length(NC.ID),1)];
HAMD_item = [table2array(MDD(:,18:34)); zeros(length(NC.ID), size(18:34,2))];
HAMA_item = [table2array(MDD(:,43:56)); zeros(length(NC.ID), size(43:56,2))];

% select participants
wanted_matrix = zeros(length(SubID_full),1);
for iSub = 1:length(SubID_full)
   if ismember(SubID_full{iSub}, SubID)
       wanted_matrix(iSub,1) = 1;
   end
end

wanted_idx = find(wanted_matrix);
HAMD_total = HAMD_total(wanted_idx);
HAMA_total = HAMA_total(wanted_idx);
HAMD_item = HAMD_item(wanted_idx,:);
HAMA_item = HAMA_item(wanted_idx,:);

% only contain MDD patients
HAMD_total = HAMD_total(Dx == 1);
HAMA_total = HAMA_total(Dx == 1);
HAMD_item = HAMD_item(Dx == 1,:);
HAMA_item = HAMA_item(Dx == 1,:);

Age = Age(Dx == 1);
Sex = Sex(Dx == 1);
Edu = Edu(Dx == 1);
HeadMotion = HeadMotion(Dx == 1);

%% do stats
% clinical_label = 'HAMD_total';
% y = HAMD_total;

% clinical_label = 'HAMA_total';
% y = HAMA_total;

% clinical_label = 'HAMD_item';
% y = HAMD_item;

clinical_label = 'HAMA_item';
y = HAMA_item;

y4NaN = sum(y,2);
wanted_matrix = ones(length(y4NaN),1);
wanted_matrix(isnan(y4NaN)) = 0;
wanted_idx = find(wanted_matrix);
y = y(wanted_idx,:);

TF_Flag = 'T';
stats_GSR = [];
for iMask = 1:length(MaskSet)
    ROI = ['A_',MaskSet{iMask}];
    for iNetwork = 1:size(network_FC_GSR.(ROI),2)
        for iItem = 1:size(y,2)
           x =  network_FC_GSR.(ROI)(:,iNetwork);
           x =  x(Dx == 1);
           AllCov = [ones(length(x),1), x, Age, Sex, Edu, HeadMotion];
           AllCov = AllCov(wanted_idx,:);
           Contrast = zeros(1,size(AllCov,2));
           Contrast(2) = 1;
           DOF = size(AllCov,1) - size(AllCov,2);
           [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y(:,iItem),AllCov,Contrast,TF_Flag);
           stats_GSR.(ROI){iItem,1}(iNetwork,1) = TF_ForContrast;
           stats_GSR.(ROI){iItem,1}(iNetwork,2) = 2*tcdf(-abs(TF_ForContrast),DOF);
           stats_GSR.(ROI){iItem,1}(iNetwork,3) = Cohen_f2;
        end
    end
end

save([output_dir,'/', clinical_label,'_GSR.mat'],'stats_GSR','-mat');