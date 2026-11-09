clear; clc;

poolsize = 30;
p = gcp('nocreate');
if isempty(p)
    parpool('local', poolsize);
elseif p.NumWorkers ~= poolsize
    delete(p);
    parpool('local', poolsize);
end

% -------------------- Paths --------------------
MDD_root = '/mnt/Data8/RfMRILab/gaoqinglin/2_TMS/DIRECT/AMN_Surf/AMN_DMN/ROISignals/AMN_DMN/ROISignal_TP230/MDD';
CAP_HC_path = '/mnt/Data8/RfMRILab/gaoqinglin/2_TMS/DIRECT/AMN_Surf/AMN_DMN/CAP/NC/CAP_Results/';
out_root = '/mnt/Data8/RfMRILab/gaoqinglin/2_TMS/DIRECT/AMN_Surf/AMN_DMN/CAP/MDD/CAP_Results/CAP_MatchedToRef';
if ~exist(out_root, 'dir'); mkdir(out_root); end

% -------------------- Params --------------------
roiN = 112;
K_range = 2:5;                  
distance_func = 'correlation';    
TR = 2.0;                        

% -------------------- Load MDD ROI time series --------------------
dataDirCandidates = {fullfile(MDD_root,'ROISignal_TP230'), MDD_root};
dataDir = '';
for i = 1:numel(dataDirCandidates)
    if exist(dataDirCandidates{i}, 'dir')
        dataDir = dataDirCandidates{i};
        break;
    end
end
assert(~isempty(dataDir), 'Cannot find MDD data dir.');

subjfiles = dir(fullfile(dataDir, '*.mat'));
if isempty(subjfiles)
    subjfiles = dir(fullfile(dataDir, '**', '*.mat')); 
end
assert(~isempty(subjfiles), 'No .mat files found under: %s', dataDir);

sub_num = numel(subjfiles);
fprintf('Found %d MDD subjects under %s\n', sub_num, dataDir);

roi_signal_all = [];
TP_sub = [];
subj_names = cell(sub_num,1);

for s = 1:sub_num
    fpath = fullfile(subjfiles(s).folder, subjfiles(s).name);
    subj_names{s} = subjfiles(s).name;

    X = load_roi_ts_mat(fpath, roiN);   % T x 112
    Xz = zscore(X);                     

    roi_signal_all = [roi_signal_all; Xz]; %#ok<AGROW>
    TP_sub = [TP_sub; repmat(s, size(Xz,1), 1)]; %#ok<AGROW>
end

save(fullfile(out_root,'MDD_roi_signal_all.mat'), ...
    'roi_signal_all','TP_sub','subj_names','TR','roiN','-v7.3');

% -------------------- Match + Metrics for each K --------------------
for K = K_range
    fprintf('\nMatching MDD frames to reference CAP states (K=%d)...\n', K);

    clust = CAP_Kmeans_corr(roi_signal_all, K, distance_func, CAP_HC_path); 

    CAP_matrices = CAP_CalculateMatrices(clust, TP_sub, K, sub_num);

   

    Metrics = struct();
    Metrics.occurrence_rate = CAP_matrices.fraction;                 
    Metrics.persistence_time_sec = CAP_matrices.persistence * TR;    
    Metrics.transition_probability = CAP_matrices.transitions_prob;  
    clust_by_subj = cell(sub_num,1);
    for s = 1:sub_num
        idx = (TP_sub == s);
        clust_by_subj{s} = clust(idx);
    end
  

    outDirK = fullfile(out_root, ['CAP' num2str(K)]);
    if ~exist(outDirK,'dir'); mkdir(outDirK); end

    save(fullfile(outDirK,'MDD_matched_clust.mat'), ...
        'clust','clust_by_subj','distance_func');

    save(fullfile(outDirK,'CAP_matrices.mat'), 'CAP_matrices', 'Metrics');

    
    [groupSummary, transMean] = summarize_metrics(Metrics);
    save(fullfile(outDirK,'GroupSummary.mat'), 'groupSummary','transMean');

    writetable(groupSummary, fullfile(outDirK,'GroupSummary.csv'));
end

fprintf('\nAll done. Results saved under: %s\n', out_root);



function [tbl, transMean] = summarize_metrics(Metrics)
% Metrics.occurrence_rate: sub x K
% Metrics.persistence_time_sec: sub x K
% Metrics.transition_probability: sub x K x K

occ = Metrics.occurrence_rate;
per = Metrics.persistence_time_sec;
[sub_num, K] = size(occ);

occ_mean = mean(occ, 1, 'omitnan');
occ_sem  = std(occ, 0, 1, 'omitnan') / sqrt(sub_num);

per_mean = mean(per, 1, 'omitnan');
per_sem  = std(per, 0, 1, 'omitnan') / sqrt(sub_num);

State = (1:K).';
tbl = table(State, occ_mean(:), occ_sem(:), per_mean(:), per_sem(:), ...
    'VariableNames', {'State','OccMean','OccSEM','PersistSecMean','PersistSecSEM'});

tp = Metrics.transition_probability;   % sub x K x K
transMean = squeeze(mean(tp, 1, 'omitnan')); % K x K
end
