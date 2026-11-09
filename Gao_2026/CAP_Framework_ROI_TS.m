function CAP_Framework_ROI_TS(working_dir, dataDir, K_Means_Para)
% CAP pipeline for ROI time series files (T x 112), all frames into k-means.
% Input: dataDir with each subject's .mat file (T x 112 matrix)
% Output: CAP clustering results and dynamic matrices

if ~exist(working_dir, 'dir'); mkdir(working_dir); end

% --------- list subjects (.mat files) ----------
subjfiles = dir(fullfile(dataDir, '*.mat'));
assert(~isempty(subjfiles), 'No .mat files found under: %s', dataDir);
sub_num = numel(subjfiles);

fprintf('Found %d subjects (.mat) in %s\n', sub_num, dataDir);

% --------- load + zscore + concat ----------
roi_signal_all = [];
TP_sub = [];                      % subject label
TP_each = zeros(sub_num, 1);       % record T per subject
subj_names = cell(sub_num, 1);

for subji = 1:sub_num
    fpath = fullfile(subjfiles(subji).folder, subjfiles(subji).name);
    subj_names{subji} = subjfiles(subji).name;

    roi_signal = load_roi_ts_mat(fpath, 112);   % -> T x 112

    % zscore(roi_signal) in author pipeline (column-wise)
    roi_signal_zscore = zscore(roi_signal);     % T x 112

    roi_signal_all = [roi_signal_all; roi_signal_zscore]; %#ok<AGROW>

    T_i = size(roi_signal_zscore, 1);
    TP_each(subji) = T_i;
    TP_sub = [TP_sub; repmat(subji, T_i, 1)];  %#ok<AGROW>
end

save(fullfile(working_dir, 'roi_signal_all.mat'), 'roi_signal_all', 'TP_sub', 'TP_each', 'subj_names', '-v7.3');
fprintf('Saved roi_signal_all.mat (frames=%d, ROI=%d)\n', size(roi_signal_all,1), size(roi_signal_all,2));

% --------- K-means clustering ----------
k_range = K_Means_Para.Range;
k_num = K_Means_Para.Num;
distance_func = K_Means_Para.DistanceFunc;
replicate_num = K_Means_Para.ReplicateNum;

% Author pipeline calls: [clust,C,sumD,D,s,eva] = CAP_Kmeans(...)
[clust, C, sumD, D, s, eva] = CAP_Kmeans(roi_signal_all, k_range, distance_func, replicate_num);

save(fullfile(working_dir, 'clust.mat'), 'clust', 'C', 'sumD', 'D', 's', 'eva', 'k_range', 'distance_func', 'replicate_num');
fprintf('Saved clust.mat\n');

% --------- Construct group CAP patterns ----------
for K_idx = 1:k_num
    K = k_range(K_idx);
    clust_i = clust(:, K_idx);

    CAP_path.group = fullfile(working_dir, ['CAP' num2str(K)], filesep);
    if ~exist(CAP_path.group, 'dir'); mkdir(CAP_path.group); end

    % ROI-only construct: save Group_CAP.mat with same variable names
    CAP_Construct_ROIonly(clust_i, K, roi_signal_all, CAP_path.group);
end

% --------- Calculate CAP matrices (dynamic metrics) ----------
for k_i = 1:k_num
    K = k_range(k_i);
    clust_i = clust(:, k_i);

    % Author pipeline calls: CAP_matrices = CAP_CalculateMatrices(clust_i,TP_sub,K,sub_num);
    CAP_matrices = CAP_CalculateMatrices(clust_i, TP_sub, K, sub_num);

    % Save the CAP matrices
    save(fullfile(working_dir, ['CAP' num2str(K)], 'CAP_matrices.mat'), 'CAP_matrices');
end

fprintf('Done. Results saved under: %s\n', working_dir);
end

