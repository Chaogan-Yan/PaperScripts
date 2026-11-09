function [CAP_matrices] = CAP_CalculateMatrices(clust, TP_sub, K, sub_num, varargin)
% CAP_CalculateMatrices (UPDATED with subnet-level spatial metrics)
%
%
% Optional Name-Value:
%   'roi_signal_all' : (T x nROI) double
%   'roiSubnetIdx'   : (nROI x 1) double, values 1..nSubnet
%   'subnetNames'    : cellstr (nSubnet x 1)

p = inputParser;
p.addParameter('roi_signal_all', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('roiSubnetIdx', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('subnetNames', {}, @(x) iscell(x) || isstring(x));
p.parse(varargin{:});

roi_signal_all = p.Results.roi_signal_all;
roiSubnetIdx   = p.Results.roiSubnetIdx;
subnetNames    = cellstr(string(p.Results.subnetNames));

doSubnet = ~isempty(roi_signal_all) && ~isempty(roiSubnetIdx);

% ---------------- original outputs (author) ----------------
fraction = zeros(K,sub_num);
persistence = zeros(K,sub_num);
counts = zeros(K,sub_num);
transitions = zeros(K,K,sub_num);
transitions_prob = zeros(K,K,sub_num);
resilience = zeros(K,sub_num);
in_degree = zeros(K,sub_num);
out_degree = zeros(K,sub_num);

% ---------------- subnet outputs (optional) ----------------
if doSubnet
    roiSubnetIdx = double(roiSubnetIdx(:));
    nROI = size(roi_signal_all,2);
    assert(numel(roiSubnetIdx)==nROI, 'roiSubnetIdx length must match roi_signal_all columns.');
    nSubnet = max(roiSubnetIdx(~isnan(roiSubnetIdx)));
    if isempty(subnetNames)
        subnetNames = arrayfun(@(i) sprintf('Subnet%d',i), 1:nSubnet, 'UniformOutput', false);
    end
    subnetNames = subnetNames(:);

    % ROI-z averaging within subnet (与你 heatmap 的做法最一致)
    subnet_z = nan(sub_num, K, nSubnet);

    % subnet-as-timeseries (可选更稳健)
    subnet_mean_ts = nan(sub_num, K, nSubnet);
    subnet_std_ts  = nan(sub_num, K, nSubnet);
    subnet_z_ts    = nan(sub_num, K, nSubnet);

    subnet_nFrames = zeros(sub_num, K);
end

% ---------------- main loop ----------------
for sub_i = 1:sub_num
    TP_sub_i_idx = find(TP_sub == sub_i);
    TP_sub_i_num = length(TP_sub_i_idx);
    clust_sub_i = clust(TP_sub_i_idx);

    clust_i_sequence_onset = cell(1, K);
    clust_i_sequence_end   = cell(1, K);
    clust_i_sequence_num   = cell(1, K);
    clust_i_sequence_duration = cell(1, K);

    for clust_i = 1:K
        clust_i_idx = find(clust_sub_i == clust_i);
        clust_i_num = length(clust_i_idx);

        clust_i_sequence = double(clust_sub_i == clust_i);
        clust_i_sequence_diff = diff([0; clust_i_sequence; 0]);

        clust_i_sequence_onset{clust_i} = find(clust_i_sequence_diff == 1);
        clust_i_sequence_end{clust_i} = find(clust_i_sequence_diff == -1) - 1;

        if length(clust_i_sequence_onset{clust_i}) ~= length(clust_i_sequence_end{clust_i})
            fprintf('Error: The length of onset and end are different !!! \n');
        end

        clust_i_sequence_num{clust_i} = length(clust_i_sequence_onset{clust_i});
        clust_i_sequence_duration{clust_i} = clust_i_sequence_end{clust_i} - clust_i_sequence_onset{clust_i} + 1;

        counts(clust_i,sub_i) = clust_i_sequence_num{clust_i};
        fraction(clust_i,sub_i) = clust_i_num/TP_sub_i_num;
        persistence(clust_i,sub_i) = mean(clust_i_sequence_duration{clust_i});
    end

    % transitions frequency
    for clust_i = 1:K
        seq_num = clust_i_sequence_num{clust_i};
        seq_end = clust_i_sequence_end{clust_i};
        for seq_i = 1:seq_num
            seq_end_i = seq_end(seq_i);
            if seq_end_i < TP_sub_i_num
                clust_next = clust_sub_i(seq_end_i + 1);
                transitions(clust_i,clust_next,sub_i) = transitions(clust_i,clust_next,sub_i) + 1;
            end
        end
    end

    % transitions probability + resilience
    end_clust = clust_sub_i(TP_sub_i_num);
    TP_num_clust = fraction(:,sub_i)*TP_sub_i_num;
    TP_num_clust(end_clust) = TP_num_clust(end_clust)-1;

    for clust_i = 1:K
        transitions_prob(clust_i,:,sub_i) = transitions(clust_i,:,sub_i) ./ TP_num_clust(clust_i);
        transitions_prob(clust_i,clust_i,sub_i) = 1 - sum(transitions_prob(clust_i,:,sub_i));
    end
    transitions_prob_sub_i = transitions_prob(:,:,sub_i);
    resilience(:,sub_i) = transitions_prob_sub_i(1:K+1:end);

    % ---------------- NEW: subnet metrics ----------------
    if doSubnet
        Xsub = roi_signal_all(TP_sub_i_idx, :); % TP_sub_i_num x nROI
        for clust_i = 1:K
            idxK = (clust_sub_i == clust_i);
            nF = sum(idxK);
            subnet_nFrames(sub_i, clust_i) = nF;

            if nF < 2
                continue;
            end

            frames = Xsub(idxK, :); % nF x nROI

            % ROI-level mean/std/z
            mu_roi = mean(frames, 1);
            sd_roi = std(frames, 0, 1);
            sd_roi(sd_roi==0) = eps;
            z_roi = mu_roi ./ sd_roi; % 1 x nROI

            nSubnet = size(subnet_z,3);
            for s = 1:nSubnet
                ridx = (roiSubnetIdx == s);
                if ~any(ridx), continue; end

                % (A) ROI z 平均（默认）
                subnet_z(sub_i, clust_i, s) = mean(z_roi(ridx), 'omitnan');

                % (B) 把子网当作“子网平均时间序列”，再算 mean/std/z
                ts = mean(frames(:,ridx), 2, 'omitnan');
                mu_ts = mean(ts, 'omitnan');
                sd_ts = std(ts, 0, 'omitnan'); if sd_ts==0, sd_ts=eps; end
                subnet_mean_ts(sub_i, clust_i, s) = mu_ts;
                subnet_std_ts(sub_i, clust_i, s)  = sd_ts;
                subnet_z_ts(sub_i, clust_i, s)    = mu_ts / sd_ts;
            end
        end
    end
end

in_degree = squeeze(sum(transitions,1));
out_degree = squeeze(sum(transitions,2));

% -------- keep original save orientation (same as your old code) :contentReference[oaicite:2]{index=2} --------
CAP_matrices.fraction = fraction';
CAP_matrices.persistence = persistence';
CAP_matrices.counts = counts';
CAP_matrices.transitions = permute(transitions,[3 1 2]);
CAP_matrices.transitions_prob = permute(transitions_prob,[3 1 2]);
CAP_matrices.resilience = resilience';
CAP_matrices.in_degree = in_degree';
CAP_matrices.out_degree = out_degree';

% -------- add new fields --------
if doSubnet
    CAP_matrices.subnet_z = subnet_z;                 % sub x K x nSubnet
    CAP_matrices.subnet_mean_ts = subnet_mean_ts;     % sub x K x nSubnet
    CAP_matrices.subnet_std_ts  = subnet_std_ts;      % sub x K x nSubnet
    CAP_matrices.subnet_z_ts    = subnet_z_ts;        % sub x K x nSubnet
    CAP_matrices.subnet_nFrames = subnet_nFrames;     % sub x K
    CAP_matrices.subnetNames    = subnetNames;        % nSubnet x 1
end
end
