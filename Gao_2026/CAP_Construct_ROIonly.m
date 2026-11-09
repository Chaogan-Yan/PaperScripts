function group_nii_clust_i_z = CAP_Construct_ROIonly(clust, K, roi_signal_all, outDir, varargin)
% CAP_Construct_ROIonly (UPDATED)

p = inputParser;
p.addParameter('roiSubnetIdx', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('subnetNames', {}, @(x) iscell(x) || isstring(x));
p.parse(varargin{:});
roiSubnetIdx = p.Results.roiSubnetIdx;
subnetNames  = cellstr(string(p.Results.subnetNames));

assert(size(roi_signal_all,1) == numel(clust), 'clust length must match roi_signal_all rows');

group_nii_clust_i = cell(1, K);
group_nii_clust_i_mean = nan(size(roi_signal_all,2), K);
group_nii_clust_i_std  = nan(size(roi_signal_all,2), K);
group_nii_clust_i_z    = nan(size(roi_signal_all,2), K);

for clust_i = 1:K
    clust_i_idx = find(clust == clust_i);

    if isempty(clust_i_idx)
        warning('CAP %d: empty cluster (no frames).', clust_i);
        group_nii_clust_i{clust_i} = [];
        continue;
    end

    group_nii_clust_i{clust_i} = roi_signal_all(clust_i_idx, :);

    mu = mean(group_nii_clust_i{clust_i}, 1);
    sd = std(group_nii_clust_i{clust_i}, 0, 1);
    sd(sd == 0) = eps;

    group_nii_clust_i_mean(:, clust_i) = mu(:);
    group_nii_clust_i_std(:, clust_i)  = sd(:);
    group_nii_clust_i_z(:, clust_i)    = (mu(:) ./ sd(:));
end

% -------- NEW: subnet-level group patterns --------
group_subnet_clust_i_z = [];
group_subnet_clust_i_mean_ts = [];
group_subnet_clust_i_std_ts  = [];
group_subnet_clust_i_z_ts    = [];

if ~isempty(roiSubnetIdx)
    roiSubnetIdx = double(roiSubnetIdx(:));
    nROI = size(roi_signal_all,2);
    assert(numel(roiSubnetIdx)==nROI, 'roiSubnetIdx length must match ROI count.');
    nSubnet = max(roiSubnetIdx(~isnan(roiSubnetIdx)));

    if isempty(subnetNames)
        subnetNames = arrayfun(@(i) sprintf('Subnet%d',i), 1:nSubnet, 'UniformOutput', false);
    end
    subnetNames = subnetNames(:);

    group_subnet_clust_i_z = nan(nSubnet, K);

    group_subnet_clust_i_mean_ts = nan(nSubnet, K);
    group_subnet_clust_i_std_ts  = nan(nSubnet, K);
    group_subnet_clust_i_z_ts    = nan(nSubnet, K);

    for clust_i = 1:K
        frames = group_nii_clust_i{clust_i};
        if isempty(frames), continue; end

        
        z_roi = group_nii_clust_i_z(:, clust_i); % 112x1
        for s = 1:nSubnet
            ridx = (roiSubnetIdx == s);
            if any(ridx)
                group_subnet_clust_i_z(s, clust_i) = mean(z_roi(ridx), 'omitnan');
            end
        end

        
        for s = 1:nSubnet
            ridx = (roiSubnetIdx == s);
            if ~any(ridx), continue; end
            ts = mean(frames(:,ridx), 2, 'omitnan');
            mu_ts = mean(ts, 'omitnan');
            sd_ts = std(ts, 0, 'omitnan'); if sd_ts==0, sd_ts=eps; end
            group_subnet_clust_i_mean_ts(s, clust_i) = mu_ts;
            group_subnet_clust_i_std_ts(s, clust_i)  = sd_ts;
            group_subnet_clust_i_z_ts(s, clust_i)    = mu_ts / sd_ts;
        end
    end
end

save(fullfile(outDir, 'Group_CAP.mat'), ...
    'group_nii_clust_i', 'group_nii_clust_i_mean', 'group_nii_clust_i_std', 'group_nii_clust_i_z', ...
    'group_subnet_clust_i_z', 'group_subnet_clust_i_mean_ts', 'group_subnet_clust_i_std_ts', 'group_subnet_clust_i_z_ts', ...
    'subnetNames');

end
