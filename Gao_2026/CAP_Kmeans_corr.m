function clust = CAP_Kmeans_corr(data, K, distance_func, CAP_HC_path)
% data: TP x ROI
% K: scalar
% distance_func: 'correlation' or 'euclidean' (see below)
% output clust: TP x 1 (labels 1..K)

TP_num = size(data, 1);

CAP_mat_path = fullfile(CAP_HC_path, ['CAP' num2str(K)], 'Group_CAP.mat');
S = load(CAP_mat_path, 'group_nii_clust_i_mean');
templates = S.group_nii_clust_i_mean.';  


switch lower(distance_func)
    case 'correlation'
        sim = 1 - pdist2(data, templates, 'correlation'); 
        [~, clust] = max(sim, [], 2);                     

    case 'euclidean'
        D = pdist2(data, templates, 'euclidean');         
        [~, clust] = min(D, [], 2);                       

    otherwise
        error('Unsupported distance_func: %s', distance_func);
end

clust = reshape(clust, TP_num, 1);
end