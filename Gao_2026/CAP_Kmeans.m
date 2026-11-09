function [clust,C,sumD,D,s,eva] = CAP_Kmeans(data,k_range,distance_func,replicate_num)

num_k = length(k_range);
clust = zeros(size(data,1), num_k);
s     = zeros(size(data,1), num_k);


opts = statset('UseParallel', true, ...
               'MaxIter', 1000, ...
               'Display', 'off');

C    = cell(1, num_k);
sumD = cell(1, num_k);
D    = cell(1, num_k);

for i = 1:num_k
    K = k_range(i);
    fprintf('Kmeans: %d clusters... \n', K);

    [clust(:,i), C{i}, sumD{i}, D{i}] = kmeans(data, K, ...
        'Distance',   distance_func, ...
        'Replicates', replicate_num, ...
        'Start',      'plus', ...
        'Options',    opts);

    
    s(:,i) = silhouette(data, clust(:,i), distance_func);
end


eva = evalclusters(data, clust, 'silhouette', 'Distance', distance_func);

end


