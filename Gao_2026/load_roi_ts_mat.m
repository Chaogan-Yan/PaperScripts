function X = load_roi_ts_mat(matFile, roiN)
% Load a numeric matrix of size T x roiN (or roiN x T) from a .mat file.
% If multiple numeric matrices exist, it picks the one matching roiN.
% check data

S = load(matFile);
vars = fieldnames(S);

cand = {};
for i = 1:numel(vars)
    v = S.(vars{i});
    if isnumeric(v) && ismatrix(v)
        sz = size(v);
        if sz(2) == roiN || sz(1) == roiN
            cand{end+1} = v; %#ok<AGROW>
        end
    end
end

if isempty(cand)
    error('No numeric matrix with ROI=%d found in %s. Vars: %s', ...
        roiN, matFile, strjoin(vars, ', '));
end

% choose the first candidate to check data 
X = cand{1};
if size(X,2) == roiN
    % T x ROI OK
elseif size(X,1) == roiN
    X = X.'; % transpose to T x ROI
else
    error('Unexpected matrix shape after selection.');
end

% basic sanity
if any(~isfinite(X(:)))
    warning('Non-finite values found in %s. You may want to clean/scrub before CAP.', matFile);
end

end
