%% run_CAP3_Recompute_SubnetMetrics.m
clear; clc;

%% ---------------- paths (your current structure) ----------------
rootEuc = '/Users/gaoqinglin/Downloads/DIRECT/AMN/CAP/Results/Correlation';


mddBase = fullfile(rootEuc, 'MDD');
mddCap3 = fullfile(mddBase, 'CAP3'); % contains CAP_matrices.mat, Group_CAP.mat, MDD_matched_clust.mat

ncBase  = fullfile(rootEuc, 'NC');
ncCap3  = fullfile(ncBase, 'CAP3');  % contains CAP_matrices.mat, Group_CAP.mat
ncClustFile = fullfile(ncBase, 'clust.mat'); % contains clust (T x 9), k_range=[2..10] etc.

labelFile = '/Users/gaoqinglin/Downloads/DIRECT/AMN/DMN/Results/settingLabel.mat';

Ktarget = 3;

%% ---------------- load labels ----------------
L = load(labelFile);
roiSubnetIdx = double(L.Node_Network(:)); % 112x1
subnetNames  = cellstr(strtrim(string(L.Network_Label(:,2)))); % 7x1
fprintf('[INFO] subnetNames:\n'); disp(subnetNames);

%% ---------------- MDD ----------------
fprintf('\n==== MDD: build TP_sub from clust_by_subj ====\n');
mddMatchedFile = fullfile(mddCap3, 'MDD_matched_clust.mat');
assert(exist(mddMatchedFile,'file')==2, 'Missing: %s', mddMatchedFile);

Sm = load(mddMatchedFile);
clust_mdd = double(Sm.clust(:));
assert(isfield(Sm,'clust_by_subj') && iscell(Sm.clust_by_subj), 'MDD_matched_clust.mat must contain clust_by_subj');
TP_sub_mdd = TP_sub_from_clust_by_subj(Sm.clust_by_subj);

assert(numel(TP_sub_mdd) == numel(clust_mdd), 'MDD: TP_sub length != clust length');
sub_num_mdd = numel(Sm.clust_by_subj);
K_mdd = max(clust_mdd);
assert(K_mdd == Ktarget, 'MDD: max(clust)=%d, expected K=%d', K_mdd, Ktarget);

[roi_mdd, subj_names_mdd] = load_roi_signal_all_any(mddBase, 'MDD');
assert(size(roi_mdd,1) == numel(clust_mdd), 'MDD: roi_signal_all rows != clust length');

CAP_matrices = CAP_CalculateMatrices(clust_mdd, TP_sub_mdd, Ktarget, sub_num_mdd, ...
    'roi_signal_all', roi_mdd, 'roiSubnetIdx', roiSubnetIdx, 'subnetNames', subnetNames);

backup_and_save(fullfile(mddCap3,'CAP_matrices.mat'), CAP_matrices);
CAP_Construct_ROIonly(clust_mdd, Ktarget, roi_mdd, mddCap3, 'roiSubnetIdx', roiSubnetIdx, 'subnetNames', subnetNames);
fprintf('[DONE] MDD updated.\n');

%% ---------------- NC ----------------
fprintf('\n==== NC: select clust column by k_range==3 ====\n');
assert(exist(ncClustFile,'file')==2, 'Missing: %s', ncClustFile);
Sn = load(ncClustFile);

clust_nc = pick_nc_clust_vector(Sn, Ktarget);

% try to find TP_sub in ncClustFile or roi_signal_all, else derive from cell roi_signal_all
TP_sub_nc = [];
if isfield(Sn,'TP_sub'), TP_sub_nc = Sn.TP_sub; end
if isempty(TP_sub_nc) && isfield(Sn,'TPsub'),  TP_sub_nc = Sn.TPsub;  end
if isempty(TP_sub_nc) && isfield(Sn,'tp_sub'), TP_sub_nc = Sn.tp_sub; end

[roi_nc, subj_names_nc, TP_sub_from_roi] = load_roi_signal_all_any(ncBase, 'NC');
if isempty(TP_sub_nc)
    TP_sub_nc = TP_sub_from_roi;
end

assert(~isempty(TP_sub_nc), ...
    ['NC: cannot find TP_sub. Fix by either:\n' ...
     '(1) put TP_sub into NC/clust.mat or NC/roi_signal_all.mat, OR\n' ...
     '(2) store roi_signal_all as cell array per subject (script will derive TP_sub).']);

TP_sub_nc = double(TP_sub_nc(:));
assert(numel(TP_sub_nc) == numel(clust_nc), 'NC: TP_sub length != clust length');
assert(size(roi_nc,1) == numel(clust_nc), 'NC: roi_signal_all rows != clust length');

sub_num_nc = max(TP_sub_nc);
K_nc = max(clust_nc);
assert(K_nc == Ktarget, 'NC: max(clust)=%d, expected K=%d', K_nc, Ktarget);

CAP_matrices = CAP_CalculateMatrices(clust_nc, TP_sub_nc, Ktarget, sub_num_nc, ...
    'roi_signal_all', roi_nc, 'roiSubnetIdx', roiSubnetIdx, 'subnetNames', subnetNames);

backup_and_save(fullfile(ncCap3,'CAP_matrices.mat'), CAP_matrices);
CAP_Construct_ROIonly(clust_nc, Ktarget, roi_nc, ncCap3, 'roiSubnetIdx', roiSubnetIdx, 'subnetNames', subnetNames);
fprintf('[DONE] NC updated.\n');

fprintf('\n[ALL DONE] Both groups CAP_matrices & Group_CAP have subnet-level outputs.\n');

%% ========================= helper functions =========================

function TP_sub = TP_sub_from_clust_by_subj(clust_by_subj)
nSub = numel(clust_by_subj);
lens = zeros(nSub,1);
for i=1:nSub
    lens(i) = numel(clust_by_subj{i});
end
TP_sub = zeros(sum(lens),1);
pos = 1;
for i=1:nSub
    TP_sub(pos:pos+lens(i)-1) = i;
    pos = pos + lens(i);
end
end



function clust_vec = pick_nc_clust_vector(Sn, Ktarget)
% NC/clust.mat screenshot shows:
% - clust: (T x 9) double
% - k_range: [2..10]

assert(isfield(Sn,'clust'), 'NC clust.mat missing variable clust');
C = Sn.clust;

if isvector(C)
    clust_vec = double(C(:));
    return;
end

% matrix case: T x nK
if isfield(Sn,'k_range')
    kr = double(Sn.k_range(:))';
    idx = find(kr == Ktarget, 1);
    if isempty(idx)
        error('NC: k_range exists but does not contain K=%d', Ktarget);
    end
    clust_vec = double(C(:,idx));
else
    % fallback: assume columns correspond to k_range = 2:10
    % then K=3 is column 2
    if size(C,2) >= 2 && Ktarget==3
        clust_vec = double(C(:,2));
    else
        error('NC: clust is matrix but k_range missing; cannot choose K column safely.');
    end
end

clust_vec = clust_vec(:);
end

function [roi_signal_all, subj_names, TP_sub] = load_roi_signal_all_any(baseDir, tag)
% try to load roi_signal_all from baseDir.
% If roi_signal_all is cell (nSub x 1), concatenate and derive TP_sub.
% If matrix, try to fetch TP_sub inside the file.
TP_sub = [];
subj_names = {};

cand = {fullfile(baseDir,'roi_signal_all.mat'), fullfile(baseDir,[lower(tag) '_roi_signal_all.mat'])};
roiFile = '';
for i=1:numel(cand)
    if exist(cand{i},'file')==2, roiFile = cand{i}; break; end
end
if isempty(roiFile)
    d = dir(fullfile(baseDir,'*roi*signal*all*.mat'));
    if ~isempty(d)
        roiFile = fullfile(baseDir, d(1).name);
    end
end
assert(~isempty(roiFile) && exist(roiFile,'file')==2, '[%s] cannot find roi_signal_all mat under %s', tag, baseDir);

S = load(roiFile);

roi_signal_all = pick_var(S, {'roi_signal_all','ROI_signal_all','MDD_roi_signal_all','NC_roi_signal_all','roi_all'});
assert(~isempty(roi_signal_all), '[%s] cannot find roi_signal_all variable in %s', tag, roiFile);

if isfield(S,'subj_names'), subj_names = S.subj_names; end

% find TP_sub if present
if isfield(S,'TP_sub'), TP_sub = S.TP_sub; end
if isempty(TP_sub) && isfield(S,'TPsub'), TP_sub = S.TPsub; end
if isempty(TP_sub) && isfield(S,'tp_sub'), TP_sub = S.tp_sub; end

% if cell: concat and derive TP_sub
if iscell(roi_signal_all)
    nSub = numel(roi_signal_all);
    lens = zeros(nSub,1);
    for i=1:nSub
        lens(i) = size(roi_signal_all{i},1);
    end
    TP_sub = zeros(sum(lens),1);
    X = zeros(sum(lens), size(roi_signal_all{1},2));
    pos = 1;
    for i=1:nSub
        li = lens(i);
        X(pos:pos+li-1,:) = roi_signal_all{i};
        TP_sub(pos:pos+li-1) = i;
        pos = pos + li;
    end
    roi_signal_all = X;
end

roi_signal_all = double(roi_signal_all);
if ~isempty(TP_sub), TP_sub = double(TP_sub(:)); end

fprintf('[%s] roi file: %s\n', tag, roiFile);
fprintf('[%s] roi_signal_all: %d x %d\n', tag, size(roi_signal_all,1), size(roi_signal_all,2));
if ~isempty(TP_sub)
    fprintf('[%s] TP_sub found/derived: nSub=%d\n', tag, max(TP_sub));
else
    fprintf('[%s] TP_sub not found/derived.\n', tag);
end
end

function v = pick_var(S, names)
v = [];
for i=1:numel(names)
    nm = names{i};
    if isfield(S,nm)
        v = S.(nm);
        return;
    end
end
end

function backup_and_save(matFile, CAP_matrices)
if exist(matFile,'file')==2
    bak = strrep(matFile,'.mat','_backup_ROIonly.mat');
    if exist(bak,'file')~=2
        copyfile(matFile, bak);
        fprintf('[INFO] backup created: %s\n', bak);
    end
end
save(matFile, 'CAP_matrices');
fprintf('[INFO] saved: %s\n', matFile);
end
