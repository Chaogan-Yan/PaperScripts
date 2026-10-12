function ana_DPABIViewerFixed(cfg)
%ANA_DPABIVIEWERFIXED Report GRF-corrected FC clusters for all study seeds.
%
% The function reads DPABI T2.nii result maps, applies the documented
% two-tailed GRF thresholds, identifies positive and negative clusters
% separately, and annotates peaks after nearest-neighbour atlas reslicing.
%
% Paths are supplied in cfg or through environment variables so the file can
% be published without workstation-specific information:
%   cfg.dataRoot    or MDD_FC_RESULTS_ROOT
%   cfg.templateDir or DPABI_TEMPLATE_DIR
%   cfg.seedManifest or MDD_FC_SEED_MANIFEST
%
% Example:
%   cfg = struct('dataRoot', '/data/fc_results', ...
%       'templateDir', '/opt/DPABI/Templates', ...
%       'seedManifest', 'config/fc_seed_manifest.csv');
%   ana_DPABIViewerFixed(cfg)

if nargin < 1 || isempty(cfg)
    cfg = struct();
end

clc;
fprintf('==== MDD_HippoAmy_GRF_Atlas_Report_Batch_DPABIViewerFixed start ====\n');

%% ========================================================================
% Public configuration
%% ========================================================================
cfg = cfg_default(cfg, 'dataRoot', getenv('MDD_FC_RESULTS_ROOT'));
cfg = cfg_default(cfg, 'templateDir', getenv('DPABI_TEMPLATE_DIR'));
cfg = cfg_default(cfg, 'seedManifest', getenv('MDD_FC_SEED_MANIFEST'));
cfg = cfg_default(cfg, 'expectedSeedCount', 4);
cfg = cfg_default(cfg, 'voxelP', 0.001);
cfg = cfg_default(cfg, 'clusterP', 0.05);
cfg = cfg_default(cfg, 'isTwoTailed', true);
cfg = cfg_default(cfg, 'connectivity', 26);
cfg = cfg_default(cfg, 'MANUAL_DF1', []);
cfg = cfg_default(cfg, 'MANUAL_DLH', []);
cfg = cfg_default(cfg, 'MANUAL_FWHM', []);
cfg = cfg_default(cfg, 'SMOOTHNESS_SOURCE', 'header_strict');
cfg = cfg_default(cfg, 'PositiveDirectionLabel', 'MDD>HC');
cfg = cfg_default(cfg, 'NegativeDirectionLabel', 'MDD<HC');

if isempty(cfg.dataRoot) || isempty(cfg.templateDir) || isempty(cfg.seedManifest)
    error(['Provide cfg.dataRoot, cfg.templateDir, and cfg.seedManifest, or set ' ...
        'MDD_FC_RESULTS_ROOT, DPABI_TEMPLATE_DIR, and MDD_FC_SEED_MANIFEST.']);
end
cfg = cfg_default(cfg, 'maskFile', ...
    fullfile(cfg.templateDir, 'BrainMask_05_97x115x97.nii'));
cfg = cfg_default(cfg, 'outputRoot', fullfile(cfg.dataRoot, ...
    sprintf('REPORT_GRF_voxP_%g_clusterP_%g_DPABIViewerFixed', ...
    cfg.voxelP, cfg.clusterP)));

% Atlas resources distributed with DPABI.
cfg.BrainnetomeNiiGz = fullfile(cfg.templateDir, 'BrainnetomeAtlas_BNA_MPM_thr25_1.25mm.nii.gz');
cfg.BrainnetomeXlsx  = fullfile(cfg.templateDir, 'BrainnetomeAtlas_BNA_subregions.xlsx');

cfg.AAL3Nii          = fullfile(cfg.templateDir, 'AAL3v1_1mm.nii');
cfg.AAL3Mat          = fullfile(cfg.templateDir, 'AAL3v1_1mm_Labels.mat');

cfg.Yeo17Nii         = fullfile(cfg.templateDir, 'Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii');

cfg.HOCortNii        = fullfile(cfg.templateDir, 'HarvardOxford-cort-maxprob-thr25-2mm_YCG.nii');
cfg.HOCortMat        = fullfile(cfg.templateDir, 'HarvardOxford-cort-maxprob-thr25-2mm_YCG_Labels.mat');
cfg.HOSubNii         = fullfile(cfg.templateDir, 'HarvardOxford-sub-maxprob-thr25-2mm_YCG.nii');
cfg.HOSubMat         = fullfile(cfg.templateDir, 'HarvardOxford-sub-maxprob-thr25-2mm_YCG_Labels.mat');

%% ========================================================================
% Dependency and input checks
%% ========================================================================
must_exist('y_ReadRPI');
must_exist('y_Write');
must_exist('y_Reslice');

if strcmpi(cfg.SMOOTHNESS_SOURCE, 'estimate') || strcmpi(cfg.SMOOTHNESS_SOURCE, 'estimate_if_missing')
    must_exist('y_Smoothest');
end

if ~(exist('bwconncomp', 'file') == 2 || exist('bwlabeln', 'file') == 2)
    error('bwconncomp/bwlabeln was not found. Confirm that Image Processing Toolbox is available.');
end

if ~exist(cfg.maskFile, 'file')
    error('Mask file does not exist: %s', cfg.maskFile);
end

if ~exist(cfg.outputRoot, 'dir')
    mkdir(cfg.outputRoot);
end

%% ========================================================================
% Seed definitions and statistic-map paths
%% ========================================================================
roiDefs = build_roi_defs(cfg);

% Preflight every statistic map so a run cannot silently complete with only
% a subset of the manuscript's four seeds.
for r = 1:numel(roiDefs)
    tFile = resolve_manifest_path(cfg.dataRoot, roiDefs(r).statFile);
    if ~exist(tFile, 'file')
        error('Statistic map does not exist for seed %s: %s', ...
            roiDefs(r).roiName, tFile);
    end
end

%% ========================================================================
% Analysis mask
%% ========================================================================
[maskData, ~, ~] = y_ReadRPI(cfg.maskFile);
maskData = squeeze(maskData);
maskBin = maskData ~= 0;
maskSize = size(maskBin);

fprintf('Mask: %s\n', cfg.maskFile);
fprintf('Mask size: [%d %d %d]\n', maskSize(1), maskSize(2), maskSize(3));

%% ========================================================================
% Atlas definitions; each atlas is resliced to the current statistic map.
%% ========================================================================
fprintf('\n==== Loading atlas definitions ====' ); fprintf('\n');
atlases = load_all_atlases(cfg);
fprintf('Atlas definitions loaded: %d atlases.\n', numel(atlases));

%% ========================================================================
% Output tables
%% ========================================================================
allDetailed = empty_detailed_table();
allReadable = empty_readable_table();
roiSummary  = empty_roi_summary_table();
failedSeeds = {};

combinedByRoiXlsx = fullfile(cfg.outputRoot, 'All_ROIs_ByROI_ClusterSheets.xlsx');
if exist(combinedByRoiXlsx, 'file')
    delete(combinedByRoiXlsx);
end

%% ========================================================================
% Seed loop
%% ========================================================================
for r = 1:numel(roiDefs)
    def = roiDefs(r);

    tFile = resolve_manifest_path(cfg.dataRoot, def.statFile);
    roiOutDir = fullfile(cfg.outputRoot, def.roiName);

    fprintf('\n============================================================\n');
    fprintf('Processing %s\n', def.roiName);
    fprintf('Input: %s\n', tFile);

    if ~exist(roiOutDir, 'dir')
        mkdir(roiOutDir);
    end

    try
        [clusterTable, readableTable, runInfo] = process_one_roi( ...
            tFile, roiOutDir, def, cfg, maskBin, maskSize, atlases);

        % Per-seed outputs
        detailCsv   = fullfile(roiOutDir, sprintf('%s_clusters_detailed.csv', def.roiName));
        readableCsv = fullfile(roiOutDir, sprintf('%s_clusters_report_ready.csv', def.roiName));
        writetable(clusterTable, detailCsv);
        writetable(readableTable, readableCsv);

        % Plain-text summary
        txtFile = fullfile(roiOutDir, sprintf('%s_report_summary.txt', def.roiName));
        write_readable_text_summary(txtFile, def.roiName, readableTable, runInfo);

        % Append to combined outputs
        allDetailed = append_table_like(allDetailed, clusterTable);
        allReadable = append_table_like(allReadable, readableTable);

        signCol = get_text_col_as_cellstr(clusterTable.Sign);
        nPos = sum(strcmp(signCol, 'Positive'));
        nNeg = sum(strcmp(signCol, 'Negative'));

        roiSummary = append_table_like(roiSummary, ...
            make_roi_summary_row(def, height(clusterTable), nPos, nNeg, 'OK'));

        % One worksheet per seed, with one row per significant cluster.
        byRoiSheetTable = make_byroi_cluster_sheet_table(clusterTable);
        write_byroi_sheet(combinedByRoiXlsx, r, def, byRoiSheetTable);

    catch ME
        warning('Error while processing %s: %s', def.roiName, ME.message);
        failedSeeds{end+1} = def.roiName; %#ok<AGROW>
        roiSummary = append_table_like(roiSummary, ...
            make_roi_summary_row(def, 0, 0, 0, ['ERROR: ' ME.message]));
        write_byroi_sheet(combinedByRoiXlsx, r, def, empty_byroi_cluster_sheet_table());
    end
end

%% ========================================================================
% Combined outputs
%% ========================================================================
combinedXlsx = fullfile(cfg.outputRoot, 'All_ROIs_GRF_Atlas_Report.xlsx');
combinedDetailedCsv = fullfile(cfg.outputRoot, 'All_ROIs_GRF_Atlas_Detailed.csv');
combinedReadableCsv = fullfile(cfg.outputRoot, 'All_ROIs_GRF_ReportReady.csv');
combinedRoiSummaryCsv = fullfile(cfg.outputRoot, 'All_ROIs_Summary.csv');

if exist(combinedXlsx, 'file')
    delete(combinedXlsx);
end

writetable(allDetailed, combinedXlsx, 'Sheet', 'Detailed');
writetable(allReadable, combinedXlsx, 'Sheet', 'ReportReady');
writetable(roiSummary, combinedXlsx, 'Sheet', 'ROISummary');

writetable(allDetailed, combinedDetailedCsv);
writetable(allReadable, combinedReadableCsv);
writetable(roiSummary, combinedRoiSummaryCsv);

% Combined text summary
allTxt = fullfile(cfg.outputRoot, 'All_ROIs_Report_Summary.txt');
write_overall_text_summary(allTxt, allReadable, roiSummary, cfg);

if ~isempty(failedSeeds)
    error('FC peak reporting failed for seed(s): %s. Review the per-seed summary.', ...
        strjoin(failedSeeds, ', '));
end

fprintf('\n==== All done ====' ); fprintf('\n');
fprintf('Output root: %s\n', cfg.outputRoot);
fprintf('Combined workbook: %s\n', combinedXlsx);
fprintf('By-ROI workbook: %s\n', combinedByRoiXlsx);

end

%% ========================================================================
function roiDefs = build_roi_defs(cfg)

roiDefs = struct( ...
    'statFile', {}, ...
    'labelID', {}, ...
    'rawName', {}, ...
    'roiName', {}, ...
    'hemisphere', {});

if ~exist(cfg.seedManifest, 'file')
    error('Seed manifest does not exist: %s', cfg.seedManifest);
end

manifest = readtable(cfg.seedManifest, 'TextType', 'string', ...
    'PreserveVariableNames', true);
required = ["SeedName", "StatFile", "Hemisphere", "LabelID", "RawName"];
missing = setdiff(required, string(manifest.Properties.VariableNames));
if ~isempty(missing)
    error('Seed manifest is missing columns: %s', strjoin(missing, ', '));
end
if height(manifest) ~= cfg.expectedSeedCount
    error('Expected %d seed rows but found %d in %s.', ...
        cfg.expectedSeedCount, height(manifest), cfg.seedManifest);
end
textFields = ["SeedName", "StatFile", "Hemisphere", "RawName"];
for field = textFields
    fieldName = char(field);
    values = strip(string(manifest.(fieldName)));
    if any(ismissing(values) | strlength(values) == 0)
        error('%s contains an empty value in %s.', cfg.seedManifest, fieldName);
    end
end
if numel(unique(manifest.SeedName)) ~= height(manifest)
    error('SeedName values must be unique in %s.', cfg.seedManifest);
end
if numel(unique(manifest.StatFile)) ~= height(manifest)
    error('StatFile values must be unique in %s.', cfg.seedManifest);
end

for i = 1:height(manifest)
    labelID = double(manifest.LabelID(i));
    if ~(isscalar(labelID) && isfinite(labelID) && labelID > 0 && labelID == fix(labelID))
        error('LabelID must be a positive integer for seed row %d.', i);
    end
    hemisphere = lower(strip(string(manifest.Hemisphere(i))));
    if ~ismember(hemisphere, ["l", "r", "left", "right"])
        error('Hemisphere must be L, R, Left, or Right for seed row %d.', i);
    end
    roiDefs(end+1).statFile   = char(manifest.StatFile(i)); %#ok<AGROW>
    roiDefs(end).labelID      = labelID;
    roiDefs(end).rawName      = char(manifest.RawName(i));
    roiDefs(end).hemisphere   = char(manifest.Hemisphere(i));
    roiDefs(end).roiName      = char(manifest.SeedName(i));
end

end

%% ========================================================================
function pathValue = resolve_manifest_path(dataRoot, suppliedPath)

suppliedPath = char(suppliedPath);
isWindowsAbsolute = numel(suppliedPath) >= 3 && isletter(suppliedPath(1)) ...
    && suppliedPath(2) == ':' && any(suppliedPath(3) == ['\' '/']);
isUnixAbsolute = ~isempty(suppliedPath) && suppliedPath(1) == '/';
if isWindowsAbsolute || isUnixAbsolute
    pathValue = suppliedPath;
else
    pathValue = fullfile(dataRoot, suppliedPath);
end

end

%% ========================================================================
function [clusterTable, readableTable, runInfo] = process_one_roi( ...
    tFile, roiOutDir, def, cfg, maskBin, maskSize, atlases)

% Statistic map
[tData, voxelSize, header] = y_ReadRPI(tFile);
tData = squeeze(double(tData));

if ndims(tData) ~= 3
    error('%s is not a 3-D volume.', tFile);
end

if ~isequal(size(tData), maskSize)
    error('%s has size [%d %d %d], which differs from mask size [%d %d %d].', ...
        tFile, size(tData,1), size(tData,2), size(tData,3), ...
        maskSize(1), maskSize(2), maskSize(3));
end

fprintf('T size: [%d %d %d]\n', size(tData,1), size(tData,2), size(tData,3));
fprintf('Voxel size: [%g %g %g]\n', voxelSize(1), voxelSize(2), voxelSize(3));

% Statistic type and degrees of freedom
[statInfo, headerExt] = get_stat_info_from_header(header, cfg);
fprintf('TestFlag = %s\n', statInfo.TestFlag);
fprintf('Df1 = %g\n', statInfo.Df1);
if ~isempty(statInfo.Df2)
    fprintf('Df2 = %g\n', statInfo.Df2);
end

% T-to-Z conversion is used only for optional smoothness estimation and reporting.
zData = stat_to_z_signed(tData, statInfo);
zData(~maskBin) = 0;

zBeforeFile = fullfile(roiOutDir, 'Z_BeforeThreshold_T2.nii');
zHeader = header;
zHeader.pinfo = [1; 0; 0];
zHeader.dt = [16, 0];
y_Write(single(zData), zHeader, zBeforeFile);
fprintf('Saved: %s\n', zBeforeFile);

% Smoothness metadata
[dLh, FWHM, smoothnessSource] = get_smoothness_dpabi_viewer_style( ...
    headerExt, cfg, zBeforeFile, cfg.maskFile, voxelSize);

fprintf('Smoothness source = %s\n', smoothnessSource);
fprintf('dLh = %.8f\n', dLh);
fprintf('FWHM = [%.4f %.4f %.4f]\n', FWHM(1), FWHM(2), FWHM(3));

% DPABI Viewer-style GRF cluster extent and statistic threshold
nMaskVox = nnz(maskBin);
[statThrd, zThrd, minClusterSize, Em, EN, Beta, clusterPUsed] = ...
    calc_dpabi_viewer_grf_threshold(nMaskVox, dLh, cfg.voxelP, cfg.clusterP, ...
    cfg.isTwoTailed, statInfo);

fprintf('Z threshold = %.6f\n', zThrd);
fprintf('Original %s threshold = %.6f\n', upper(statInfo.TestFlag), statThrd);
fprintf('Minimum cluster size = %d voxels\n', minClusterSize);
fprintf('nMaskVox = %d\n', nMaskVox);

% Apply voxel and cluster thresholds to the original statistic map.
[statCorr, keepMask] = apply_dpabi_viewer_grf_to_stat( ...
    tData, maskBin, statThrd, minClusterSize, cfg.isTwoTailed, cfg.connectivity);
zCorr = zData .* double(keepMask);

zCorrFile = fullfile(roiOutDir, 'Z_ClusterThresholded_T2.nii');
tCorrFile = fullfile(roiOutDir, 'ClusterThresholded_T2.nii');

y_Write(single(zCorr), zHeader, zCorrFile);

tHeader = header;
tHeader.pinfo = [1; 0; 0];
tHeader.dt = [16, 0];
y_Write(single(statCorr), tHeader, tCorrFile);

fprintf('Saved: %s\n', zCorrFile);
fprintf('Saved: %s\n', tCorrFile);

% Reslice atlases to the current T2.nii grid using nearest-neighbour interpolation.
reslicedAtlases = reslice_all_atlases_to_current_stat(atlases, tFile, voxelSize, size(tData));

% Detailed cluster table
clusterTable = build_cluster_table_with_atlas( ...
    def, statCorr, zCorr, voxelSize, header, statThrd, zThrd, cfg.connectivity, ...
    reslicedAtlases, cfg, statInfo, minClusterSize, dLh, FWHM, Em, EN, Beta, clusterPUsed, smoothnessSource);

% Report-ready cluster table
readableTable = make_readable_table(clusterTable);

% Run metadata
runInfo = struct();
runInfo.Df1 = statInfo.Df1;
runInfo.Df2 = statInfo.Df2;
runInfo.TestFlag = statInfo.TestFlag;
runInfo.VoxelP = cfg.voxelP;
runInfo.ClusterP = cfg.clusterP;
runInfo.ClusterPUsed = clusterPUsed;
runInfo.StatThreshold = statThrd;
runInfo.ZThreshold = zThrd;
runInfo.MinClusterSize = minClusterSize;
runInfo.dLh = dLh;
runInfo.FWHM = FWHM;
runInfo.Em = Em;
runInfo.EN = EN;
runInfo.Beta = Beta;
runInfo.SmoothnessSource = smoothnessSource;

end

%% ========================================================================
function [statInfo, headerExt] = get_stat_info_from_header(header, cfg)

headerExt = header;

if exist('w_ReadDF', 'file') == 2
    try
        headerExt = w_ReadDF(headerExt);
    catch ME
        warning('MDD:DPABIViewer:ReadDFFailed', ...
            'w_ReadDF failed; trying existing header fields or MANUAL_DF1: %s', ME.message);
    end
end

statInfo = struct();
statInfo.TestFlag = 'T';
statInfo.Df1 = [];
statInfo.Df2 = [];

if isfield(headerExt, 'TestFlag') && ~isempty(headerExt.TestFlag)
    statInfo.TestFlag = upper(strtrim(char(headerExt.TestFlag)));
end

if ~isempty(cfg.MANUAL_DF1)
    statInfo.Df1 = cfg.MANUAL_DF1;
elseif isfield(headerExt, 'Df') && ~isempty(headerExt.Df)
    statInfo.Df1 = double(headerExt.Df);
elseif isfield(headerExt, 'DF') && ~isempty(headerExt.DF)
    statInfo.Df1 = double(headerExt.DF);
end

if isfield(headerExt, 'Df2') && ~isempty(headerExt.Df2)
    statInfo.Df2 = double(headerExt.Df2);
elseif isfield(headerExt, 'DF2') && ~isempty(headerExt.DF2)
    statInfo.Df2 = double(headerExt.DF2);
end

if strcmpi(statInfo.TestFlag, 'T') || strcmpi(statInfo.TestFlag, 'R') || strcmpi(statInfo.TestFlag, 'F')
    if isempty(statInfo.Df1) || ~isscalar(statInfo.Df1) || ~isfinite(statInfo.Df1) || statInfo.Df1 <= 0
        error(['Df1 could not be read from the T2.nii header and cfg.MANUAL_DF1 is empty.' newline ...
               'Set cfg.MANUAL_DF1 in the configuration section.']);
    end
end

if strcmpi(statInfo.TestFlag, 'F')
    if isempty(statInfo.Df2) || ~isscalar(statInfo.Df2) || ~isfinite(statInfo.Df2) || statInfo.Df2 <= 0
        error('An F map requires Df2, but no valid value is available in the header or configuration.');
    end
end

end

%% ========================================================================
function [dLh, FWHM, sourceText] = get_smoothness_dpabi_viewer_style( ...
    headerExt, cfg, zBeforeFile, maskFile, voxelSize)

sourceText = '';
dLh = [];
FWHM = [];

if ~isempty(cfg.MANUAL_DLH)
    dLh = double(cfg.MANUAL_DLH);
    sourceText = 'manual dLh';
end

if ~isempty(cfg.MANUAL_FWHM)
    FWHM = double(cfg.MANUAL_FWHM(:))';
    if numel(FWHM) ~= 3
        error('cfg.MANUAL_FWHM must be a 1-by-3 array.');
    end
    if isempty(dLh)
        dLh = 1 / prod(FWHM);
        sourceText = 'manual FWHM -> dLh=1/prod(FWHM)';
    end
end

if isempty(dLh) || isempty(FWHM)
    if exist('w_ReadFWHM', 'file') == 2
        try
            headerExt = w_ReadFWHM(headerExt);
        catch ME
            warning('MDD:DPABIViewer:ReadFWHMFailed', ...
                'w_ReadFWHM failed: %s', ME.message);
        end
    end
end

if isempty(dLh)
    if isfield(headerExt, 'dLh') && ~isempty(headerExt.dLh)
        dLh = double(headerExt.dLh);
        sourceText = 'header dLh via w_ReadFWHM';
    elseif isfield(headerExt, 'DLH') && ~isempty(headerExt.DLH)
        dLh = double(headerExt.DLH);
        sourceText = 'header DLH';
    end
end

if isempty(FWHM)
    if isfield(headerExt, 'FWHM') && numel(headerExt.FWHM) >= 3
        FWHM = double(headerExt.FWHM(1:3));
    elseif isfield(headerExt, 'FWHMx') && isfield(headerExt, 'FWHMy') && isfield(headerExt, 'FWHMz')
        FWHM = [double(headerExt.FWHMx), double(headerExt.FWHMy), double(headerExt.FWHMz)];
    end
end

if isempty(dLh) && ~isempty(FWHM) && all(isfinite(FWHM)) && all(FWHM > 0)
    dLh = 1 / prod(FWHM);
    if isempty(sourceText)
        sourceText = 'header FWHM -> dLh=1/prod(FWHM)';
    end
end

needEstimate = false;
if strcmpi(cfg.SMOOTHNESS_SOURCE, 'estimate')
    needEstimate = true;
elseif strcmpi(cfg.SMOOTHNESS_SOURCE, 'estimate_if_missing')
    needEstimate = isempty(dLh) || isempty(FWHM) || ~isfinite(dLh) || dLh <= 0;
end

if needEstimate
    [dLh2, ~, FWHM2, ~] = y_Smoothest(zBeforeFile, maskFile, '', voxelSize);
    dLh = double(dLh2);
    FWHM = double(FWHM2(:))';
    sourceText = 'estimated by y_Smoothest on Z map';
end

if isempty(dLh) || ~isscalar(dLh) || ~isfinite(dLh) || dLh <= 0
    error(['No valid dLh is available. In header_strict mode, smoothness is not re-estimated.' newline ...
           'Confirm that the T2.nii header contains dLh/FWHM, set cfg.MANUAL_DLH, or use estimate_if_missing.']);
end

if isempty(FWHM) || numel(FWHM) < 3 || any(~isfinite(FWHM(1:3))) || any(FWHM(1:3) <= 0)
    FWHM = [NaN NaN NaN];
end
FWHM = FWHM(1:3);

if isempty(sourceText)
    sourceText = 'header';
end

end

%% ========================================================================
function [statThrd, zThrd, minClusterSize, Em, EN, Beta, clusterPUsed] = ...
    calc_dpabi_viewer_grf_threshold(nVoxels, dLh, voxelP, clusterP, isTwoTailed, statInfo)

if isTwoTailed
    zThrd = norminv_local(1 - voxelP/2);
    clusterPUsed = clusterP/2;
    denominator = 2;
else
    zThrd = norminv_local(1 - voxelP);
    clusterPUsed = clusterP;
    denominator = 1;
end

[minClusterSize, Em, EN, Beta] = calc_min_cluster_size(nVoxels, dLh, zThrd, clusterPUsed);

switch upper(statInfo.TestFlag)
    case 'T'
        statThrd = tinv_safe(1 - voxelP/denominator, statInfo.Df1);
    case 'R'
        T = tinv_safe(1 - voxelP/denominator, statInfo.Df1);
        statThrd = sqrt(T^2/(statInfo.Df1 + T^2));
    case 'F'
        statThrd = finv_safe(1 - voxelP/denominator, statInfo.Df1, statInfo.Df2);
    case 'Z'
        statThrd = norminv_local(1 - voxelP/denominator);
    otherwise
        error('Unsupported TestFlag: %s', statInfo.TestFlag);
end

end

%% ========================================================================
function [statCorr, keepMask] = apply_dpabi_viewer_grf_to_stat( ...
    statData, maskBin, statThrd, minClusterSize, isTwoTailed, connectivity)

% Threshold positive and negative tails independently. Combining both signs
% before connected-component labelling can incorrectly merge touching clusters
% of opposite direction.
if isTwoTailed
    positiveMask = (statData >= statThrd) & maskBin;
    negativeMask = (statData <= -statThrd) & maskBin;
    positiveKeep = keep_large_clusters(positiveMask, minClusterSize, connectivity);
    negativeKeep = keep_large_clusters(negativeMask, minClusterSize, connectivity);
    keepMask = positiveKeep | negativeKeep;
else
    positiveMask = (statData >= statThrd) & maskBin;
    keepMask = keep_large_clusters(positiveMask, minClusterSize, connectivity);
end

statCorr = zeros(size(statData));
statCorr(keepMask) = statData(keepMask);

end

%% ========================================================================
function atlases = load_all_atlases(cfg)

atlases = repmat(make_empty_atlas_struct(), 0, 1);

bnNii = ensure_nifti_file(cfg.BrainnetomeNiiGz);
if ~exist(bnNii, 'file')
    error('Brainnetome NIfTI does not exist: %s', bnNii);
end
if ~exist(cfg.BrainnetomeXlsx, 'file')
    error('Brainnetome label table does not exist: %s', cfg.BrainnetomeXlsx);
end
atlases(end+1) = make_atlas_struct('Brainnetome', 'BN', bnNii, load_brainnetome_label_map(cfg.BrainnetomeXlsx));

if ~exist(cfg.AAL3Nii, 'file')
    error('AAL3 NIfTI does not exist: %s', cfg.AAL3Nii);
end
if ~exist(cfg.AAL3Mat, 'file')
    error('AAL3 label MAT file does not exist: %s', cfg.AAL3Mat);
end
atlases(end+1) = make_atlas_struct('AAL3', 'AAL3', cfg.AAL3Nii, load_mat_label_map(cfg.AAL3Mat));

if ~exist(cfg.Yeo17Nii, 'file')
    error('Yeo17 NIfTI does not exist: %s', cfg.Yeo17Nii);
end
atlases(end+1) = make_atlas_struct('Yeo17', 'Yeo17', cfg.Yeo17Nii, build_yeo17_label_map());

if ~exist(cfg.HOCortNii, 'file')
    error('HarvardOxford cortical NIfTI does not exist: %s', cfg.HOCortNii);
end
if ~exist(cfg.HOCortMat, 'file')
    error('HarvardOxford cortical Labels.mat does not exist: %s', cfg.HOCortMat);
end
atlases(end+1) = make_atlas_struct('HarvardOxford_Cortical', 'HO_Cort', cfg.HOCortNii, load_mat_label_map(cfg.HOCortMat));

if ~exist(cfg.HOSubNii, 'file')
    error('HarvardOxford subcortical NIfTI does not exist: %s', cfg.HOSubNii);
end
if ~exist(cfg.HOSubMat, 'file')
    error('HarvardOxford subcortical Labels.mat does not exist: %s', cfg.HOSubMat);
end
atlases(end+1) = make_atlas_struct('HarvardOxford_Subcortical', 'HO_Sub', cfg.HOSubNii, load_mat_label_map(cfg.HOSubMat));

% Print label-map counts as a lightweight atlas-loading audit.
for a = 1:numel(atlases)
    fprintf('Atlas %d: %s | %s | labels=%d\n', a, atlases(a).Alias, atlases(a).File, atlases(a).LabelMap.Count);
end

end

%% ========================================================================
function atlas = make_empty_atlas_struct()
atlas = struct('Alias', '', 'ShortName', '', 'File', '', 'LabelMap', [], 'Data', []);
end

%% ========================================================================
function atlas = make_atlas_struct(alias, shortName, niiFile, labelMap)
atlas = make_empty_atlas_struct();
atlas.Alias = alias;
atlas.ShortName = shortName;
atlas.File = niiFile;
atlas.LabelMap = labelMap;
atlas.Data = [];
end

%% ========================================================================
function reslicedAtlases = reslice_all_atlases_to_current_stat(atlases, refFile, refVox, refSize)

reslicedAtlases = atlases;

for a = 1:numel(atlases)
    fprintf('Reslicing atlas to current T2 space: %s\n', atlases(a).Alias);
    data = y_Reslice(atlases(a).File, '', refVox, 0, refFile);
    data = round(squeeze(double(data)));

    if ~isequal(size(data), refSize)
        error('Resliced atlas %s has size [%s], which differs from T2 size [%s].', ...
            atlases(a).Alias, mat2str(size(data)), mat2str(refSize));
    end

    reslicedAtlases(a).Data = data;
end

end

%% ========================================================================
function clusterTable = build_cluster_table_with_atlas( ...
    def, statCorr, zCorr, voxelSize, header, statThrd, zThrd, connectivity, ...
    atlases, cfg, statInfo, minClusterSize, dLh, FWHM, Em, EN, Beta, clusterPUsed, smoothnessSource)

rows = repmat(make_empty_detail_struct_row(), 0, 1);

if ~any(statCorr(:))
    clusterTable = empty_detailed_table();
    return;
end

% Keep opposite effect directions as distinct clusters even if their voxels
% touch under the requested neighbourhood definition.
positiveLists = connected_pixel_lists(statCorr > 0, connectivity);
negativeLists = connected_pixel_lists(statCorr < 0, connectivity);
pixelLists = [positiveLists(:); negativeLists(:)];
numClusters = numel(pixelLists);

for c = 1:numClusters
    idx = pixelLists{c};
    if isempty(idx)
        continue;
    end

    clusterSizeVoxel = numel(idx);
    clusterSizeMM3 = clusterSizeVoxel * prod(abs(double(voxelSize(:)')));

    clusterT = statCorr(idx);
    clusterZ = zCorr(idx);

    [~, loc] = max(abs(clusterT));
    peakIdx = idx(loc);
    peakT = statCorr(peakIdx);
    peakZ = zCorr(peakIdx);

    if peakT >= 0
        signName = 'Positive';
        directionLabel = cfg.PositiveDirectionLabel;
    else
        signName = 'Negative';
        directionLabel = cfg.NegativeDirectionLabel;
    end

    [i, j, k] = ind2sub(size(statCorr), peakIdx);
    peakXYZmm = voxel_to_mm(header, [i, j, k]);

    meanT = mean(clusterT);
    meanZ = mean(clusterZ);

    atlasResults = get_all_atlas_results(atlases, idx, [i, j, k], clusterSizeVoxel);

    bn   = atlas_result_by_shortname(atlasResults, 'BN');
    aal  = atlas_result_by_shortname(atlasResults, 'AAL3');
    yeo  = atlas_result_by_shortname(atlasResults, 'Yeo17');
    hoc  = atlas_result_by_shortname(atlasResults, 'HO_Cort');
    hos  = atlas_result_by_shortname(atlasResults, 'HO_Sub');

    needBoundaryCheck = 0;
    if (bn.PeakID ~= 0 && bn.DominantID ~= 0 && bn.PeakID ~= bn.DominantID) || ...
       (bn.DominantPct > 0 && bn.DominantPct < 60) || ...
       (bn.CoveragePct > 0 && bn.CoveragePct < 70)
        needBoundaryCheck = 1;
    end

    [reportTitle, reportSentence] = build_report_text( ...
        def.roiName, directionLabel, c, clusterSizeVoxel, clusterSizeMM3, ...
        peakXYZmm, peakT, peakZ, bn, aal, yeo, hoc, hos, needBoundaryCheck);

    row = make_empty_detail_struct_row();
    row.ROI              = def.roiName;
    row.Hemisphere       = def.hemisphere;
    row.SeedLabelID      = def.labelID;
    row.SeedRawName      = def.rawName;

    row.Sign             = signName;
    row.Direction        = directionLabel;
    row.ClusterID        = c;
    row.ClusterSizeVoxel = clusterSizeVoxel;
    row.ClusterSizeMM3   = clusterSizeMM3;

    row.PeakT            = peakT;
    row.PeakZ            = peakZ;
    row.MeanT            = meanT;
    row.MeanZ            = meanZ;

    row.PeakI            = i;
    row.PeakJ            = j;
    row.PeakK            = k;
    row.PeakXmm          = peakXYZmm(1);
    row.PeakYmm          = peakXYZmm(2);
    row.PeakZmm          = peakXYZmm(3);

    row.Peak_BN_ID        = bn.PeakID;
    row.Peak_BN_Label     = bn.PeakLabel;
    row.Dominant_BN_ID    = bn.DominantID;
    row.Dominant_BN_Label = bn.DominantLabel;
    row.BN_CoveragePct    = bn.CoveragePct;
    row.BN_DominantPct    = bn.DominantPct;
    row.BN_Top3           = bn.Top3Text;

    row.Peak_AAL3_ID        = aal.PeakID;
    row.Peak_AAL3_Label     = aal.PeakLabel;
    row.Dominant_AAL3_ID    = aal.DominantID;
    row.Dominant_AAL3_Label = aal.DominantLabel;
    row.AAL3_CoveragePct    = aal.CoveragePct;
    row.AAL3_DominantPct    = aal.DominantPct;
    row.AAL3_Top3           = aal.Top3Text;

    row.Peak_Yeo17_ID        = yeo.PeakID;
    row.Peak_Yeo17_Label     = yeo.PeakLabel;
    row.Dominant_Yeo17_ID    = yeo.DominantID;
    row.Dominant_Yeo17_Label = yeo.DominantLabel;
    row.Yeo17_CoveragePct    = yeo.CoveragePct;
    row.Yeo17_DominantPct    = yeo.DominantPct;
    row.Yeo17_Top3           = yeo.Top3Text;

    row.Peak_HO_Cort_ID        = hoc.PeakID;
    row.Peak_HO_Cort_Label     = hoc.PeakLabel;
    row.Dominant_HO_Cort_ID    = hoc.DominantID;
    row.Dominant_HO_Cort_Label = hoc.DominantLabel;
    row.HO_Cort_CoveragePct    = hoc.CoveragePct;
    row.HO_Cort_DominantPct    = hoc.DominantPct;
    row.HO_Cort_Top3           = hoc.Top3Text;

    row.Peak_HO_Sub_ID        = hos.PeakID;
    row.Peak_HO_Sub_Label     = hos.PeakLabel;
    row.Dominant_HO_Sub_ID    = hos.DominantID;
    row.Dominant_HO_Sub_Label = hos.DominantLabel;
    row.HO_Sub_CoveragePct    = hos.CoveragePct;
    row.HO_Sub_DominantPct    = hos.DominantPct;
    row.HO_Sub_Top3           = hos.Top3Text;

    row.NeedBoundaryCheck = needBoundaryCheck;
    row.ReportTitle       = reportTitle;
    row.ReportSentence    = reportSentence;

    row.VoxelP         = cfg.voxelP;
    row.ClusterP       = cfg.clusterP;
    row.ClusterPUsed   = clusterPUsed;
    row.IsTwoTailed    = double(cfg.isTwoTailed);
    row.TestFlag       = statInfo.TestFlag;
    row.Df1            = statInfo.Df1;
    if isempty(statInfo.Df2)
        row.Df2 = NaN;
    else
        row.Df2 = statInfo.Df2;
    end
    row.StatThreshold  = statThrd;
    row.ZThreshold     = zThrd;
    row.MinClusterSize = minClusterSize;
    row.dLh            = dLh;
    row.FWHMx          = FWHM(1);
    row.FWHMy          = FWHM(2);
    row.FWHMz          = FWHM(3);
    row.SmoothnessSource = smoothnessSource;
    row.Em             = Em;
    row.EN             = EN;
    row.Beta           = Beta;

    rows(end+1,1) = row; %#ok<AGROW>
end

if isempty(rows)
    clusterTable = empty_detailed_table();
    return;
end

clusterTable = struct2table(rows);
clusterTable = normalize_table_text_columns_to_cell(clusterTable);

[~, ord] = sort(abs(clusterTable.PeakT), 'descend');
clusterTable = clusterTable(ord, :);

% Renumber after sorting to match the order used in the report.
clusterTable.ClusterID = (1:height(clusterTable))';

end

%% ========================================================================
function atlasResults = get_all_atlas_results(atlases, clusterIdx, peakIJK, clusterSizeVoxel)

atlasResults = repmat(make_empty_atlas_result(), numel(atlases), 1);
for a = 1:numel(atlases)
    atlasResults(a) = summarize_cluster_with_resliced_atlas(atlases(a), clusterIdx, peakIJK, clusterSizeVoxel);
end

end

%% ========================================================================
function result = make_empty_atlas_result()
result = struct();
result.Alias = '';
result.ShortName = '';
result.PeakID = 0;
result.PeakLabel = 'Unlabeled';
result.DominantID = 0;
result.DominantLabel = 'Unlabeled';
result.CoveragePct = 0;
result.DominantPct = 0;
result.Top3Text = 'Unlabeled';
end

%% ========================================================================
function result = atlas_result_by_shortname(atlasResults, shortName)
result = make_empty_atlas_result();
for i = 1:numel(atlasResults)
    if strcmpi(atlasResults(i).ShortName, shortName)
        result = atlasResults(i);
        return;
    end
end
end

%% ========================================================================
function result = summarize_cluster_with_resliced_atlas(atlas, clusterIdx, peakIJK, clusterSizeVoxel)

result = make_empty_atlas_result();
result.Alias = atlas.Alias;
result.ShortName = atlas.ShortName;

if isempty(atlas.Data)
    return;
end

try
    peakID = round(double(atlas.Data(peakIJK(1), peakIJK(2), peakIJK(3))));
catch
    peakID = 0;
end

result.PeakID = peakID;
result.PeakLabel = id_to_label(atlas.LabelMap, peakID);

ids = round(double(atlas.Data(clusterIdx)));
ids = ids(:);

nAll = clusterSizeVoxel;
if nAll <= 0 || isempty(ids)
    return;
end

valid = ids(ids > 0);
result.CoveragePct = 100 * numel(valid) / nAll;

if isempty(valid)
    return;
end

[u, ~, ic] = unique(valid);
counts = accumarray(ic, 1);
[countsSorted, ord] = sort(counts, 'descend');
uSorted = u(ord);

dominantID = uSorted(1);
dominantCount = countsSorted(1);

result.DominantID = dominantID;
result.DominantLabel = id_to_label(atlas.LabelMap, dominantID);

% DPABI Viewer reports label percentages relative to the full cluster size.
result.DominantPct = 100 * dominantCount / nAll;

nTop = min(3, numel(uSorted));
parts = cell(nTop, 1);
for i = 1:nTop
    id = uSorted(i);
    pct = 100 * countsSorted(i) / nAll;
    parts{i} = sprintf('[%d] %s %.1f%%', id, id_to_label(atlas.LabelMap, id), pct);
end
result.Top3Text = strjoin(parts, '; ');

end

%% ========================================================================
function [titleText, sentenceText] = build_report_text( ...
    roiName, directionLabel, clusterID, clusterSizeVoxel, clusterSizeMM3, ...
    peakXYZmm, peakT, peakZ, bn, aal, yeo, hoc, hos, needBoundaryCheck)

coordText = sprintf('[%d, %d, %d]', round(peakXYZmm(1)), round(peakXYZmm(2)), round(peakXYZmm(3)));

peakMain = choose_best_peak_anatomy(bn.PeakLabel, aal.PeakLabel, hoc.PeakLabel, hos.PeakLabel);
dominantMain = choose_best_peak_anatomy(bn.DominantLabel, aal.DominantLabel, hoc.DominantLabel, hos.DominantLabel);

titleText = sprintf('%s | %s | Cluster %d | peak %s', ...
    roiName, directionLabel, clusterID, peakMain);

sentenceText = sprintf([ ...
    'For seed %s, significant %s cluster %d peaks at %s. ' ...
    'Peak T = %.3f, peak Z = %.3f, and cluster size = %d voxels (%.1f mm^3). ' ...
    'The peak Brainnetome label is %s and the dominant cluster label is %s. ' ...
    'The AAL3 peak label is %s; the Yeo17 peak network is %s; ' ...
    'the HarvardOxford cortical peak label is %s; and the HarvardOxford subcortical peak label is %s.'], ...
    roiName, directionLabel, clusterID, coordText, ...
    peakT, peakZ, clusterSizeVoxel, clusterSizeMM3, ...
    bn.PeakLabel, dominantMain, aal.PeakLabel, yeo.PeakLabel, hoc.PeakLabel, hos.PeakLabel);

if ~strcmpi(bn.Top3Text, 'Unlabeled')
    sentenceText = [sentenceText, sprintf(' The three largest Brainnetome overlaps are: %s.', bn.Top3Text)];
end

if needBoundaryCheck
    sentenceText = [sentenceText, ' This cluster may cross an atlas boundary; inspect its full extent in DPABI Viewer.'];
end

end

%% ========================================================================
function best = choose_best_peak_anatomy(varargin)
best = 'Unlabeled';
for i = 1:nargin
    v = varargin{i};
    if ischar(v) && ~isempty(strtrim(v)) && ~strcmpi(strtrim(v), 'Unlabeled')
        best = strtrim(v);
        return;
    end
end
end

%% ========================================================================
function readableTable = make_readable_table(clusterTable)

if isempty(clusterTable)
    readableTable = empty_readable_table();
    return;
end

n = height(clusterTable);

ROI = cell(n,1);
Direction = cell(n,1);
ClusterID = zeros(n,1);
PeakMNI = cell(n,1);
PeakT = zeros(n,1);
PeakZ = zeros(n,1);
ClusterSizeVoxel = zeros(n,1);
ClusterSizeMM3 = zeros(n,1);
PeakAnatomy = cell(n,1);
DominantAnatomy = cell(n,1);
PeakAAL3 = cell(n,1);
PeakYeo17 = cell(n,1);
DominantYeo17 = cell(n,1);
Peak_HO_Cort = cell(n,1);
Peak_HO_Sub = cell(n,1);
SpatialNote = cell(n,1);
Summary = cell(n,1);

for r = 1:n
    coordText = sprintf('[%d, %d, %d]', ...
        round(clusterTable.PeakXmm(r)), ...
        round(clusterTable.PeakYmm(r)), ...
        round(clusterTable.PeakZmm(r)));

    targetAnatomy = choose_best_peak_anatomy( ...
        table_text(clusterTable, 'Peak_BN_Label', r), ...
        table_text(clusterTable, 'Peak_AAL3_Label', r), ...
        table_text(clusterTable, 'Peak_HO_Cort_Label', r), ...
        table_text(clusterTable, 'Peak_HO_Sub_Label', r));

    dominantAnatomy = choose_best_peak_anatomy( ...
        table_text(clusterTable, 'Dominant_BN_Label', r), ...
        table_text(clusterTable, 'Dominant_AAL3_Label', r), ...
        table_text(clusterTable, 'Dominant_HO_Cort_Label', r), ...
        table_text(clusterTable, 'Dominant_HO_Sub_Label', r));

    spatialNote = 'Peak and dominant cluster labels are consistent';
    if clusterTable.NeedBoundaryCheck(r) == 1
        spatialNote = 'Possible boundary crossing; inspect the full cluster extent in Viewer';
    end

    ROI{r} = table_text(clusterTable, 'ROI', r);
    Direction{r} = table_text(clusterTable, 'Direction', r);
    ClusterID(r) = clusterTable.ClusterID(r);
    PeakMNI{r} = coordText;
    PeakT(r) = clusterTable.PeakT(r);
    PeakZ(r) = clusterTable.PeakZ(r);
    ClusterSizeVoxel(r) = clusterTable.ClusterSizeVoxel(r);
    ClusterSizeMM3(r) = clusterTable.ClusterSizeMM3(r);

    PeakAnatomy{r} = targetAnatomy;
    DominantAnatomy{r} = dominantAnatomy;
    PeakAAL3{r} = table_text(clusterTable, 'Peak_AAL3_Label', r);
    PeakYeo17{r} = table_text(clusterTable, 'Peak_Yeo17_Label', r);
    DominantYeo17{r} = table_text(clusterTable, 'Dominant_Yeo17_Label', r);
    Peak_HO_Cort{r} = table_text(clusterTable, 'Peak_HO_Cort_Label', r);
    Peak_HO_Sub{r} = table_text(clusterTable, 'Peak_HO_Sub_Label', r);
    SpatialNote{r} = spatialNote;
    Summary{r} = table_text(clusterTable, 'ReportSentence', r);
end

readableTable = table(ROI, Direction, ClusterID, PeakMNI, PeakT, PeakZ, ...
    ClusterSizeVoxel, ClusterSizeMM3, PeakAnatomy, DominantAnatomy, ...
    PeakAAL3, PeakYeo17, DominantYeo17, Peak_HO_Cort, Peak_HO_Sub, SpatialNote, Summary);

readableTable = rename_table_vars(readableTable, { ...
    'ROI', 'Direction', 'ClusterID', 'PeakMNI', 'PeakT', 'PeakZ', ...
    'ClusterSizeVoxel', 'ClusterSizeMM3', 'PeakAnatomy', 'DominantAnatomy', ...
    'PeakAAL3', 'PeakYeo17', 'DominantYeo17', 'Peak_HO_Cort', 'Peak_HO_Sub', 'SpatialNote', 'Summary'});

end

%% ========================================================================
function T = make_byroi_cluster_sheet_table(clusterTable)

if isempty(clusterTable)
    T = empty_byroi_cluster_sheet_table();
    return;
end

n = height(clusterTable);
ROI = cell(n,1);
Direction = cell(n,1);
ClusterID = zeros(n,1);
PeakMNI = cell(n,1);
PeakT = zeros(n,1);
PeakZ = zeros(n,1);
ClusterSizeVoxel = zeros(n,1);
ClusterSizeMM3 = zeros(n,1);
PeakBrainnetome = cell(n,1);
PeakAAL3 = cell(n,1);
PeakYeo17 = cell(n,1);
PeakHOCortical = cell(n,1);
PeakHOSubcortical = cell(n,1);
PeakAnatomyBest = cell(n,1);
DominantBrainnetome = cell(n,1);
DominantAAL3 = cell(n,1);
DominantHOCortical = cell(n,1);
DominantHOSubcortical = cell(n,1);
BNTop3 = cell(n,1);
AAL3Top3 = cell(n,1);
Yeo17Top3 = cell(n,1);
HOCorticalTop3 = cell(n,1);
HOSubcorticalTop3 = cell(n,1);
NeedBoundaryCheck = zeros(n,1);

for r = 1:n
    ROI{r} = table_text(clusterTable, 'ROI', r);
    Direction{r} = table_text(clusterTable, 'Direction', r);
    ClusterID(r) = clusterTable.ClusterID(r);
    PeakMNI{r} = sprintf('[%d, %d, %d]', round(clusterTable.PeakXmm(r)), round(clusterTable.PeakYmm(r)), round(clusterTable.PeakZmm(r)));
    PeakT(r) = clusterTable.PeakT(r);
    PeakZ(r) = clusterTable.PeakZ(r);
    ClusterSizeVoxel(r) = clusterTable.ClusterSizeVoxel(r);
    ClusterSizeMM3(r) = clusterTable.ClusterSizeMM3(r);

    PeakBrainnetome{r} = table_text(clusterTable, 'Peak_BN_Label', r);
    PeakAAL3{r} = table_text(clusterTable, 'Peak_AAL3_Label', r);
    PeakYeo17{r} = table_text(clusterTable, 'Peak_Yeo17_Label', r);
    PeakHOCortical{r} = table_text(clusterTable, 'Peak_HO_Cort_Label', r);
    PeakHOSubcortical{r} = table_text(clusterTable, 'Peak_HO_Sub_Label', r);
    PeakAnatomyBest{r} = choose_best_peak_anatomy(PeakBrainnetome{r}, PeakAAL3{r}, PeakHOCortical{r}, PeakHOSubcortical{r});

    DominantBrainnetome{r} = table_text(clusterTable, 'Dominant_BN_Label', r);
    DominantAAL3{r} = table_text(clusterTable, 'Dominant_AAL3_Label', r);
    DominantHOCortical{r} = table_text(clusterTable, 'Dominant_HO_Cort_Label', r);
    DominantHOSubcortical{r} = table_text(clusterTable, 'Dominant_HO_Sub_Label', r);

    BNTop3{r} = table_text(clusterTable, 'BN_Top3', r);
    AAL3Top3{r} = table_text(clusterTable, 'AAL3_Top3', r);
    Yeo17Top3{r} = table_text(clusterTable, 'Yeo17_Top3', r);
    HOCorticalTop3{r} = table_text(clusterTable, 'HO_Cort_Top3', r);
    HOSubcorticalTop3{r} = table_text(clusterTable, 'HO_Sub_Top3', r);
    NeedBoundaryCheck(r) = clusterTable.NeedBoundaryCheck(r);
end

T = table(ROI, Direction, ClusterID, PeakMNI, PeakT, PeakZ, ...
    ClusterSizeVoxel, ClusterSizeMM3, PeakAnatomyBest, PeakBrainnetome, PeakAAL3, PeakYeo17, ...
    PeakHOCortical, PeakHOSubcortical, DominantBrainnetome, DominantAAL3, ...
    DominantHOCortical, DominantHOSubcortical, BNTop3, AAL3Top3, Yeo17Top3, ...
    HOCorticalTop3, HOSubcorticalTop3, NeedBoundaryCheck);

end

%% ========================================================================
function T = empty_byroi_cluster_sheet_table()
T = table(cell(0,1), cell(0,1), zeros(0,1), cell(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), cell(0,1), cell(0,1), cell(0,1), cell(0,1), ...
    cell(0,1), cell(0,1), cell(0,1), cell(0,1), cell(0,1), cell(0,1), ...
    cell(0,1), cell(0,1), cell(0,1), cell(0,1), cell(0,1), zeros(0,1));
T = rename_table_vars(T, {'ROI','Direction','ClusterID','PeakMNI','PeakT','PeakZ', ...
    'ClusterSizeVoxel','ClusterSizeMM3','PeakAnatomyBest','PeakBrainnetome','PeakAAL3','PeakYeo17', ...
    'PeakHOCortical','PeakHOSubcortical','DominantBrainnetome','DominantAAL3', ...
    'DominantHOCortical','DominantHOSubcortical','BNTop3','AAL3Top3','Yeo17Top3', ...
    'HOCorticalTop3','HOSubcorticalTop3','NeedBoundaryCheck'});
end

%% ========================================================================
function write_byroi_sheet(xlsxFile, roiIndex, def, T)
sheetName = make_excel_sheet_name(roiIndex, def.roiName);
writetable(T, xlsxFile, 'Sheet', sheetName);
end

%% ========================================================================
function sheetName = make_excel_sheet_name(roiIndex, roiName)
sheetName = sprintf('ROI%02d_%s', roiIndex, roiName);
sheetName = regexprep(sheetName, '[\[\]\:\*\?\/\\]', '_');
sheetName = regexprep(sheetName, '\s+', '_');
if numel(sheetName) > 31
    sheetName = sheetName(1:31);
end
end

%% ========================================================================
function write_readable_text_summary(txtFile, roiName, readableTable, runInfo)

fid = fopen(txtFile, 'w');
if fid < 0
    warning('Cannot write text file: %s', txtFile);
    return;
end

fprintf(fid, 'ROI: %s\n', roiName);
fprintf(fid, 'TestFlag = %s\n', runInfo.TestFlag);
fprintf(fid, 'Df1 = %g\n', runInfo.Df1);
if ~isempty(runInfo.Df2)
    fprintf(fid, 'Df2 = %g\n', runInfo.Df2);
end
fprintf(fid, 'VoxelP = %g\n', runInfo.VoxelP);
fprintf(fid, 'ClusterP = %g\n', runInfo.ClusterP);
fprintf(fid, 'ClusterPUsed = %g\n', runInfo.ClusterPUsed);
fprintf(fid, 'StatThreshold = %.6f\n', runInfo.StatThreshold);
fprintf(fid, 'ZThreshold = %.6f\n', runInfo.ZThreshold);
fprintf(fid, 'MinClusterSize = %d voxels\n', runInfo.MinClusterSize);
fprintf(fid, 'dLh = %.8f\n', runInfo.dLh);
fprintf(fid, 'FWHM = [%.4f %.4f %.4f]\n', runInfo.FWHM(1), runInfo.FWHM(2), runInfo.FWHM(3));
fprintf(fid, 'SmoothnessSource = %s\n', runInfo.SmoothnessSource);
fprintf(fid, '\n');

if isempty(readableTable)
    fprintf(fid, 'No significant clusters.\n');
    fclose(fid);
    return;
end

for i = 1:height(readableTable)
    fprintf(fid, '------------------------------------------------------------\n');
    fprintf(fid, '%s\n', table_text(readableTable, 'Summary', i));
    fprintf(fid, '\n');
end

fclose(fid);

end

%% ========================================================================
function write_overall_text_summary(txtFile, allReadable, roiSummary, cfg)

fid = fopen(txtFile, 'w');
if fid < 0
    warning('Cannot write overall text summary: %s', txtFile);
    return;
end

fprintf(fid, 'All ROIs summary\n');
fprintf(fid, 'Data root: %s\n', cfg.dataRoot);
fprintf(fid, 'VoxelP = %g\n', cfg.voxelP);
fprintf(fid, 'ClusterP = %g\n', cfg.clusterP);
fprintf(fid, 'TwoTailed = %d\n', double(cfg.isTwoTailed));
fprintf(fid, 'SmoothnessSource = %s\n', cfg.SMOOTHNESS_SOURCE);
fprintf(fid, '\n');

fprintf(fid, 'ROI-level summary\n');
for i = 1:height(roiSummary)
    fprintf(fid, '%s | Total=%d | Positive=%d | Negative=%d | %s\n', ...
        table_text(roiSummary, 'ROI', i), ...
        roiSummary.TotalClusters(i), ...
        roiSummary.PositiveClusters(i), ...
        roiSummary.NegativeClusters(i), ...
        table_text(roiSummary, 'Status', i));
end

fprintf(fid, '\n\nDetailed readable summaries\n');
if isempty(allReadable)
    fprintf(fid, 'No significant clusters found.\n');
else
    for i = 1:height(allReadable)
        fprintf(fid, '------------------------------------------------------------\n');
        fprintf(fid, '%s\n', table_text(allReadable, 'Summary', i));
        fprintf(fid, '\n');
    end
end

fclose(fid);

end

%% ========================================================================
function row = make_roi_summary_row(def, totalClusters, nPos, nNeg, statusText)

row = struct();
row.ROI = def.roiName;
row.Hemisphere = def.hemisphere;
row.SeedLabelID = def.labelID;
row.SeedRawName = def.rawName;
row.TotalClusters = totalClusters;
row.PositiveClusters = nPos;
row.NegativeClusters = nNeg;
row.Status = statusText;

end

%% ========================================================================
function T = append_table_like(T, addition)

if isempty(addition)
    return;
end

if isstruct(addition)
    addition = struct2table(addition);
elseif ~istable(addition)
    error('append_table_like accepts only a table or struct.');
end

addition = normalize_table_text_columns_to_cell(addition);

if isempty(T)
    T = addition;
    return;
end

T = normalize_table_text_columns_to_cell(T);

if ~isequal(T.Properties.VariableNames, addition.Properties.VariableNames)
    error('append_table_like: variable names differ between the two tables.');
end

T = [T; addition];

end

%% ========================================================================
function T = normalize_table_text_columns_to_cell(T)

vars = T.Properties.VariableNames;

for i = 1:numel(vars)
    col = T.(vars{i});

    if ischar(col)
        if isempty(col)
            T.(vars{i}) = cell(0,1);
        else
            T.(vars{i}) = cellstr(col);
        end

    elseif iscell(col)
        isTextLike = true(size(col));
        for k = 1:numel(col)
            v = col{k};
            isTextLike(k) = ischar(v) || isempty(v);
        end
        if all(isTextLike)
            tmp = cell(size(col));
            for k = 1:numel(col)
                v = col{k};
                if isempty(v)
                    tmp{k} = '';
                else
                    tmp{k} = char(v);
                end
            end
            T.(vars{i}) = tmp;
        end

    elseif iscategorical(col)
        T.(vars{i}) = cellstr(col);
    end
end

end

%% ========================================================================
function c = get_text_col_as_cellstr(col)

if iscell(col)
    c = col;
elseif ischar(col)
    c = cellstr(col);
elseif iscategorical(col)
    c = cellstr(col);
else
    n = numel(col);
    c = cell(n,1);
    for i = 1:n
        c{i} = safe_char(col(i));
    end
end

end

%% ========================================================================
function txt = table_text(T, varName, rowIdx)

txt = '';

if ~ismember(varName, T.Properties.VariableNames)
    return;
end

col = T.(varName);

if iscell(col)
    txt = safe_char(col{rowIdx});
elseif ischar(col)
    txt = safe_char(col(rowIdx,:));
elseif iscategorical(col)
    txt = safe_char(char(col(rowIdx)));
else
    v = col(rowIdx);
    if isnumeric(v) || islogical(v)
        txt = num2str(v);
    else
        txt = safe_char(v);
    end
end

end

%% ========================================================================
function s = safe_char(v)

if ischar(v)
    s = v;
elseif isnumeric(v) || islogical(v)
    if isscalar(v)
        s = num2str(v);
    else
        s = mat2str(v);
    end
elseif iscell(v) && isscalar(v)
    s = safe_char(v{1});
else
    try
        s = char(v);
    catch
        s = '';
    end
end

s = strtrim(s);

end

%% ========================================================================
function row = make_empty_detail_struct_row()

row = struct( ...
    'ROI', '', ...
    'Hemisphere', '', ...
    'SeedLabelID', 0, ...
    'SeedRawName', '', ...
    'Sign', '', ...
    'Direction', '', ...
    'ClusterID', 0, ...
    'ClusterSizeVoxel', 0, ...
    'ClusterSizeMM3', 0, ...
    'PeakT', 0, ...
    'PeakZ', 0, ...
    'MeanT', 0, ...
    'MeanZ', 0, ...
    'PeakI', 0, ...
    'PeakJ', 0, ...
    'PeakK', 0, ...
    'PeakXmm', 0, ...
    'PeakYmm', 0, ...
    'PeakZmm', 0, ...
    'Peak_BN_ID', 0, ...
    'Peak_BN_Label', '', ...
    'Dominant_BN_ID', 0, ...
    'Dominant_BN_Label', '', ...
    'BN_CoveragePct', 0, ...
    'BN_DominantPct', 0, ...
    'BN_Top3', '', ...
    'Peak_AAL3_ID', 0, ...
    'Peak_AAL3_Label', '', ...
    'Dominant_AAL3_ID', 0, ...
    'Dominant_AAL3_Label', '', ...
    'AAL3_CoveragePct', 0, ...
    'AAL3_DominantPct', 0, ...
    'AAL3_Top3', '', ...
    'Peak_Yeo17_ID', 0, ...
    'Peak_Yeo17_Label', '', ...
    'Dominant_Yeo17_ID', 0, ...
    'Dominant_Yeo17_Label', '', ...
    'Yeo17_CoveragePct', 0, ...
    'Yeo17_DominantPct', 0, ...
    'Yeo17_Top3', '', ...
    'Peak_HO_Cort_ID', 0, ...
    'Peak_HO_Cort_Label', '', ...
    'Dominant_HO_Cort_ID', 0, ...
    'Dominant_HO_Cort_Label', '', ...
    'HO_Cort_CoveragePct', 0, ...
    'HO_Cort_DominantPct', 0, ...
    'HO_Cort_Top3', '', ...
    'Peak_HO_Sub_ID', 0, ...
    'Peak_HO_Sub_Label', '', ...
    'Dominant_HO_Sub_ID', 0, ...
    'Dominant_HO_Sub_Label', '', ...
    'HO_Sub_CoveragePct', 0, ...
    'HO_Sub_DominantPct', 0, ...
    'HO_Sub_Top3', '', ...
    'NeedBoundaryCheck', 0, ...
    'ReportTitle', '', ...
    'ReportSentence', '', ...
    'VoxelP', 0, ...
    'ClusterP', 0, ...
    'ClusterPUsed', 0, ...
    'IsTwoTailed', 0, ...
    'TestFlag', '', ...
    'Df1', 0, ...
    'Df2', NaN, ...
    'StatThreshold', 0, ...
    'ZThreshold', 0, ...
    'MinClusterSize', 0, ...
    'dLh', 0, ...
    'FWHMx', 0, ...
    'FWHMy', 0, ...
    'FWHMz', 0, ...
    'SmoothnessSource', '', ...
    'Em', 0, ...
    'EN', 0, ...
    'Beta', 0);

end

%% ========================================================================
function T = empty_detailed_table()
rows = repmat(make_empty_detail_struct_row(), 0, 1);
T = struct2table(rows);
T = normalize_table_text_columns_to_cell(T);
end

%% ========================================================================
function T = empty_readable_table()

T = table( ...
    cell(0,1), cell(0,1), zeros(0,1), cell(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), cell(0,1), cell(0,1), cell(0,1), cell(0,1), ...
    cell(0,1), cell(0,1), cell(0,1), cell(0,1), cell(0,1));

T = rename_table_vars(T, { ...
    'ROI', 'Direction', 'ClusterID', 'PeakMNI', 'PeakT', 'PeakZ', ...
    'ClusterSizeVoxel', 'ClusterSizeMM3', 'PeakAnatomy', 'DominantAnatomy', ...
    'PeakAAL3', 'PeakYeo17', 'DominantYeo17', 'Peak_HO_Cort', 'Peak_HO_Sub', 'SpatialNote', 'Summary'});

end

%% ========================================================================
function T = empty_roi_summary_table()

T = table( ...
    cell(0,1), cell(0,1), zeros(0,1), cell(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), cell(0,1));

T = rename_table_vars(T, { ...
    'ROI', 'Hemisphere', 'SeedLabelID', 'SeedRawName', ...
    'TotalClusters', 'PositiveClusters', 'NegativeClusters', 'Status'});

end

%% ========================================================================
function T = rename_table_vars(T, names)

for i = 1:numel(names)
    T.Properties.VariableNames{i} = names{i};
end

end

%% ========================================================================
function must_exist(fname)
if exist(fname, 'file') ~= 2
    error('Required function was not found: %s', fname);
end
end

%% ========================================================================
function z = stat_to_z_signed(statData, statInfo)

switch upper(statInfo.TestFlag)
    case 'T'
        z = t_to_z_signed(statData, statInfo.Df1);
    case 'Z'
        z = statData;
    case 'R'
        df = statInfo.Df1;
        t = statData .* sqrt(df ./ max(1 - statData.^2, eps));
        z = t_to_z_signed(t, df);
    case 'F'
        % F statistics are unsigned; convert their upper-tail probability to Z.
        p = fcdf_safe(statData, statInfo.Df1, statInfo.Df2);
        eps0 = 1e-300;
        p(p < eps0) = eps0;
        p(p > 1 - eps0) = 1 - eps0;
        z = norminv_local(p);
    otherwise
        error('Unsupported TestFlag: %s', statInfo.TestFlag);
end

z(~isfinite(statData)) = NaN;

end

%% ========================================================================
function z = t_to_z_signed(t, df)

if ~(isscalar(df) && isfinite(df) && df > 0)
    error('Invalid Df1: %g', df);
end

p = tcdf_local(t, df);

eps0 = 1e-300;
p(p < eps0) = eps0;
p(p > 1 - eps0) = 1 - eps0;

z = norminv_local(p);
z(~isfinite(t)) = NaN;

end

%% ========================================================================
function p = tcdf_local(t, nu)

x = nu ./ (nu + t.^2);
ib = betainc(x, nu/2, 0.5);

p = zeros(size(t));

neg = t < 0;
zeroMask = (t == 0);
pos = t > 0;

p(neg) = 0.5 .* ib(neg);
p(zeroMask) = 0.5;
p(pos) = 1 - 0.5 .* ib(pos);

p(~isfinite(p)) = NaN;

end

%% ========================================================================
function x = norminv_local(p)
x = -sqrt(2) .* erfcinv(2 .* p);
end

%% ========================================================================
function p = normcdf_local(x)
p = 0.5 .* erfc(-x ./ sqrt(2));
end

%% ========================================================================
function x = tinv_safe(p, nu)
if exist('tinv', 'file') == 2
    x = tinv(p, nu);
    return;
end

% Scalar bisection fallback when the Statistics Toolbox is unavailable.
if ~(isscalar(p) && isscalar(nu))
    error('The tinv_safe fallback supports scalar inputs only.');
end
lo = -100;
hi = 100;
for iter = 1:200
    mid = (lo + hi)/2;
    if tcdf_local(mid, nu) < p
        lo = mid;
    else
        hi = mid;
    end
end
x = (lo + hi)/2;
end

%% ========================================================================
function x = finv_safe(p, df1, df2)
if exist('finv', 'file') == 2
    x = finv(p, df1, df2);
    return;
end
error('finv was not found for an F map. Install Statistics Toolbox or handle the conversion manually.');
end

%% ========================================================================
function p = fcdf_safe(x, df1, df2)
if exist('fcdf', 'file') == 2
    p = fcdf(x, df1, df2);
    return;
end
p = betainc((df1 .* x) ./ (df1 .* x + df2), df1/2, df2/2);
end

%% ========================================================================
function [minClusterSize, Em, EN, Beta] = calc_min_cluster_size(nVoxels, dLh, zThrd, clusterPUsed)

D = 3;

if zThrd^2 <= 1
    error('zThrd^2 <= 1; the voxel threshold is too low for a stable calculation.');
end

Em = nVoxels * (2*pi)^(-(D+1)/2) * dLh * (zThrd*zThrd - 1)^((D-1)/2) * exp(-zThrd*zThrd/2);
EN = nVoxels * (1 - normcdf_local(zThrd));

if Em <= 0 || EN <= 0 || ~isfinite(Em) || ~isfinite(EN)
    error('Invalid GRF parameters: Em=%g, EN=%g', Em, EN);
end

Beta = ((gamma(D/2 + 1) * Em) / EN)^(2 / D);

if ~isfinite(Beta) || Beta <= 0
    error('Invalid Beta value: %g', Beta);
end

pTemp = 1;
minClusterSize = 0;

while pTemp >= clusterPUsed
    minClusterSize = minClusterSize + 1;
    pTemp = 1 - exp(-Em * exp(-Beta * minClusterSize^(2 / D)));

    if minClusterSize > 1e7
        error('Minimum-cluster-size calculation did not converge; check the parameters.');
    end
end

end

%% ========================================================================
function pixelLists = connected_pixel_lists(binMask, connectivity)
% Return one linear-index vector per connected component.
pixelLists = cell(0, 1);
if ~any(binMask(:))
    return;
end
if exist('bwconncomp', 'file') == 2
    CC = bwconncomp(binMask, connectivity);
    pixelLists = CC.PixelIdxList(:);
else
    [labels, count] = bwlabeln(binMask, connectivity);
    pixelLists = cell(count, 1);
    for component = 1:count
        pixelLists{component} = find(labels == component);
    end
end
end

%% ========================================================================
function keepMask = keep_large_clusters(binMask, minClusterSize, connectivity)

keepMask = false(size(binMask));

if ~any(binMask(:))
    return;
end

if exist('bwconncomp', 'file') == 2
    CC = bwconncomp(binMask, connectivity);
    for i = 1:CC.NumObjects
        idx = CC.PixelIdxList{i};
        if numel(idx) >= minClusterSize
            keepMask(idx) = true;
        end
    end
else
    [L, num] = bwlabeln(binMask, connectivity);
    for i = 1:num
        idx = find(L == i);
        if numel(idx) >= minClusterSize
            keepMask(idx) = true;
        end
    end
end

end

%% ========================================================================
function xyzmm = voxel_to_mm(header, ijk)

A = get_affine_from_header(header);

v = [ijk(:); 1];
m = A * v;
xyzmm = m(1:3).';

end

%% ========================================================================
function A = get_affine_from_header(header)

A = [];

if isstruct(header)
    if isfield(header, 'mat') && isequal(size(header.mat), [4 4])
        A = header.mat;
    elseif isfield(header, 'Transform') && isstruct(header.Transform) ...
            && isfield(header.Transform, 'T') && isequal(size(header.Transform.T), [4 4])
        A = header.Transform.T;
    elseif isfield(header, 'srow_x') && isfield(header, 'srow_y') && isfield(header, 'srow_z')
        A = [header.srow_x; header.srow_y; header.srow_z; 0 0 0 1];
    end
end

if isempty(A)
    error('A 4-by-4 affine could not be obtained from the NIfTI header.');
end

end

%% ========================================================================
function niiFile = ensure_nifti_file(fileOrGz)

niiFile = fileOrGz;

if exist(fileOrGz, 'file')
    [p, n, e] = fileparts(fileOrGz);
    if strcmpi(e, '.gz')
        niiFile = fullfile(p, n);
        if ~exist(niiFile, 'file')
            fprintf('Unzipping atlas: %s\n', fileOrGz);
            gunzip(fileOrGz, p);
        end
    end
else
    error('File does not exist: %s', fileOrGz);
end

end

%% ========================================================================
function labelMap = build_yeo17_label_map()

labelMap = containers.Map('KeyType', 'char', 'ValueType', 'char');
labelMap('0')  = 'Unlabeled';
labelMap('1')  = 'VisualA';
labelMap('2')  = 'VisualB';
labelMap('3')  = 'SomMotA';
labelMap('4')  = 'SomMotB';
labelMap('5')  = 'DorsAttnA';
labelMap('6')  = 'DorsAttnB';
labelMap('7')  = 'SalVentAttnA';
labelMap('8')  = 'SalVentAttnB';
labelMap('9')  = 'LimbicA';
labelMap('10') = 'LimbicB';
labelMap('11') = 'ContA';
labelMap('12') = 'ContB';
labelMap('13') = 'ContC';
labelMap('14') = 'DefaultA';
labelMap('15') = 'DefaultB';
labelMap('16') = 'DefaultC';
labelMap('17') = 'TempPar';

end

%% ========================================================================
function labelMap = load_brainnetome_label_map(xlsxFile)

labelMap = containers.Map('KeyType', 'char', 'ValueType', 'char');
labelMap('0') = 'Unlabeled';

T = readtable(xlsxFile, 'PreserveVariableNames', true);
vars = T.Properties.VariableNames;
varsNorm = lower(regexprep(vars, '[^a-zA-Z0-9]', ''));

isIdCol = false(size(vars));
for i = 1:numel(vars)
    v = varsNorm{i};
    if ~isempty(strfind(v, 'id')) %#ok<STREMP>
        isIdCol(i) = true;
    end
end

prefTextPatterns = { ...
    'leftandrighthemisphere', ...
    'region', ...
    'labelname', ...
    'subregion', ...
    'name', ...
    'anatomy', ...
    'gyrus', ...
    'lobe'};

textColOrder = [];
for p = 1:numel(prefTextPatterns)
    hit = find(~cellfun('isempty', strfind(varsNorm, prefTextPatterns{p})));
    textColOrder = [textColOrder, hit(:)']; %#ok<AGROW>
end
textColOrder = unique(textColOrder, 'stable');

if isempty(textColOrder)
    for i = 1:numel(vars)
        col = T.(vars{i});
        if iscellstr_col(col) || isstring_col(col) || iscategorical(col)
            textColOrder(end+1) = i; %#ok<AGROW>
        end
    end
end

if ~any(isIdCol)
    error('No ID column was identified in the Brainnetome workbook.');
end

for c = find(isIdCol)
    colName = vars{c};
    colNorm = varsNorm{c};
    colData = T.(colName);

    if iscell(colData)
        tmp = nan(height(T), 1);
        for i = 1:height(T)
            if isnumeric(colData{i}) && isscalar(colData{i})
                tmp(i) = double(colData{i});
            end
        end
        colData = tmp;
    end

    if ~isnumeric(colData)
        continue;
    end

    hemiPrefix = '';
    if ~isempty(strfind(colNorm, 'left')) || (~isempty(colNorm) && colNorm(end) == 'l') %#ok<STREMP>
        hemiPrefix = 'Left_';
    elseif ~isempty(strfind(colNorm, 'right')) || (~isempty(colNorm) && colNorm(end) == 'r') %#ok<STREMP>
        hemiPrefix = 'Right_';
    end

    for i = 1:height(T)
        id = colData(i);
        if ~(isfinite(id) && id >= 0 && abs(id-round(id)) < 1e-8)
            continue;
        end
        id = round(id);

        label = '';
        for tc = textColOrder
            txt = get_table_text_cell(T.(vars{tc}), i);
            if ~isempty(txt)
                label = txt;
                break;
            end
        end

        if isempty(label)
            label = sprintf('ID_%d', id);
        end

        label = strtrim(label);
        label = regexprep(label, '\s+', ' ');
        label = strrep(label, '/', '_');

        if ~isempty(hemiPrefix) && id ~= 0
            lowlab = lower(label);
            if isempty(strfind(lowlab, 'left_')) && isempty(strfind(lowlab, 'right_')) %#ok<STREMP>
                label = [hemiPrefix label]; %#ok<AGROW>
            end
        end

        key = num2str(id);
        if ~isKey(labelMap, key)
            labelMap(key) = label;
        end
    end
end

if labelMap.Count <= 1
    error('Failed to parse the Brainnetome label table.');
end

end

%% ========================================================================
function labelMap = load_mat_label_map(matFile)

labelMap = containers.Map('KeyType', 'char', 'ValueType', 'char');
labelMap('0') = 'Unlabeled';

S = load(matFile);

% DPABI Labels.mat files normally store labels in Reference(:,1) and IDs in Reference(:,2).
if isfield(S, 'Reference')
    labelMap = parse_reference_cell_to_label_map(S.Reference, labelMap);
    if labelMap.Count > 1
        return;
    end
end

% Fallback parser for alternative MAT label structures.
labelMap2 = try_parse_label_map_from_struct(S);
if labelMap2.Count > 0
    ks = labelMap2.keys;
    for i = 1:numel(ks)
        labelMap(ks{i}) = labelMap2(ks{i});
    end
end

if labelMap.Count <= 1
    error('Failed to parse label MAT file %s. Expected a Reference variable or an ID/name table.', matFile);
end

end

%% ========================================================================
function labelMap = parse_reference_cell_to_label_map(Reference, labelMap)

if ~iscell(Reference)
    return;
end

[m, n] = size(Reference);
if m == 0 || n < 2
    return;
end

bestNameCol = [];
bestIdCol = [];
bestScore = -inf;

for nameCol = 1:n
    for idCol = 1:n
        if nameCol == idCol
            continue;
        end
        score = 0;
        for r = 1:m
            nm = Reference{r, nameCol};
            id = Reference{r, idCol};
            if ischar(nm) && ~isempty(strtrim(nm)) && isnumeric(id) && isscalar(id) && isfinite(id)
                score = score + 1;
            end
        end
        if score > bestScore
            bestScore = score;
            bestNameCol = nameCol;
            bestIdCol = idCol;
        end
    end
end

if isempty(bestNameCol) || bestScore == 0
    return;
end

for r = 1:m
    nm = Reference{r, bestNameCol};
    id = Reference{r, bestIdCol};
    if ischar(nm) && isnumeric(id) && isscalar(id) && isfinite(id)
        id = round(double(id));
        nm = strtrim(nm);
        if isempty(nm)
            nm = sprintf('ID_%d', id);
        end
        labelMap(num2str(id)) = nm;
    end
end

end

%% ========================================================================
function labelMap = try_parse_label_map_from_struct(S)

labelMap = containers.Map('KeyType', 'char', 'ValueType', 'char');

if ~isstruct(S)
    return;
end

fns = fieldnames(S);
numCandidates = {};
numNames = {};
txtCandidates = {};
txtNames = {};

for i = 1:numel(fns)
    v = S.(fns{i});

    if iscell(v)
        tmp = try_parse_label_map_from_cell(v);
        if tmp.Count > 0
            labelMap = tmp;
            return;
        end
    end

    if isnumeric(v) && isvector(v) && numel(v) > 1
        numCandidates{end+1} = double(v(:)); %#ok<AGROW>
        numNames{end+1} = lower(fns{i}); %#ok<AGROW>
    elseif iscell(v) && isvector(v) && numel(v) > 1
        alltxt = true;
        for k = 1:numel(v)
            alltxt = alltxt && (ischar(v{k}) || isempty(v{k}));
        end
        if alltxt
            txtCandidates{end+1} = v(:); %#ok<AGROW>
            txtNames{end+1} = lower(fns{i}); %#ok<AGROW>
        end
    elseif iscategorical(v) && isvector(v) && numel(v) > 1
        txtCandidates{end+1} = cellstr(v(:)); %#ok<AGROW>
        txtNames{end+1} = lower(fns{i}); %#ok<AGROW>
    end
end

bestScore = -inf;
bestNum = [];
bestTxt = [];

for i = 1:numel(numCandidates)
    ids = numCandidates{i};
    for j = 1:numel(txtCandidates)
        names = txtCandidates{j};
        if numel(ids) ~= numel(names)
            continue;
        end
        score = 0;
        if ~isempty(strfind(numNames{i}, 'id')) || ~isempty(strfind(numNames{i}, 'index'))
            score = score + 2;
        end
        if ~isempty(strfind(txtNames{j}, 'label')) || ~isempty(strfind(txtNames{j}, 'name')) || ~isempty(strfind(txtNames{j}, 'nom'))
            score = score + 2;
        end
        if score > bestScore
            bestScore = score;
            bestNum = ids;
            bestTxt = names;
        end
    end
end

if ~isempty(bestNum)
    for k = 1:numel(bestNum)
        id = bestNum(k);
        if ~(isfinite(id) && id >= 0 && abs(id-round(id)) < 1e-8)
            continue;
        end
        id = round(id);
        name = strtrim(bestTxt{k});
        if isempty(name)
            name = sprintf('ID_%d', id);
        end
        labelMap(num2str(id)) = name;
    end
end

end

%% ========================================================================
function labelMap = try_parse_label_map_from_cell(C)

labelMap = containers.Map('KeyType', 'char', 'ValueType', 'char');

if ~iscell(C)
    return;
end

[m, n] = size(C);
if m == 0 || n == 0
    return;
end

for idCol = 1:n
    for txtCol = 1:n
        if txtCol == idCol
            continue;
        end

        ids = nan(m,1);
        okId = false(m,1);
        txts = cell(m,1);
        okTxt = false(m,1);

        for i = 1:m
            a = C{i, idCol};
            b = C{i, txtCol};

            if isnumeric(a) && isscalar(a) && isfinite(a)
                ids(i) = double(a);
                okId(i) = true;
            end

            if ischar(b)
                txts{i} = b;
                okTxt(i) = ~isempty(strtrim(b));
            end
        end

        good = okId & okTxt;
        if nnz(good) >= 5
            for i = find(good(:))'
                id = round(ids(i));
                if id >= 0
                    labelMap(num2str(id)) = strtrim(txts{i});
                end
            end
            return;
        end
    end
end

end

%% ========================================================================
function label = id_to_label(labelMap, id)

key = num2str(round(double(id)));
if isKey(labelMap, key)
    label = labelMap(key);
else
    if id == 0
        label = 'Unlabeled';
    else
        label = sprintf('ID_%d', round(double(id)));
    end
end

end

%% ========================================================================
function tf = iscellstr_col(col)
tf = iscell(col) && all(cellfun(@(x) ischar(x) || isempty(x), col(:)));
end

%% ========================================================================
function tf = isstring_col(col)
try
    tf = isstring(col) && isvector(col);
catch
    tf = false;
end
end

%% ========================================================================
function txt = get_table_text_cell(col, rowIdx)

txt = '';

try
    if iscell(col)
        v = col{rowIdx};
        if ischar(v)
            txt = strtrim(v);
        end
    elseif iscategorical(col)
        txt = strtrim(char(col(rowIdx)));
    elseif ischar(col)
        txt = strtrim(col(rowIdx,:));
    end
catch
    txt = '';
end
end

%% ========================================================================
function cfg = cfg_default(cfg, fieldName, defaultValue)
% Set a configuration default without overriding an explicitly supplied value.
if ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
    cfg.(fieldName) = defaultValue;
end
end
