function extract_falff_individual_rois(manifestCsv, roiDefinitionCsv, outputCsv, dpabiRoot)
%EXTRACT_FALFF_INDIVIDUAL_ROIS Extract mean fALFF from subject-specific ROIs.
%
% Required manifest columns:
%   SubjectID   unique participant identifier
%   FalffFile   participant fALFF NIfTI path; relative paths are resolved
%               from the manifest directory
%
% Required ROI-definition columns:
%   ROIName     public output label
%   MaskFile    mask path or pattern; the token {SubjectID} is replaced
%               with the current participant identifier, and relative paths
%               are resolved from the ROI-definition directory
%
% Example:
%   extract_falff_individual_rois(...
%       'config/falff_manifest.csv', ...
%       'config/roi_definitions.csv', ...
%       'derived/falff_roi_values.csv', ...
%       getenv('DPABI_ROOT'));
%
arguments
    manifestCsv (1,:) char
    roiDefinitionCsv (1,:) char
    outputCsv (1,:) char
    dpabiRoot (1,:) char
end

if isempty(dpabiRoot) || ~isfolder(dpabiRoot)
    error('Set DPABI_ROOT or pass a valid DPABI installation directory.');
end
addpath(genpath(dpabiRoot));
require_function('y_ReadRPI');
require_function('y_Reslice');

manifest = readtable(manifestCsv, 'TextType', 'string', 'PreserveVariableNames', true);
roiDefs = readtable(roiDefinitionCsv, 'TextType', 'string', 'PreserveVariableNames', true);
require_columns(manifest, ["SubjectID", "FalffFile"], manifestCsv);
require_columns(roiDefs, ["ROIName", "MaskFile"], roiDefinitionCsv);

if numel(unique(manifest.SubjectID)) ~= height(manifest)
    error('SubjectID must be unique in %s.', manifestCsv);
end
if numel(unique(roiDefs.ROIName)) ~= height(roiDefs)
    error('ROIName must be unique in %s.', roiDefinitionCsv);
end

outputDir = fileparts(outputCsv);
if isempty(outputDir)
    outputDir = pwd;
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end
resliceDir = fullfile(outputDir, 'resliced_roi_masks');
if ~isfolder(resliceDir)
    mkdir(resliceDir);
end

nSubjects = height(manifest);
nRois = height(roiDefs);
values = nan(nSubjects, nRois);
voxelCounts = zeros(nSubjects, nRois);

for subjectIndex = 1:nSubjects
    subjectID = manifest.SubjectID(subjectIndex);
    falffFile = resolve_table_path(manifestCsv, manifest.FalffFile(subjectIndex));
    if ~isfile(falffFile)
        warning('Missing fALFF file for %s: %s', subjectID, falffFile);
        continue;
    end

    [falffData, falffVoxelSize, ~] = y_ReadRPI(char(falffFile));
    falffData = squeeze(double(falffData));
    if ndims(falffData) ~= 3
        error('Expected a 3-D fALFF map for %s.', subjectID);
    end

    for roiIndex = 1:nRois
        maskPattern = roiDefs.MaskFile(roiIndex);
        maskFile = strrep(maskPattern, "{SubjectID}", subjectID);
        maskFile = resolve_table_path(roiDefinitionCsv, maskFile);
        if ~isfile(maskFile)
            warning('Missing mask for %s / %s: %s', ...
                subjectID, roiDefs.ROIName(roiIndex), maskFile);
            continue;
        end

        safeSubject = matlab.lang.makeValidName(char(subjectID));
        safeRoi = matlab.lang.makeValidName(char(roiDefs.ROIName(roiIndex)));
        reslicedMask = fullfile(resliceDir, sprintf('%s__%s.nii', safeSubject, safeRoi));
        y_Reslice(char(maskFile), reslicedMask, falffVoxelSize, 0, char(falffFile));
        [maskData, ~, ~] = y_ReadRPI(reslicedMask);
        maskData = squeeze(maskData);
        if ~isequal(size(maskData), size(falffData))
            error('Resliced mask dimensions do not match fALFF for %s / %s.', ...
                subjectID, roiDefs.ROIName(roiIndex));
        end

        valid = maskData > 0.5 & isfinite(falffData);
        voxelCounts(subjectIndex, roiIndex) = nnz(valid);
        if any(valid(:))
            values(subjectIndex, roiIndex) = mean(falffData(valid), 'omitnan');
        end
    end
end

variableNames = matlab.lang.makeUniqueStrings(...
    matlab.lang.makeValidName(cellstr(roiDefs.ROIName)));
valueTable = array2table(values, 'VariableNames', variableNames);
valueTable = addvars(valueTable, manifest.SubjectID, 'Before', 1, 'NewVariableNames', 'SubjectID');
writetable(valueTable, outputCsv);

qcTable = array2table(voxelCounts, 'VariableNames', variableNames);
qcTable = addvars(qcTable, manifest.SubjectID, 'Before', 1, 'NewVariableNames', 'SubjectID');
[qcDir, qcStem, ~] = fileparts(outputCsv);
writetable(qcTable, fullfile(qcDir, [qcStem '_voxel_counts.csv']));
end


function pathValue = resolve_table_path(tableFile, suppliedPath)
suppliedPath = char(suppliedPath);
isWindowsAbsolute = numel(suppliedPath) >= 3 && isletter(suppliedPath(1)) ...
    && suppliedPath(2) == ':' && any(suppliedPath(3) == ['\' '/']);
isUnixAbsolute = ~isempty(suppliedPath) && suppliedPath(1) == '/';
if isWindowsAbsolute || isUnixAbsolute
    pathValue = suppliedPath;
else
    tableDir = fileparts(char(tableFile));
    if isempty(tableDir)
        tableDir = pwd;
    end
    pathValue = fullfile(tableDir, suppliedPath);
end
end


function require_columns(tableValue, required, source)
missing = setdiff(required, string(tableValue.Properties.VariableNames));
if ~isempty(missing)
    error('%s is missing columns: %s', source, strjoin(missing, ', '));
end
end


function require_function(functionName)
if exist(functionName, 'file') ~= 2
    error('Required DPABI function is unavailable: %s', functionName);
end
end
