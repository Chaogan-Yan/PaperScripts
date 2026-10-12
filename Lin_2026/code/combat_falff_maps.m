function combat_falff_maps(manifestCsv, brainMaskFile, outputDir, dpabiRoot)
%COMBAT_FALFF_MAPS Apply parametric ComBat to participant fALFF maps.
%
% This is a portable reference implementation of the documented model. The
% reported study workflow used DPABI Harmonization through its graphical
% interface; the filled GUI provenance record remains authoritative.
%
% Required manifest columns (relative FalffFile paths are resolved from the
% manifest directory):
%   SubjectID, FalffFile, Site, Age, Sex, Motion, Diagnosis
%
% Site is the batch variable. Age, sex, mean FD (Motion), and diagnosis are
% preserved in the ComBat design matrix, matching the manuscript. The DPABI
% distribution must provide y_ReadAll, y_ReadRPI, y_Write, and combat.
%
% This function intentionally accepts paths as arguments and contains no
% workstation-specific directory names. For large datasets, memory use is
% comparable to DPABI Harmonization because the empirical-Bayes priors are
% estimated jointly over all retained voxels.

arguments
    manifestCsv (1,:) char
    brainMaskFile (1,:) char
    outputDir (1,:) char
    dpabiRoot (1,:) char
end

if isempty(dpabiRoot) || ~isfolder(dpabiRoot)
    error('Set DPABI_ROOT or pass a valid DPABI installation directory.');
end
addpath(genpath(dpabiRoot));
require_function('y_ReadAll');
require_function('y_ReadRPI');
require_function('y_Write');
require_function('combat');

manifest = readtable(manifestCsv, 'TextType', 'string', 'PreserveVariableNames', true);
required = ["SubjectID", "FalffFile", "Site", "Age", "Sex", "Motion", "Diagnosis"];
require_columns(manifest, required, manifestCsv);
if numel(unique(manifest.SubjectID)) ~= height(manifest)
    error('SubjectID must be unique.');
end
if any(ismissing(manifest(:, required)), 'all')
    error('ComBat manifest contains missing required values.');
end
resolvedFiles = strings(height(manifest), 1);
for subjectIndex = 1:height(manifest)
    resolvedFiles(subjectIndex) = resolve_table_path(...
        manifestCsv, manifest.FalffFile(subjectIndex));
end
if any(~isfile(resolvedFiles))
    bad = resolvedFiles(find(~isfile(resolvedFiles), 1));
    error('Missing fALFF file: %s', bad);
end
if ~isfile(brainMaskFile)
    error('Missing brain mask: %s', brainMaskFile);
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

imageFiles = cellstr(resolvedFiles);
[allMaps, ~, ~, header] = y_ReadAll(imageFiles);
[maskData, ~, ~] = y_ReadRPI(brainMaskFile);
maskData = logical(squeeze(maskData));

mapSize = size(allMaps);
if numel(mapSize) ~= 4 || mapSize(4) ~= height(manifest)
    error('Expected one 3-D fALFF map per manifest row.');
end
if ~isequal(mapSize(1:3), size(maskData))
    error('Brain-mask dimensions do not match the fALFF maps.');
end

flatMaps = reshape(double(allMaps), [], mapSize(4));
maskIndex = find(maskData(:));
features = flatMaps(maskIndex, :);
eligible = all(isfinite(features), 2) & std(features, 0, 2) > 0;
if ~any(eligible)
    error('No finite, nonconstant in-mask voxels are available for ComBat.');
end

batch = findgroups(categorical(manifest.Site));
sexDesign = dummyvar(categorical(manifest.Sex));
diagnosisDesign = dummyvar(categorical(manifest.Diagnosis));
mod = [double(manifest.Age), sexDesign(:, 2:end), double(manifest.Motion), ...
    diagnosisDesign(:, 2:end)];

% The final argument 1 selects parametric empirical-Bayes adjustment.
harmonized = combat(features(eligible, :), batch, mod, 1);
features(eligible, :) = harmonized;
flatMaps(maskIndex, :) = features;
allMapsHarmonized = reshape(flatMaps, mapSize);

for subjectIndex = 1:height(manifest)
    [~, stem, ~] = fileparts(resolvedFiles(subjectIndex));
    outputFile = fullfile(outputDir, sprintf('%s__%s_ComBat.nii', ...
        manifest.SubjectID(subjectIndex), stem));
    outputHeader = header;
    outputHeader.pinfo = [1; 0; 0];
    outputHeader.dt = [16, 0];
    y_Write(single(allMapsHarmonized(:, :, :, subjectIndex)), outputHeader, outputFile);
end
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
