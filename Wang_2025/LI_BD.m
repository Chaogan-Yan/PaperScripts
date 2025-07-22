%This script calculates the region-level  lateralization index (LI) of SC and FC 
%Each subject will receive 454 LI values corresponding to each brain region

%%
%LI SC
clear; clc;

data_dir = uigetdir('', 'Select SC folder');
if data_dir == 0, error('No folder selected'); end

threshold = 0.2;
output_dir = fullfile(data_dir, 'LI_Results');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

LH = [1:200, 401:427]; RH = [201:400, 428:454];
files = dir(fullfile(data_dir, '*.mat'));
if isempty(files), error('No .mat files'); end

all_LI = zeros(454, length(files));

for i = 1:length(files)
    file = fullfile(data_dir, files(i).name);
    D = load(file, 'Data');
    SC = D.Data;

    if ~isequal(size(SC), [454, 454])
        warning('Skipped: %s', files(i).name); continue;
    end

    SC(eye(454) == 1) = 0;
    LI = zeros(454, 1);

    for r = 1:454
        intra = SC(r, ismember(1:454, r <= 227 && LH || RH));
        inter = SC(r, ismember(1:454, r <= 227 && RH || LH));
        intra(intra <= threshold) = NaN;
        inter(inter <= threshold) = NaN;
        li = nanmean(inter) - nanmean(intra);
        LI(r) = ifelse(isnan(li), 0, li); % 替换为标准if
    end

    [~, name] = fileparts(files(i).name);
    save(fullfile(output_dir, ['LI_', name, '.mat']), 'LI');
    all_LI(:, i) = LI;
    fprintf('Done: %s (%d/%d)\n', name, i, length(files));
end

save(fullfile(output_dir, 'AllSubjects_LI.mat'), 'all_LI', 'files');

figure; plot(all_LI(:, 1)); grid on;
xlabel('Region'); ylabel('LI'); title('Subject 1');

%%
%LIFC
clear; clc;

data_dir = uigetdir('', 'Select FC folder');
if data_dir == 0, error('No folder selected'); end

threshold = 0.2;
output_dir = fullfile(data_dir, 'LI_Results');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

LH = [1:200, 401:427]; RH = [201:400, 428:454];
files = dir(fullfile(data_dir, '*.mat'));
if isempty(files), error('No .mat files'); end

all_LI = zeros(454, length(files));

for i = 1:length(files)
    file = fullfile(data_dir, files(i).name);
    D = load(file, 'Data');
    FC = D.Data;

    if ~isequal(size(FC), [454, 454])
        warning('Skipped: %s', files(i).name); continue;
    end

    FC(eye(454) == 1) = 0;
    LI = zeros(454, 1);

    for r = 1:454
        intra = FC(r, ismember(1:454, r <= 227 && LH || RH));
        inter = FC(r, ismember(1:454, r <= 227 && RH || LH));
        intra(intra <= threshold) = NaN;
        inter(inter <= threshold) = NaN;
        li = nanmean(inter) - nanmean(intra);
        LI(r) = ifelse(isnan(li), 0, li); % 替换为标准if
    end

    [~, name] = fileparts(files(i).name);
    save(fullfile(output_dir, ['LI_', name, '.mat']), 'LI');
    all_LI(:, i) = LI;
    fprintf('Done: %s (%d/%d)\n', name, i, length(files));
end

save(fullfile(output_dir, 'AllSubjects_LI.mat'), 'all_LI', 'files');

figure; plot(all_LI(:, 1)); grid on;
xlabel('Region'); ylabel('LI'); title('Subject 1');


