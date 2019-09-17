function f_DFC_W(data_path, tpl_path, output_path, mask_path, WinLength, WinStep, WinType)
% Compute Kendall's W of dynamic functional connectivity (DFC) with windows
% as rater, for one subject. The W quantifies functinal stability of the 
% brain organization (from a dynamic perspective, ~1 min).

% read 4D data and tranform to 2D with rows for each time point
[data_all, ~] = y_ReadAll(data_path);
tp = size(data_all, 4);
data_all = reshape(data_all, [], tp)';
% confined analysis to available grey matter voxels
[mask, Header] = y_Read(mask_path);
idx_GM = find(mask(:)>0 & mean(data_all)'~=0);
if length(idx_GM) < length(find(mask>0))*0.9
    fprintf('-------only %d voxels is covered within the brain.-------\n', length(idx_GM))
    return
end
dataP1 = data_all(:, idx_GM);


% conduct a voxel-voxel or voxel-atlas analysis, depending on if tpl_path exists
if isempty(tpl_path)
    dataP2 = dataP1;
elseif isnumeric(tpl_path)
    dataP2 = tpl_path;
else
    idx_WB = mean(data_all)'~=0;
    data_all = data_all(:, idx_WB);
    [tpl, ~] = y_Read(tpl_path);
    tpl = tpl(idx_WB);
    for j = 1:max(tpl(:))
        tpl_ts(:, j) = mean(data_all(:, tpl==j), 2);
    end
    miss_area = find(isnan(mean(tpl_ts)));
    tpl_ts(:, miss_area) = [];
    if ~isempty(miss_area)
        fprintf('data in the area %s of template are missing \n', num2str(miss_area))
    end
    dataP2 = tpl_ts;
end
clear data_all

%% calculate dynamic FC and W across voxels for each segment
if ~isempty(WinType)
    eval(['WinType=', WinType, '(WinLength);'])
else
    WinType = rectwin(WinLength);
end

if ~isempty(tpl_path), CutNumber = 24; else CutNumber = 3000; end
SegmentLength = ceil(size(dataP1,2) / CutNumber);
CutNumber = ceil(size(dataP1,2) / SegmentLength);
parfor iCut = 1:CutNumber
    if iCut~=CutNumber
        Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
    else
        Segment = (iCut-1)*SegmentLength+1 : size(dataP1,2);
    end
    dataP1_Segment = dataP1(:, Segment);
    WinMultiplierP1 = repmat(WinType, 1, length(Segment));
    WinMultiplierP2 = repmat(WinType, 1, size(dataP2, 2));
    dfc_Segment = zeros(length(Segment), size(dataP2,2), fix((tp-WinLength+1)/WinStep));
    for iWin = 1 : WinStep : (tp-WinLength+1)
        dataP1_win = dataP1_Segment(iWin : iWin+WinLength-1, :);
        dataP2_win = dataP2(iWin : iWin+WinLength-1, :);
        dataP1_win = dataP1_win .* WinMultiplierP1;
        dataP2_win = dataP2_win .* WinMultiplierP2;
        
        dataP1_win = (dataP1_win-repmat(mean(dataP1_win),size(dataP1_win,1),1))./repmat(std(dataP1_win),size(dataP1_win,1),1);
        dataP2_win = (dataP2_win-repmat(mean(dataP2_win),size(dataP2_win,1),1))./repmat(std(dataP2_win),size(dataP2_win,1),1);
        fc_Segment = dataP1_win' * dataP2_win / (WinLength-1);
        dfc_Segment(:, :, ceil(iWin/WinStep)) = fc_Segment; %%% voxel x template x window
    end
    dfc_rankset = tiedrank(permute(dfc_Segment, [2 3 1]));  %%% template x window x voxel
    W_Segment=zeros(SegmentLength, 1);
    if rand(1)<0.05; fprintf('there are %d windows and %d voxels \n', size(dfc_rankset, 2), size(dfc_rankset, 3)); end
    for iVoxel = 1:length(Segment)
        W_Segment(iVoxel) = f_kendall(dfc_rankset(:, :, iVoxel));
    end
    W_GM(:, iCut) = W_Segment;
end
W_GM = W_GM(:);
W_img = zeros(size(mask));
W_img(idx_GM) = W_GM(1:length(idx_GM));
zW_img = zeros(size(mask));
zW_img(idx_GM) = zscore(W_GM(1:length(idx_GM)));
%% write result to nii file
[~, subname, ~] = fileparts(data_path);
Header.pinfo = [1;0;0];
Header.dt    =[16,0];
Header.fname = [output_path, '/W_', subname,  '.nii'];
y_Write(W_img, Header, [output_path, '/W_', subname,  '.nii'])
Header.fname = [output_path, '/zW_', subname,  '.nii'];
y_Write(zW_img, Header, [output_path, '/zW_', subname,  '.nii'])
return
