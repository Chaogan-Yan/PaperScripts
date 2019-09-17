clear
dataDir = '/mnt/Data/RfMRILab/Lile/DFC/data_movie/S3_FunImgRCWFS';
roiDir = '/mnt/Data/RfMRILab/Lile/DFCnew/movieGmd/DFCstd/ROI_MV-RS_W';
outputDir = '/mnt/Data/RfMRILab/Lile/DFCnew/movieGmd/DFCstd/std_RS2';
mkdir(outputDir)
WinLength = 80; WinStep = 5; WinType = 'rectwin';
% 
sublist = dir([dataDir, '/sub*']);
roiPath = spm_select('FPList', roiDir, '.nii$');
[mask, hdr] = y_Read('/mnt/Data/RfMRILab/Lile/DFCnew/movie/GroupMask_GM.nii');
idxGM = find(mask>0);
hdr.pinfo = [1;0;0];
hdr.dt    =[16,0];
tic
for rn = 1:size(roiPath, 1)
    [roiMask, ~] = y_Read(deblank(roiPath(rn,:)));
    roiMask = roiMask>0;
    roiGM = roiMask(idxGM);
    [~, roiName, ~] = fileparts(deblank(roiPath(rn,:)));
    si = strfind(roiName, '_');
    roiName = roiName(1:si(2)-1);
    saveDir = [outputDir, '/', roiName];
    mkdir(saveDir)
    parfor sn = 1:length(sublist)
        data = y_ReadAll([dataDir, '/' sublist(sn).name]);
        tp(sn) = size(data,4);
        data = reshape(data, [], tp(sn))';
        % confined analysis to GM
        data = data(:, idxGM);
        roiTS = mean(data(:, roiGM),2);
        % calculate seed-based dfc
        dfc = [];
        for iWin = 1 : WinStep : (tp(sn)-WinLength+1)
            data_win = data(iWin : iWin+WinLength-1, :);
            roi_win = roiTS(iWin : iWin+WinLength-1, :);
            data_win = (data_win-repmat(mean(data_win),size(data_win,1),1))./repmat(std(data_win),size(data_win,1),1);
            roi_win = (roi_win-repmat(mean(roi_win),size(roi_win,1),1))./repmat(std(roi_win),size(roi_win,1),1);
            fc = data_win' * roi_win / (WinLength-1);
            dfc = cat(2, dfc, fc);
        end
        dfcStd = std(dfc, 1, 2);
        dfcStdImg = zeros(size(mask));
        dfcStdImg(idxGM) = dfcStd;
        Header = hdr;
        Header.fname = sprintf('%s/V_%s.nii', saveDir, sublist(sn).name);
        y_Write(dfcStdImg, Header, Header.fname)
        fprintf('dfc std computation is done for %s \n\n', sublist(sn).name)
    end
    toc
end
disp('done------done------done------done------done------done')