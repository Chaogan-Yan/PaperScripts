VarPath = '/mnt/Data/RfMRILab/Lile/DFCnew/movie/movie_subinfo.xlsx';
DataDir='/mnt/Data/RfMRILab/Lile/DFC/data_movie/FunImgRCWFS';
OutDir='//mnt/Data/RfMRILab/Lile/940/movie_gmd20/KCC_MV2';
MaskFile = '/mnt/Data/RfMRILab/Lile/DFCnew/movieGmd/GroupMask_gmd20.nii';
tpl_path = '';

VarAll = xlsread(VarPath, 1, 'b:k');
[~, SubID] = xlsread(VarPath, 1, 'a2:a1000');

Sex = VarAll(:,1);
Age = VarAll(:,2);
FD = VarAll(:,3);
VoxelNum = VarAll(:,4);
QCscore = VarAll(:,6);
ExcessiveHM = VarAll(:,5);
WantedSubMatrix = ones(length(SubID), 1);

% %%%1. Deal With Age
% AgeRange = [16 32];
% AgeTemp = (Age>=AgeRange(1)) .* (Age<=AgeRange(2));
% WantedSubMatrix = WantedSubMatrix.*AgeTemp;

%%%2. Deal with head motion FD
FDTemp = (FD<(mean(FD)+2*std(FD)))*1;
WantedSubMatrix = WantedSubMatrix.*FDTemp;

%%%3. Deal with bad coverage
VoxelNumTemp = (VoxelNum>60000)*1;
WantedSubMatrix=WantedSubMatrix .* VoxelNumTemp;

%%%4. Deal with low QC score
QCTemp = (QCscore>=4)*1;
WantedSubMatrix=WantedSubMatrix .* QCTemp;

%%%5. Deal with head motion > 3mm
HMTemp = (ExcessiveHM==1)*1;
WantedSubMatrix=WantedSubMatrix .* HMTemp;

%%% Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);


WinLength = 80; WinStep = 5; WinType = 'rectangle';
if ~exist(OutDir, 'dir'); mkdir(OutDir); end
cd(DataDir)
tic
MissSub = [];
SubID = dir([DataDir, '/sub*']);
for sn = 1:numel(SubID)
    SubDir = [DataDir, '/' SubID(sn).name];
    if exist(SubDir, 'dir') 
%         f_DFC_W(SubDir, tpl_path, OutDir, MaskFile, WinLength, WinStep, WinType);
%         f_DFC_W_sect(SubDir, tpl_path, OutDir, MaskFile, WinLength, WinStep, WinType, 1:350);
        f_DFC_W_sect(SubDir, tpl_path, OutDir, MaskFile, WinLength, WinStep, WinType, 376:725);
        fprintf('%s finished, total elapsed time is %.2f mins, and there remains %d/%d subjects.\n', ...
            SubID(sn).name, toc/60, numel(SubID)-sn, numel(SubID))
    else
        fprintf('\n\n %s does not exist ============== \n', SubDir)
        MissSub = [MissSub, sn];
    end
end
cd(OutDir)
disp('done------done------done------done------done------done');
