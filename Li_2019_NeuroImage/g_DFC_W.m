VarPath = '/mnt/Data/RfMRILab/Lile/DFCnew/SWU/subinfo_SWU.xlsx';
DataDir='/mnt/Data/RfMRILab/Lile/DFC/data_CORR/S2_FunImgARCWFS';
OutDir='/mnt/Data/RfMRILab/Lile/940/SWU_gmd20/KCC_S2';
MaskFile = '/mnt/Data/RfMRILab/Lile/DFCnew/SWU/GroupMask_gmd20.nii';
tpl_path = '/mnt/Data/RfMRILab/Lile/DFC/CC200ROI_333.nii';

VarAll = xlsread(VarPath, 1, 'b2:f1000');
[~, SubID] = xlsread(VarPath, 1, 'a2:a1000');

Sex = VarAll(:,1);
Age = VarAll(:,2);
FD = VarAll(:,3);
VoxelNum = VarAll(:,4);
QCscore = VarAll(:,5);
WantedSubMatrix = ones(length(SubID), 1);

%%%1. Deal With Age
AgeRange = [16 32];
AgeTemp = (Age>=AgeRange(1)) .* (Age<=AgeRange(2));
WantedSubMatrix = WantedSubMatrix.*AgeTemp;

%%%2. Deal with head motion FD
FDTemp = (FD<(mean(FD)+2*std(FD)))*1;
% WantedSubMatrix = WantedSubMatrix.*FDTemp;

%%%3. Deal with bad coverage
VoxelNumTemp = (VoxelNum>60000)*1;
WantedSubMatrix=WantedSubMatrix .* VoxelNumTemp;

%%%4. Deal with low QC score
QCTemp = (QCscore>=4)*1;
WantedSubMatrix=WantedSubMatrix .* QCTemp;

%%% Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Age=Age(WantedSubIndex);
FD=FD(WantedSubIndex);


WinLength = 32; WinStep = 2; WinType = 'rectwin';
if ~exist(OutDir, 'dir'); mkdir(OutDir); end
cd(DataDir)
tic
for sn = 1:numel(SubID)
    SubDir = [DataDir, '/' SubID{sn}];
    if exist(SubDir, 'dir') 
        f_DFC_W(SubDir, '', OutDir, MaskFile, WinLength, WinStep, WinType);
        fprintf('%s finished, total elapsed time is %.2f mins, and there remains %d/%d subjects.\n', ...
            SubID{sn}, toc/60, numel(SubID)-sn, numel(SubID))
    else
        fprintf('\n\n %s does not exist ============== \n', SubDir)
    end
end
cd(OutDir)
disp('done------done------done------done------done------done');
