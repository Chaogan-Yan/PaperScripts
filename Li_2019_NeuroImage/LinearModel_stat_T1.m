VarPath = '/mnt/Data/RfMRILab/Lile/DFCnew/SWU/subinfo_SWU.xlsx';
DataDir='/mnt/Data/RfMRILab/Lile/DFCnew/SWUv2v/KCC_gmd20_S1';
OutDir='/mnt/Data/RfMRILab/Lile/DFCnew/SWUv2v/T1_gmd20_S1';
MaskFile = '/mnt/Data/RfMRILab/Lile/DFCnew/SWUv2v/GroupMask_gmd20.nii';

VarAll = xlsread(VarPath, 1, 'b2:f1000');
[~, SubID] = xlsread(VarPath, 1, 'a2:a1000');

Sex = VarAll(:,1);
Age = VarAll(:,2);
FD = VarAll(:,3);       % session1
% FD = VarAll(:,11);      % session2
VoxelNum = VarAll(:,4);
QCscore = VarAll(:,5);
ExcessiveHM = VarAll(:,6);
WantedSubMatrix = ones(length(FD), 1);

%%% Linear Model: 
AllCov = [ones(length(SubID),1), Sex, Age, FD];
%Centering: Let the first column (constant) have the mean effect.
AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
Contrast=zeros(1,size(AllCov,2));
Contrast(1)=1;

FileList=[];
[Mask, ~] = y_Read(MaskFile);
ImgSum = zeros(size(Mask));
for iSub=1:length(SubID)
    FileList{iSub,1}=[DataDir, '/zW_', SubID{iSub},'.nii'];
    [ImgSub, ~] = y_Read(FileList{iSub});
    ImgSum = ImgSum+ImgSub;
end

OutputName=[OutDir, '/T1.nii'];
if ~exist(fileparts(OutputName), 'dir'); mkdir(fileparts(OutputName)); end
[b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(FileList,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
y_Write(ImgSum/iSub, Header, [OutputName(1:end-4), '_mean.nii'])
disp('done------done------done------done------done------done');
