clear
VarPath = '/mnt/Data/RfMRILab/Lile/DFCnew/movie/movie_subinfo.xlsx';
DataDirP1 = '/mnt/Data/RfMRILab/Lile/DFCnew/HBNv2v/KCC_MV1';
DataDirP2 = '/mnt/Data/RfMRILab/Lile/DFCnew/HBNv2v/KCC_RS1';
OutputName= '/mnt/Data/RfMRILab/Lile/DFCnew/HBNv2v/PT_gmd20/PT_MV1-RS1.nii';

MaskFile = '/mnt/Data/RfMRILab/Lile/DFCnew/HBNv2v/GroupMask_gmd20.nii';

VarAll = xlsread(VarPath, 1, 'd:j');
[~, SubID] = xlsread(VarPath, 1, 'b2:b50');

Sex = VarAll(:,1);
Age = VarAll(:,2);
FD = VarAll(:,3);
FD2 = VarAll(:,5);

%%% Linear Model: 
nSub = length(SubID);
AllCov = [ones(nSub,1);-1*ones(nSub,1)];
for i=1:nSub
    SubjectRegressors(:,i) = zeros(nSub*2,1);
    SubjectRegressors(i:nSub:nSub*2,i) = 1;
end
Cov = [FD; FD2];
AllCov = [AllCov,SubjectRegressors, Cov];

Contrast = zeros(1,size(AllCov,2));
Contrast(1) = 1;

FileList=[];
for iSub=1:length(SubID)
    FileList{iSub,1}=[DataDirP1, '/zW_', SubID{iSub},'.nii'];
    FileList{iSub+nSub,1}=[DataDirP2, '/zW_', SubID{iSub},'.nii'];
end

if ~exist(fileparts(OutputName), 'dir'); mkdir(fileparts(OutputName)); end
[b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = ...
    y_GroupAnalysis_Image(FileList,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
disp('done------done------done------done------done------done');
