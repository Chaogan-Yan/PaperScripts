

% X=[ones(size(DurationScoreBinarized)),DurationScoreBinarized,Age,Sex,Edu,Motion];
% Z={ones(size(DurationScoreBinarized)),DurationScoreBinarized};
% G={Site,Site};


% Convert into 2D
YMatrix=reshape(NetworkCorrSet,[],size(NetworkCorrSet,3))';

SubIDAll=SubID;
DxAll=Dx;
AgeAll=Age;
SexAll=Sex;
EduAll=Edu;
MotionAll=Motion;
DurationScoreBinarizedAll=DurationScoreBinarized;

SiteIndex = unique(Site);

TMatrix=zeros(length(SiteIndex),size(YMatrix,2));
N1=zeros(length(SiteIndex),1);
N2=zeros(length(SiteIndex),1);
for i=1:length(SiteIndex)
    SubID=SubIDAll(find(Site==SiteIndex(i)));
    Dx=DxAll(find(Site==SiteIndex(i)));
    Age=AgeAll(find(Site==SiteIndex(i)));
    Sex=SexAll(find(Site==SiteIndex(i)));
    Edu=EduAll(find(Site==SiteIndex(i)));
    Motion=MotionAll(find(Site==SiteIndex(i)));
    DurationScoreBinarized=DurationScoreBinarizedAll(find(Site==SiteIndex(i)));
    AllCov = [ones(length(SubID),1), DurationScoreBinarized,Age,Sex,Edu,Motion];
    %Centering: Let the first column (constant) have the mean effect.
    AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
    Contrast=zeros(1,size(AllCov,2));
    Contrast(2)=1;

    for iY=1:size(YMatrix,2)
        y=squeeze(YMatrix(find(Site==SiteIndex(i)),iY));
        [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(y,AllCov,Contrast,'T');
        TMatrix(i,iY)=TF_ForContrast;
    end
    N1(i,1)=length(find(DurationScoreBinarized==1));
    N2(i,1)=length(find(DurationScoreBinarized==-1));
end

TVal=TMatrix;
nTests=size(YMatrix,2);
OutputName='/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_848_794/Network/Edge/Meta/Temp/Temp';
[Path, fileN, extn] = fileparts(OutputName);
MatNameForR=fullfile(Path,[fileN,'_ForR.mat']);
save(MatNameForR,'TVal','N1','N2','nTests');
MatNameRResults=fullfile(Path,[fileN,'_RResults.mat']);

[ProgramPath] = fileparts(which('y_Meta_Image_CallR.m'));
Expression = sprintf('!Rscript %s%sR_Cal_Meta.R %s %s', ProgramPath, filesep,MatNameForR,MatNameRResults);
eval(Expression);
load(MatNameRResults);

ZMatrix=reshape(Z,[size(NetworkCorrSet,1) size(NetworkCorrSet,2)]);
PMatrix=reshape(P,[size(NetworkCorrSet,1) size(NetworkCorrSet,2)]);





%FDR
TriuMat = triu(ones(size(PMatrix)),0)';

PVector = PMatrix(find(TriuMat));

addpath /mnt/Data/RfMRILab/Yan/YAN_Program/gretna

[pID,pN] = FDR(PVector,0.05);

PSig = PMatrix<=pID;




