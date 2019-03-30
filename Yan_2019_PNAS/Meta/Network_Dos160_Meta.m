
addpath(fullfile(matlabroot,'toolbox','stats','stats'))

%Select subjects
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_943_846/DataQC/CorrSet.mat
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/Info_Final1789_943_846.mat

ReHoGood = (CorrSet_All(:,3) >= 0.6); %Exclude ReHo Correlation < 0.6

%Exclude Site with N<10
SubjectNumberPerSite=[];
SiteIndex = unique(Site);
WantedSubMatrix=ones(length(SubID),1);
for i=1:length(SiteIndex)
    DxTemp=Dx(find((Site==SiteIndex(i)).*ReHoGood)); %DxTemp=Dx(find(Site==SiteIndex(i)));
    SubjectNumberPerSite(i,:)=[SiteIndex(i),length(find(DxTemp==1)),length(find(DxTemp==-1))];
    if (length(find(DxTemp==1))<10)||(length(find(DxTemp==-1))<10)
        WantedSubMatrix(find(Site==SiteIndex(i)))=0;
    end
end

WantedSubMatrix = WantedSubMatrix.*ReHoGood;

%Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);


load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/NetworkAnalysis/zROICorr/ROISignals_FunImgARCWF/Dos160_ROICorrelation_FisherZ_Set.mat
ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);



%Check for networks
load /mnt/Data/RfMRILab/Yan/YAN_Program/Atlas/Dos160_WithName.mat
Network=zeros(160,1);
for i=1:length(Dos160_WithName)
    switch Dos160_WithName{i,5}
        case 'occipital'
            Network(i)=1;
        case 'sensorimotor'
            Network(i)=2;
        case 'default'
            Network(i)=3;
        case 'fronto-parietal'
            Network(i)=4;
        case 'cingulo-opercular'
            Network(i)=5;
        case 'cerebellum'
            Network(i)=6;
    end
end

Network142=Network;
Network142(find(Network==6))=[];
Dos142_WithName=Dos160_WithName;
Dos142_WithName(find(Network==6),:)=[];

Network142_Yeo=Dos160_YeoNetwork_YanModified; %Network142_Yeo=Dos160_YeoNetwork;
Network142_Yeo(find(Network==6))=[];

Network142_Yeo(find(Network142_Yeo==10))=5; %Change subcortical to 5





ROICorrelation_FisherZ_Set(find(Network==6),:,:)=[];
ROICorrelation_FisherZ_Set(:,find(Network==6),:)=[];



X=[ones(size(Dx)),Dx,Age,Sex,Edu,Motion];
Z={ones(size(Dx)),Dx};
G={Site,Site};









%Average FC as between network
%For Yeo
MergeLabel = Network142_Yeo;
LabelIndex = unique(MergeLabel);
NetworkCorr = zeros(length(LabelIndex),length(LabelIndex));
NetworkCorrSet = zeros(length(LabelIndex),length(LabelIndex),size(ROICorrelation_FisherZ_Set,3));
FullMatrix=ones(size(ROICorrelation_FisherZ_Set,1),size(ROICorrelation_FisherZ_Set,2))-eye(size(ROICorrelation_FisherZ_Set,1),size(ROICorrelation_FisherZ_Set,2));

for iSub=1:size(ROICorrelation_FisherZ_Set,3)
    CountSet_Full = zeros(length(LabelIndex),length(LabelIndex));
    for j=1:length(LabelIndex)
        for k=1:length(LabelIndex)
            A=double(MergeLabel==LabelIndex(j));
            B=double(MergeLabel==LabelIndex(k));
            Matrix = A*B';
            MatrixIndex = find(Matrix);
            CorrZ = ROICorrelation_FisherZ_Set(:,:,iSub);
            NetworkCorr(j,k) = sum(CorrZ(MatrixIndex));
            CountSet_Full(j,k) = sum(FullMatrix(MatrixIndex));
        end
    end
    NetworkCorr=NetworkCorr./CountSet_Full;
    NetworkCorrSet(:,:,iSub)=NetworkCorr;
end





% Convert into 2D
YMatrix=reshape(NetworkCorrSet,[],size(NetworkCorrSet,3))';

SubIDAll=SubID;
DxAll=Dx;
AgeAll=Age;
SexAll=Sex;
EduAll=Edu;
MotionAll=Motion;

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
    AllCov = [ones(length(SubID),1), Dx,Age,Sex,Edu,Motion];
    %Centering: Let the first column (constant) have the mean effect.
    AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
    Contrast=zeros(1,size(AllCov,2));
    Contrast(2)=1;

    for iY=1:size(YMatrix,2)
        y=squeeze(YMatrix(find(Site==SiteIndex(i)),iY));
        [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(y,AllCov,Contrast,'T');
        TMatrix(i,iY)=TF_ForContrast;
    end
    N1(i,1)=length(find(Dx==1));
    N2(i,1)=length(find(Dx==-1));
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





















%Nodal Degree
TMatrix = zeros(size(ROICorrelation_FisherZ_Set,2),1);
PMatrix = ones(size(ROICorrelation_FisherZ_Set,2),1);

%DegreeSet=sum(ROICorrelation_FisherZ_Set,1);
DegreeSet=sum(ROICorrelation_FisherZ_Set.*(ROICorrelation_FisherZ_Set>0.25),1);

for ii=1:size(DegreeSet,2)
    y=squeeze(DegreeSet(1,ii,:));
    lme = fitlmematrix(X,y,Z,G);
    TMatrix(ii)=lme.Coefficients{2,4}; %
    PMatrix(ii)=lme.Coefficients{2,6};
end



PVector = PMatrix;

addpath /mnt/Data/RfMRILab/Yan/YAN_Program/gretna

[pID,pN] = FDR(PVector,0.05);

PSig = PMatrix<=pID;







%For edges

% Convert into 2D
YMatrix=reshape(ROICorrelation_FisherZ_Set,[],size(ROICorrelation_FisherZ_Set,3))';

SubIDAll=SubID;
DxAll=Dx;
AgeAll=Age;
SexAll=Sex;
EduAll=Edu;
MotionAll=Motion;

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
    AllCov = [ones(length(SubID),1), Dx,Age,Sex,Edu,Motion];
    %Centering: Let the first column (constant) have the mean effect.
    AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
    Contrast=zeros(1,size(AllCov,2));
    Contrast(2)=1;

    for iY=1:size(YMatrix,2)
        y=squeeze(YMatrix(find(Site==SiteIndex(i)),iY));
        [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(y,AllCov,Contrast,'T');
        TMatrix(i,iY)=TF_ForContrast;
    end
    N1(i,1)=length(find(Dx==1));
    N2(i,1)=length(find(Dx==-1));
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

ZMatrix=reshape(Z,[size(ROICorrelation_FisherZ_Set,1) size(ROICorrelation_FisherZ_Set,2)]);
PMatrix=reshape(P,[size(ROICorrelation_FisherZ_Set,1) size(ROICorrelation_FisherZ_Set,2)]);









%Check for networks
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_848_794/Network/Edge/Meta/Edge_Meta.mat

%FDR
TriuMat = triu(ones(size(PMatrix)),1)';

PVector = PMatrix(find(TriuMat));

addpath /mnt/Data/RfMRILab/Yan/YAN_Program/gretna

[pID,pN] = FDR(PVector,0.05);

PSig = PMatrix<=pID;

MergeLabel = Network142_Yeo;

PSurviveP = (PMatrix<=pID).*(ZMatrix>0);
PSurviveN = (PMatrix<=pID).*(ZMatrix<0);

PSurviveCount=PSurviveP;

LabelIndex = unique(MergeLabel);
CountSet = zeros(length(LabelIndex),length(LabelIndex));
CountSet_Full = zeros(length(LabelIndex),length(LabelIndex));
FullMatrix=ones(size(PSurviveCount))-eye(size(PSurviveCount));
for j=1:length(LabelIndex)
    for k=1:length(LabelIndex)
        A=double(MergeLabel==LabelIndex(j));
        B=double(MergeLabel==LabelIndex(k));
        Matrix = A*B';
        MatrixIndex = find(Matrix);
        CountSet(j,k) = sum(PSurviveCount(MatrixIndex));
        CountSet_Full(j,k) = sum(FullMatrix(MatrixIndex));
    end
end

CountSet=CountSet./(eye(size(CountSet))+ones(size(CountSet)));
CountSet_Full=CountSet_Full./(eye(size(CountSet_Full))+ones(size(CountSet_Full)));
CountSetPercent=CountSet./CountSet_Full;





%Get connection name
Table=[];
for j=1:length(LabelIndex)
    for k=j:length(LabelIndex)
        A=double(MergeLabel==LabelIndex(j));
        B=double(MergeLabel==LabelIndex(k));
        Matrix = A*B';
        
        PSurviveCountMatrix = PSurviveCount.*Matrix;
        PSurviveCountMatrixIndex = find(PSurviveCountMatrix);

        for iInd=1:length(PSurviveCountMatrixIndex)
            [II,JJ] = ind2sub(size(PSurviveCount),PSurviveCountMatrixIndex(iInd));
            if ~(j==k && JJ>=II)
                Row={ZMatrix(II,JJ),j,k,II,JJ,Dos142_WithName{II,4},Dos142_WithName{JJ,4}};
                Table=[Table;Row]
            end
        end
    end
end






%Check for Systems
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_848_794/Network/Edge/Edge_LME.mat

%Restore to symetric
TMatrix=TMatrix+TMatrix';
TriuMat = triu(ones(size(PMatrix)),1);
PMatrix(find(TriuMat))=0;
PMatrix=PMatrix+PMatrix';
PMatrix(find(eye(size(PMatrix))))=1;

MergeLabel = Network142;

PSurviveP = (PMatrix<=pID).*(TMatrix>0);
PSurviveN = (PMatrix<=pID).*(TMatrix<0);

PSurviveCount=PSurviveP;

LabelIndex = unique(MergeLabel);

CountSet = zeros(length(LabelIndex),length(LabelIndex));
CountSet_Full = zeros(length(LabelIndex),length(LabelIndex));
FullMatrix=ones(size(PSurviveCount))-eye(size(PSurviveCount));
for j=1:length(LabelIndex)
    for k=1:length(LabelIndex)
        A=double(MergeLabel==LabelIndex(j));
        B=double(MergeLabel==LabelIndex(k));
        Matrix = A*B';
        MatrixIndex = find(Matrix);
        CountSet(j,k) = sum(PSurviveCount(MatrixIndex));
        CountSet_Full(j,k) = sum(FullMatrix(MatrixIndex));
    end
end

CountSet=CountSet./(eye(size(CountSet))+ones(size(CountSet)));
CountSet_Full=CountSet_Full./(eye(size(CountSet_Full))+ones(size(CountSet_Full)));
CountSetPercent=CountSet./CountSet_Full;












%Test if site random effect significant
X=[ones(size(Dx)),Dx,Age,Sex,Edu,Motion];
ZReduced={ones(size(Dx))};
GReduced={Site};

NetworkCorrSet(5,:,:)=[];
NetworkCorrSet(:,5,:)=[];
TMatrix = zeros(size(NetworkCorrSet,2),size(NetworkCorrSet,2));
PMatrix = ones(size(NetworkCorrSet,2),size(NetworkCorrSet,2));
for ii=1:size(NetworkCorrSet,2)
    for jj=1:size(NetworkCorrSet,2)
        y=squeeze(NetworkCorrSet(ii,jj,:));
        lme = fitlmematrix(X,y,Z,G);
        lme_Reduced = fitlmematrix(X,y,ZReduced,GReduced);
        ModelCompare=compare(lme_Reduced,lme);
        PMatrix(ii,jj)=ModelCompare.pValue(2);
    end
    fprintf('%d\n',ii)
end






%Check Site Effect

SubIDAll=SubID;
DxAll=Dx;
AgeAll=Age;
SexAll=Sex;
EduAll=Edu;
MotionAll=Motion;

SiteIndex = unique(Site);
TMatrix = zeros(size(NetworkCorrSet,2),size(NetworkCorrSet,2),length(SiteIndex));
PMatrix = ones(size(NetworkCorrSet,2),size(NetworkCorrSet,2),length(SiteIndex));
for i=1:length(SiteIndex)
    
    SubID=SubIDAll(find(Site==SiteIndex(i)));
    Dx=DxAll(find(Site==SiteIndex(i)));
    Age=AgeAll(find(Site==SiteIndex(i)));
    Sex=SexAll(find(Site==SiteIndex(i)));
    Edu=EduAll(find(Site==SiteIndex(i)));
    Motion=MotionAll(find(Site==SiteIndex(i)));
    
    
    AllCov = [ones(length(SubID),1), Dx,Age,Sex,Edu,Motion];
    
    %Centering: Let the first column (constant) have the mean effect.
    AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
    
    Contrast=zeros(1,size(AllCov,2));
    Contrast(2)=1;
    
    for ii=1:size(NetworkCorrSet,2)
        for jj=1:size(NetworkCorrSet,2)
            y=squeeze(NetworkCorrSet(ii,jj,find(Site==SiteIndex(i))));
            [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,AllCov,Contrast,'T');
            
            TMatrix(ii,jj,i)= T(2); %
            PMatrix(ii,jj,i)=2*(1-tcdf(abs(T(2)), size(AllCov,1)-size(AllCov,2)));
        end
        fprintf('%d\n',ii)
    end
    
end







%Correlation with HAMD
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/HAMD17.mat

%Get the correponding HAMD scroe
HAMD17Score=zeros(length(SubID),1);
for i=1:length(SubID)
    for j=1:size(HAMD17,1)
        if strcmpi(SubID{i},HAMD17{j,1})
            if isempty(HAMD17{j,2})
                HAMD17Score(i)=0;
            else
                HAMD17Score(i)=HAMD17{j,2};
            end
        end
    end
end

WantedSubMatrix=ones(length(SubID),1);
WantedSubMatrix(find(HAMD17Score<7))=0;
%Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);
HAMD17Score=HAMD17Score(WantedSubIndex);

ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);




%Exclude Site with N<10
SubjectNumberPerSite=[];
SiteIndex = unique(Site);
WantedSubMatrix=ones(length(SubID),1);
for i=1:length(SiteIndex)
    DxTemp=Dx(find((Site==SiteIndex(i)))); %DxTemp=Dx(find(Site==SiteIndex(i)));
    SubjectNumberPerSite(i,:)=[SiteIndex(i),length(find(DxTemp==1)),length(find(DxTemp==-1))];
    if length(find(DxTemp==1))<10
        WantedSubMatrix(find(Site==SiteIndex(i)))=0;
    end
end

%Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);
HAMD17Score=HAMD17Score(WantedSubIndex);


ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);



X=[ones(size(HAMD17Score)),HAMD17Score,Age,Sex,Edu,Motion];
Z={ones(size(HAMD17Score)),HAMD17Score};
G={Site,Site};















%Correlation with Illness Druation
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/IllnessDuration.mat

%Get the correponding HAMD scroe
DurationScore=zeros(length(SubID),1);
for i=1:length(SubID)
    for j=1:size(IllnessDuration,1)
        if strcmpi(SubID{i},IllnessDuration{j,1})
            if isempty(IllnessDuration{j,2})
                DurationScore(i)=0;
            else
                DurationScore(i)=IllnessDuration{j,2};
            end
        end
    end
end

WantedSubMatrix=ones(length(SubID),1);
WantedSubMatrix(find(DurationScore<0.1))=0;
%Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);
DurationScore=DurationScore(WantedSubIndex);

ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);




%Exclude Site with N<10
SubjectNumberPerSite=[];
SiteIndex = unique(Site);
WantedSubMatrix=ones(length(SubID),1);
for i=1:length(SiteIndex)
    DxTemp=Dx(find((Site==SiteIndex(i)))); %DxTemp=Dx(find(Site==SiteIndex(i)));
    SubjectNumberPerSite(i,:)=[SiteIndex(i),length(find(DxTemp==1)),length(find(DxTemp==-1))];
    if length(find(DxTemp==1))<10
        WantedSubMatrix(find(Site==SiteIndex(i)))=0;
    end
end

%Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);
DurationScore=DurationScore(WantedSubIndex);


ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);



X=[ones(size(DurationScore)),DurationScore,Age,Sex,Edu,Motion];
Z={ones(size(DurationScore)),DurationScore};
G={Site,Site};





