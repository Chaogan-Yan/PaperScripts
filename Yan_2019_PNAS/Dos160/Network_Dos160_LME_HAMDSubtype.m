clear all
clc


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







%HAMD Subtype
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/HAMD17_Updated180903.mat

HAMD17SubItem{678,10}=4; % Original value is 10, which is impossible

HAMD17SubNum=HAMD17SubItem(:,2:end);
for i=1:size(HAMD17SubNum,1)
    for j=1:size(HAMD17SubNum,2)
        if isempty(HAMD17SubNum{i,j})
            HAMD17SubNum{i,j}=0;
        end
    end
end
HAMD17SubNum=cell2mat(HAMD17SubNum);

CD=-1*ones(size(HAMD17SubItem,1),1);
ANX=-1*ones(size(HAMD17SubItem,1),1);
NVSM=-1*ones(size(HAMD17SubItem,1),1);

CD(find( ((HAMD17SubNum(:,1)>=3) .* (HAMD17SubNum(:,7)>=3)) )) = 1;
ANX(find( sum(HAMD17SubNum(:,[9,10,11,15]),2) >= 6)) = 1; 
NVSM(find( ((HAMD17SubNum(:,6)==1)+(HAMD17SubNum(:,6)==2)) .* ((HAMD17SubNum(:,12)==1)+(HAMD17SubNum(:,12)==2)) )) = 1;

CD(find(sum(HAMD17SubNum,2)<7))=0;
ANX(find(sum(HAMD17SubNum,2)<7))=0;
NVSM(find(sum(HAMD17SubNum,2)<7))=0;


%Subtype1300 = CD; %Test for subtype CD+ vs. CD-
%Subtype1300 = ANX; %Test for subtype ANX+ vs. ANX-
%Subtype1300 = NVSM; %Test for subtype NVSM+ vs. NVSM-

% %CD vs ANX
% CDnoANX=CD;CDnoANX(find(ANX==1))=0;
% ANXnoCD=ANX;ANXnoCD(find(CD==1))=0;
% Subtype1300=zeros(1300,1);
% Subtype1300(find(CDnoANX))=1;
% Subtype1300(find(ANXnoCD))=-1;

% %CD vs NVSM
% CDnoNVSM=CD;CDnoNVSM(find(NVSM==1))=0;
% NVSMnoCD=NVSM;NVSMnoCD(find(CD==1))=0;
% Subtype1300=zeros(1300,1);
% Subtype1300(find(CDnoNVSM))=1;
% Subtype1300(find(NVSMnoCD))=-1;

%ANX vs NVSM
ANXnoNVSM=ANX;ANXnoNVSM(find(NVSM==1))=0;
NVSMnoANX=NVSM;NVSMnoANX(find(ANX==1))=0;
Subtype1300=zeros(1300,1);
Subtype1300(find(ANXnoNVSM))=1;
Subtype1300(find(NVSMnoANX))=-1;


%Get the correponding HAMD scroe
Subtype=zeros(length(SubID),1);
for i=1:length(SubID)
    for j=1:size(HAMD17SubItem,1)
        if strcmpi(SubID{i},HAMD17SubItem{j,1})
            Subtype(i)=Subtype1300(j);
        end
    end
end

WantedSubMatrix=ones(length(SubID),1);
WantedSubMatrix(find(Subtype==0))=0;
%Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);
Subtype=Subtype(WantedSubIndex);

ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);



%Exclude Site with N<10
SubjectNumberPerSite=[];
SiteIndex = unique(Site);
WantedSubMatrix=ones(length(SubID),1);

for i=1:length(SiteIndex)
    SubtypeTemp=Subtype(find((Site==SiteIndex(i)))); %DxTemp=Dx(find(Site==SiteIndex(i)));
    SubjectNumberPerSite(i,:)=[SiteIndex(i),length(find(SubtypeTemp==1)),length(find(SubtypeTemp==-1))];
    if (length(find(SubtypeTemp==1))<10)||(length(find(SubtypeTemp==-1))<10)
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
Subtype=Subtype(WantedSubIndex);


ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);



X=[ones(size(Subtype)),Subtype,Age,Sex,Edu,Motion];
Z={ones(size(Subtype)),Subtype};
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




%For Average FC
TMatrix = zeros(size(NetworkCorrSet,2),size(NetworkCorrSet,2));
PMatrix = ones(size(NetworkCorrSet,2),size(NetworkCorrSet,2));
for ii=1:size(NetworkCorrSet,2)
    for jj=1:size(NetworkCorrSet,2)
        y=squeeze(NetworkCorrSet(ii,jj,:));
        lme = fitlmematrix(X,y,Z,G);
        TMatrix(ii,jj)=lme.Coefficients{2,4}; %
        PMatrix(ii,jj)=lme.Coefficients{2,6};
    end
    fprintf('%d\n',ii)
end


%FDR
TriuMat = triu(ones(size(PMatrix)),0)';

PVector = PMatrix(find(TriuMat));

addpath /mnt/Data/RfMRILab/Yan/YAN_Program/gretna

[pID,pN] = FDR(PVector,0.05);

PSig = PMatrix<=pID;

%edited by ChenXiao
%To extract the phenotype Score










