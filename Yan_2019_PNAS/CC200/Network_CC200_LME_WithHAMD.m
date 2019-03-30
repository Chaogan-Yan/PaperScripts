
clear all;clc

%Select subjects
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_943_846/DataQC/CorrSet.mat
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/Info_Final1789_943_846.mat
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/NetworkAnalysis/zROICorr/ROISignals_FunImgARCWF/Crad200_ROICorrelation_FisherZ_Set.mat
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
ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);


load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/FirstEpisodeDrugNaive.mat
%Get the correponding FEDN scroe
FirstEpisodeScore=zeros(length(SubID),1);
DrugUseScore=zeros(length(SubID),1);
for i=1:length(SubID)
    for j=1:size(FEDNTalbe,1)
        if strcmpi(SubID{i},FEDNTalbe{j,1})
            FirstEpisodeScore(i)=FirstEpisode(j);
            DrugUseScore(i)=DrugUse(j);
        end
    end
end

%Correlation with HAMD
%load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/HAMD17.mat
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/HAMD17_Updated180903.mat

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
FirstEpisodeScore=FirstEpisodeScore(WantedSubIndex);
DrugUseScore=DrugUseScore(WantedSubIndex);
HAMD17Score=HAMD17Score(WantedSubIndex);
ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);


WantedSubMatrix=ones(length(SubID),1);
%%% FEDN
WantedSubMatrix(find( (Dx==1) .* ((FirstEpisodeScore==1).*((DrugUseScore==-1))==0) ))=0;


% %%% Recurrent
% WantedSubMatrix(find( (Dx==1) .* ((FirstEpisodeScore==-1)==0) ))=0; %%Recurrent

%Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);
FirstEpisodeScore=FirstEpisodeScore(WantedSubIndex);
DrugUseScore=DrugUseScore(WantedSubIndex);
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
FirstEpisodeScore=FirstEpisodeScore(WantedSubIndex);
DrugUseScore=DrugUseScore(WantedSubIndex);
HAMD17Score=HAMD17Score(WantedSubIndex);
ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);


%Check for networks
load /mnt/Data/RfMRILab/Lile/MDD/DMN/CC200_To_Yeo7_Label.mat


X=[ones(size(HAMD17Score)),HAMD17Score,Age,Sex,Edu,Motion];
Z={ones(size(HAMD17Score)),HAMD17Score};
G={Site,Site};



%Average FC as between network
%For Yeo

MergeLabel = YeoNetwork;
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
            NetworkCorr(j,k) = sum(CorrZ(MatrixIndex), 'omitnan');
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
        fprintf('%d NaN \n', length(find(isnan(y))) )
        lme = fitlmematrix(X,y,Z,G);
        TMatrix(ii,jj)=lme.Coefficients{2,4}; %
        PMatrix(ii,jj)=lme.Coefficients{2,6};
    end
    fprintf('%d\n\n',ii)
end


%FDR
TriuMat = triu(ones(size(PMatrix)),0)';

PVector = PMatrix(find(TriuMat));

addpath /mnt/Data/RfMRILab/Yan/YAN_Program/gretna

[pID,pN] = FDR(PVector,0.05);

if ~isempty(pID); PSig = PMatrix<=pID; end
