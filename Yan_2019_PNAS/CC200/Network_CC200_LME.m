
clear all;clc

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

load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/IllnessDuration.mat
%Get the correponding IllnessDuration scroe
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

load /mnt/Data/RfMRILab/Lile/MDD/DMN/Episode.mat
EpisodeNum=zeros(length(SubID),1);
for i=1:length(SubID)
    for j=1:size(Episode,1)
        if strcmpi(SubID{i},Episode{j,1})
            if isempty(Episode{j,2})
                EpisodeNum(i)=0;
            else
                EpisodeNum(i)=Episode{j,2};
            end
        end
    end
end

% load /mnt/Data/RfMRILab/Lile/MDD/DMN/HAMD17.mat
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/HAMD17_Updated180903.mat

HAMD=zeros(length(SubID),1);
for i=1:length(SubID)
    for j=1:size(HAMD17,1)
        if strcmpi(SubID{i},HAMD17{j,1})
            if isempty(HAMD17{j,2})
                HAMD(i)=0;
            else
                HAMD(i)=HAMD17{j,2};
            end
        end
    end
end

% %%% FEDN vs NC
% WantedSubMatrix(find( (Dx==1) .* ((FirstEpisodeScore==1).*((DrugUseScore==-1))==0) ))=0;


% %%% Recurrent vs NC
% WantedSubMatrix(find( (Dx==1) .* ((FirstEpisodeScore==-1)==0) ))=0; %%Recurrent


% %%% FEDN vs Recurrent
% WantedSubMatrix(find( ((FirstEpisodeScore==1).*((DrugUseScore==-1))+(FirstEpisodeScore==-1)) ==0 )) = 0;
% Dx = FirstEpisodeScore;


% %%% DrugUse vs NC
% WantedSubMatrix(DrugUseScore~=1 & DrugUseScore~=0) = 0; %%Exclude noDrugUse and unknown


% %%% FEDN vs DrugUse
% SelectedSub = (FirstEpisodeScore==1)&(DrugUseScore==-1) | (DrugUseScore==1);
% WantedSubMatrix( SelectedSub==0 ) =0;
% Dx = zeros(length(Dx),1);
% Dx( (FirstEpisodeScore==1)&(DrugUseScore==-1) ) = 1;
% Dx(DrugUseScore==1) = -1;


% %%% low vs high duration score
% WantedSubMatrix( DurationScore==0 | (DurationScore>6 & DurationScore<24) ) = 0;
% Dx = zeros(length(Dx),1);
% Dx(DurationScore<=6 & DurationScore>0) = 1;
% Dx(DurationScore>=24) = -1;


% %%% low vs high duration score within FEDN
% WantedSubMatrix( DurationScore==0 | (DurationScore>3 & DurationScore<12) ) = 0;
% WantedSubMatrix(find( (Dx==1) .* ((FirstEpisodeScore==1).*((DrugUseScore==-1))==0) ))=0;
% Dx = zeros(length(Dx),1);
% Dx(DurationScore<=3 & DurationScore>0) = 1;
% Dx(DurationScore>=12) = -1;


% AgeTemp = Age >30;
% WantedSubMatrix = WantedSubMatrix.*AgeTemp;
% SexTemp = Sex==2;
% WantedSubMatrix = WantedSubMatrix.*SexTemp;

%Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);
EpisodeNum=EpisodeNum(WantedSubIndex);
FirstEpisodeScore=FirstEpisodeScore(WantedSubIndex);
DrugUseScore=DrugUseScore(WantedSubIndex);
DurationScore=DurationScore(WantedSubIndex);
HAMD=HAMD(WantedSubIndex);


WantedSubMatrix1789=WantedSubMatrix;
WantedSubMatrix1789Index=find(WantedSubMatrix1789);

%Exclude Site with N<10
SubjectNumberPerSite=[];
SiteIndex = unique(Site);
WantedSubMatrix=ones(length(SubID),1);
for i=1:length(SiteIndex)
    DxTemp=Dx(find((Site==SiteIndex(i)))); %DxTemp=Dx(find(Site==SiteIndex(i)));
    SubjectNumberPerSite(i,:)=[SiteIndex(i),length(find(DxTemp==1)),length(find(DxTemp==-1))];
    if (length(find(DxTemp==1))<10)||(length(find(DxTemp==-1))<10)
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
EpisodeNum=EpisodeNum(WantedSubIndex);
FirstEpisodeScore=FirstEpisodeScore(WantedSubIndex);
DrugUseScore=DrugUseScore(WantedSubIndex);
DurationScore=DurationScore(WantedSubIndex);
HAMD=HAMD(WantedSubIndex);

WantedSubMatrix1789Index=WantedSubMatrix1789Index(WantedSubIndex);

WantedSubMatrix1789=zeros(1789,1);
WantedSubMatrix1789(WantedSubMatrix1789Index)=1;


load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/NetworkAnalysis/zROICorr/ROISignals_FunImgARCWF/Crad200_ROICorrelation_FisherZ_Set.mat
ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,find(WantedSubMatrix1789));




%Check for networks
load /mnt/Data/RfMRILab/Lile/MDD/DMN/CC200_To_Yeo7_Label.mat


X=[ones(size(Dx)),Dx,Age,Sex,Edu,Motion];
Z={ones(size(Dx)),Dx};
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


SiteID = unique(Site);
yy = zeros(length(y), 1);
for sn = 1:length(SiteID)
    idx = Site==SiteID(sn);
    cov = X(idx,:);
    cov(:,2:end) = (cov(:,2:end)-repmat(mean(cov(:,2:end)),size(cov(:,2:end),1),1));
    [b,r,SSE,SSR, T] = y_regress_ss(y(idx), cov);
    yy(idx) = y(idx) - cov(:,3:6)*b(3:6);
%     yy=y;
    
    [H,P,CI,STATS] = ttest2(yy(idx&Dx==1), yy(idx&Dx==-1));
    G1MeanSite(sn,1) = mean(yy(idx & Dx==1));
    G2MeanSite(sn,1) = mean(yy(idx & Dx==-1));
    tSite(sn,1) = STATS.tstat;
    pSite(sn,1) = P;
end
SubInfo = [yy, Dx, Site];
disp([SiteID, G1MeanSite, G2MeanSite, tSite, pSite])