%This is the final version
clear;clc;
addpath('/mnt/Data/RfMRILab/Yan/YAN_Program');
addpath('/mnt/Data/RfMRILab/Yan/YAN_Program/BCT_20120814');

%Main analysis, cc200 and scrubbing
%Main analysis
metricpath = '/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SmallWorld/SmallWorld/Dos160SmallWorld/AllSM_ReRunWeighted/GTA_WieghtedUndirected.mat';
statspath = '/mnt/Data/RfMRILab/ChenX/Small_World/Weighted_BugFree';
figspath = '/mnt/Data/RfMRILab/ChenX/Small_World/Figures';
load /mnt/Data/RfMRILab/ChenX/Small_World/SubLoc_1627_1789.mat
load /mnt/Data/RfMRILab/ChenX/Small_World/Info_Final1789_943_846.mat
SubIDAll = SubID(SubLoc_1627_1789);

% %scrubbing
% metricpath = '/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SmallWorld/SmallWorld/Dos160SmallWorld/AllSM_Scrubbing/GTA_WieghtedUndirected.mat';
% load /mnt/Data/RfMRILab/ChenX/Small_World/SubID_1625.mat
% SubIDAll = SubID_1625;

%Select subjects
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_943_846/DataQC/CorrSet.mat
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/Info_Final1789_943_846.mat
ReHoGood = (CorrSet_All(:,3) >= 0.6); %Exclude ReHo Correlation < 0.6
WantedSubMatrix=ones(length(SubID),1);

    
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

% %%%% MDD vs NC
% SavePath = '/mnt/Data/RfMRILab/ChenX/Small_World/Weighted_BugFree_scrubbing/MDD_NC/AucStat_MDD_NC.mat';
% mkdir('/mnt/Data/RfMRILab/ChenX/Small_World/Weighted_BugFree_scrubbing/MDD_NC');

% %%%% FEDN vs NC
% WantedSubMatrix(find( (Dx==1) .* ((FirstEpisodeScore==1).*((DrugUseScore==-1))==0) ))=0;
% SavePath = '/mnt/Data/RfMRILab/ChenX/Small_World/Weighted_BugFree_scrubbing/FEDN_NC/AucStat_FEDN_NC.mat';
% mkdir('/mnt/Data/RfMRILab/ChenX/Small_World/Weighted_BugFree_scrubbing/FEDN_NC');

% %%%% Recurrent vs NC
% WantedSubMatrix(find( (Dx==1) .* ((FirstEpisodeScore==-1)==0) ))=0;
% SavePath = '/mnt/Data/RfMRILab/ChenX/Small_World/Weighted_BugFree_scrubbing/Recurrent_NC/AucStat_Recurrent_NC.mat';
% mkdir('/mnt/Data/RfMRILab/ChenX/Small_World/Weighted_BugFree_scrubbing/Recurrent_NC');


%%% FEDN vs Recurrent
WantedSubMatrix(find( ((FirstEpisodeScore==1).*((DrugUseScore==-1))+(FirstEpisodeScore==-1)) ==0 )) = 0;
Dx = FirstEpisodeScore;
SavePath = [statspath,'/FEDN_Recurrent/AucStat_FEDN_Recurrent.mat'];
mkdir([statspath,'/FEDN_Recurrent']);
figSavePath = [figspath,'/FEDN_Recurrent'];

%Exclude ReHo
WantedSubMatrix = WantedSubMatrix.*ReHoGood;

%Exclude Site with N<10
SubjectNumberPerSite=[];
SiteIndex = unique(Site);
for i=1:length(SiteIndex)
    DxTemp=Dx(find((Site==SiteIndex(i)).*ReHoGood.*WantedSubMatrix)); %DxTemp=Dx(find(Site==SiteIndex(i)));
    SubjectNumberPerSite(i,:)=[SiteIndex(i),length(find(DxTemp==1)),length(find(DxTemp==-1))];
    if (length(find(DxTemp==1))<10)||(length(find(DxTemp==-1))<10)
        WantedSubMatrix(find(Site==SiteIndex(i)))=0;
    end
end

%Exclude Site 4
WantedSubMatrix(find(Site==4))=0;


%Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);
FirstEpisodeScore=FirstEpisodeScore(WantedSubIndex,:);
DrugUseScore=DrugUseScore(WantedSubIndex,:);


load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/NetworkAnalysis/zROICorr/ROISignals_FunImgARCWF/Dos160_ROICorrelation_FisherZ_Set.mat
ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);

%Check NaN
TriuMat = triu(ones(size(ROICorrelation_FisherZ_Set,1),size(ROICorrelation_FisherZ_Set,2)),1);
X=zeros(size(ROICorrelation_FisherZ_Set,3),1);
for i=1:size(ROICorrelation_FisherZ_Set,3)
    Temp=ROICorrelation_FisherZ_Set(:,:,i);
    X(i)=sum(Temp(find(TriuMat)));
end


%Get rid of those with NaNs
WantedSubMatrix=ones(length(SubID),1);
WantedSubMatrix(find(isnan(X)))=0;

%Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);
FirstEpisodeScore=FirstEpisodeScore(WantedSubIndex,:);
DrugUseScore=DrugUseScore(WantedSubIndex,:);

ROICorrelation_FisherZ_Set = ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);

% length(find((FirstEpisodeScore==-1).*(DrugUseScore==-1)))
% length(find((FirstEpisodeScore==-1).*(DrugUseScore==1)))
% length(find((FirstEpisodeScore==-1)))

% exclude those not belong to the small worldness matrix
[~, WantedSubIndex] = ismember(SubIDAll, SubID);
WantedSubIndex(find(WantedSubIndex == 0)) = [];
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);
FirstEpisodeScore=FirstEpisodeScore(WantedSubIndex,:);
DrugUseScore=DrugUseScore(WantedSubIndex,:);
ROICorrelation_FisherZ_Set = ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);



%1. load metrics' matrix
MAT = load(metricpath);
MAT.SparsityRange=[0.01:0.01:0.5]';
%MAT.SparsityRange = [0.02:0.02:1]';
DeltaSparsity=0.01;
AUCRange = [10:34];
Index=[1:size(MAT.CpSet,1)];
%SparsityRange=MAT.SparsityRange;
NodeNumber = size(MAT.DegreeSet,3);

% AllSet_Cp_AUC = zeros(length(ConditionList),1);
% AllSet_Lp_AUC = zeros(length(ConditionList),1);
% AllSet_Gamma_AUC = zeros(length(ConditionList),1);
% AllSet_Lambda_AUC = zeros(length(ConditionList),1);
% AllSet_Sigma_AUC = zeros(length(ConditionList),1);
% AllSet_Eloc_AUC = zeros(length(ConditionList),1);
% AllSet_Eglob_AUC = zeros(length(ConditionList),1);
% AllSet_Assortativity_AUC = zeros(length(ConditionList),1);
% AllSet_Modularity_AUC = zeros(length(ConditionList),1);

%%%AUC
[~, SubLoc] = ismember(SubID, SubIDAll);
SubLoc(find(SubLoc == 0)) = [];

AllSet_Degree_AUC = zeros(length(SubLoc),NodeNumber);
AllSet_NodalEfficiency_AUC = zeros(length(SubLoc),NodeNumber);
AllSet_Betweenness_AUC = zeros(length(SubLoc),NodeNumber);
AllSet_ClusteringCoefficient_AUC = zeros(length(SubLoc),NodeNumber);
AllSet_ParticipantCoefficient_AUC = zeros(length(SubLoc),NodeNumber);
AllSet_SubgraphCentrality_AUC = zeros(length(SubLoc),NodeNumber);
AllSet_EigenvectorCentrality_AUC = zeros(length(SubLoc),NodeNumber);
AllSet_PageRankCentrality_AUC = zeros(length(SubLoc),NodeNumber);

%Temp=(mean(MAT.CpSet(Index,AUCRange),2));
Temp=(sum(MAT.CpSet(Index,AUCRange),2) - sum(MAT.CpSet(Index,AUCRange([1 end])),2)/2)*DeltaSparsity;
AllSet_Cp_AUC=Temp(SubLoc);
Temp=(sum(MAT.LpSet(Index,AUCRange),2) - sum(MAT.LpSet(Index,AUCRange([1 end])),2)/2)*DeltaSparsity;
AllSet_Lp_AUC=Temp(SubLoc);
Temp=(sum(MAT.GammaSet(Index,AUCRange),2) - sum(MAT.GammaSet(Index,AUCRange([1 end])),2)/2)*DeltaSparsity;
AllSet_Gamma_AUC=Temp(SubLoc);
Temp=(sum(MAT.LambdaSet(Index,AUCRange),2) - sum(MAT.LambdaSet(Index,AUCRange([1 end])),2)/2)*DeltaSparsity;
AllSet_Lambda_AUC=Temp(SubLoc);
Temp=(sum(MAT.SigmaSet(Index,AUCRange),2) - sum(MAT.SigmaSet(Index,AUCRange([1 end])),2)/2)*DeltaSparsity;
AllSet_Sigma_AUC=Temp(SubLoc);
Temp=(sum(MAT.ElocSet(Index,AUCRange),2) - sum(MAT.ElocSet(Index,AUCRange([1 end])),2)/2)*DeltaSparsity;
AllSet_Eloc_AUC=Temp(SubLoc);
Temp=(sum(MAT.EglobSet(Index,AUCRange),2) - sum(MAT.EglobSet(Index,AUCRange([1 end])),2)/2)*DeltaSparsity;
AllSet_Eglob_AUC=Temp(SubLoc);
Temp=(sum(MAT.AssortativitySet(Index,AUCRange),2) - sum(MAT.AssortativitySet(Index,AUCRange([1 end])),2)/2)*DeltaSparsity;
AllSet_Assortativity_AUC=Temp(SubLoc);
Temp=(sum(MAT.ModularitySet(Index,AUCRange),2) - sum(MAT.ModularitySet(Index,AUCRange([1 end])),2)/2)*DeltaSparsity;
AllSet_Modularity_AUC=Temp(SubLoc);

for iNode = 1:size(MAT.DegreeSet,3)
    Temp=(sum(MAT.DegreeSet(Index,AUCRange,iNode),2) - sum(MAT.DegreeSet(Index,AUCRange([1 end]),iNode),2)/2)*DeltaSparsity;
    AllSet_Degree_AUC(:,iNode)=Temp(SubLoc);
    Temp=(sum(MAT.NodalEfficiencySet(Index,AUCRange,iNode),2) - sum(MAT.NodalEfficiencySet(Index,AUCRange([1 end]),iNode),2)/2)*DeltaSparsity;
    AllSet_NodalEfficiency_AUC(:,iNode)=Temp(SubLoc);
    Temp=(sum(MAT.BetweennessSet(Index,AUCRange,iNode),2) - sum(MAT.BetweennessSet(Index,AUCRange([1 end]),iNode),2)/2)*DeltaSparsity;
    AllSet_Betweenness_AUC(:,iNode)=Temp(SubLoc);
    Temp=(sum(MAT.ClusteringCoefficientSet(Index,AUCRange,iNode),2) - sum(MAT.ClusteringCoefficientSet(Index,AUCRange([1 end]),iNode),2)/2)*DeltaSparsity;
    AllSet_ClusteringCoefficient_AUC(:,iNode)=Temp(SubLoc);
    Temp=(sum(MAT.ParticipantCoefficientSet(Index,AUCRange,iNode),2) - sum(MAT.ParticipantCoefficientSet(Index,AUCRange([1 end]),iNode),2)/2)*DeltaSparsity;
    AllSet_ParticipantCoefficient_AUC(:,iNode)=Temp(SubLoc);
    Temp=(sum(MAT.SubgraphCentralitySet(Index,AUCRange,iNode),2) - sum(MAT.SubgraphCentralitySet(Index,AUCRange([1 end]),iNode),2)/2)*DeltaSparsity;
    AllSet_SubgraphCentrality_AUC(:,iNode)=Temp(SubLoc);
    Temp=(sum(MAT.EigenvectorCentralitySet(Index,AUCRange,iNode),2) - sum(MAT.EigenvectorCentralitySet(Index,AUCRange([1 end]),iNode),2)/2)*DeltaSparsity;
    AllSet_EigenvectorCentrality_AUC(:,iNode)=Temp(SubLoc);
    Temp=(sum(MAT.PageRankCentralitySet(Index,AUCRange,iNode),2) - sum(MAT.PageRankCentralitySet(Index,AUCRange([1 end]),iNode),2)/2)*DeltaSparsity;
    AllSet_PageRankCentrality_AUC(:,iNode)=Temp(SubLoc);
end

X=[ones(size(Dx)),Dx,Age,Sex,Edu,Motion];
Z={ones(size(Dx)),Dx};
G={Site,Site};

y=AllSet_Cp_AUC;
lme = fitlmematrix(X,y,Z,G);
T_Cp=lme.Coefficients{2,4}; %
P_Cp=lme.Coefficients{2,6};

y=AllSet_Lp_AUC;
lme = fitlmematrix(X,y,Z,G);
T_Lp=lme.Coefficients{2,4}; %
P_Lp=lme.Coefficients{2,6};

y=AllSet_Gamma_AUC;
lme = fitlmematrix(X,y,Z,G);
T_Gamma=lme.Coefficients{2,4}; %
P_Gamma=lme.Coefficients{2,6};

y=AllSet_Lambda_AUC;
lme = fitlmematrix(X,y,Z,G);
T_Lambda=lme.Coefficients{2,4}; %
P_Lambda=lme.Coefficients{2,6};

y=AllSet_Sigma_AUC;
lme = fitlmematrix(X,y,Z,G);
T_Sigma=lme.Coefficients{2,4}; %
P_Sigma=lme.Coefficients{2,6};

y=AllSet_Eloc_AUC;
lme = fitlmematrix(X,y,Z,G);
T_Eloc=lme.Coefficients{2,4}; %
P_Eloc=lme.Coefficients{2,6};

y=AllSet_Eglob_AUC;
lme = fitlmematrix(X,y,Z,G);
T_Eglob=lme.Coefficients{2,4}; %
P_Eglob=lme.Coefficients{2,6};

y=AllSet_Assortativity_AUC;
lme = fitlmematrix(X,y,Z,G);
T_Assortativity=lme.Coefficients{2,4}; %
P_Assortativity=lme.Coefficients{2,6};

y=AllSet_Modularity_AUC;
lme = fitlmematrix(X,y,Z,G);
T_Modularity=lme.Coefficients{2,4}; %
P_Modularity=lme.Coefficients{2,6};



%For Node
TMatrixDegree = zeros(size(AllSet_Degree_AUC,2),1);
PMatrixDegree = ones(size(AllSet_Degree_AUC,2),1);
TMatrixNodalEfficiency = zeros(size(AllSet_Degree_AUC,2),1);
PMatrixNodalEfficiency = ones(size(AllSet_Degree_AUC,2),1);
TMatrixBetweenness = zeros(size(AllSet_Degree_AUC,2),1);
PMatrixBetweenness = ones(size(AllSet_Degree_AUC,2),1);
TMatrixClusteringCoefficient = zeros(size(AllSet_Degree_AUC,2),1);
PMatrixClusteringCoefficient = ones(size(AllSet_Degree_AUC,2),1);
TMatrixParticipantCoefficient = zeros(size(AllSet_Degree_AUC,2),1);
PMatrixParticipantCoefficient = ones(size(AllSet_Degree_AUC,2),1);
TMatrixSubgraphCentrality = zeros(size(AllSet_Degree_AUC,2),1);
PMatrixSubgraphCentrality = ones(size(AllSet_Degree_AUC,2),1);
TMatrixEigenvectorCentrality = zeros(size(AllSet_Degree_AUC,2),1);
PMatrixEigenvectorCentrality = ones(size(AllSet_Degree_AUC,2),1);
TMatrixPageRankCentrality = zeros(size(AllSet_Degree_AUC,2),1);
PMatrixPageRankCentrality = ones(size(AllSet_Degree_AUC,2),1);
h = waitbar(0,'nodal level calculating...');
for ii=1:size(AllSet_Degree_AUC,2)
    
    y=AllSet_Degree_AUC(:,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrixDegree(ii,1)=lme.Coefficients{2,4}; %
    PMatrixDegree(ii,1)=lme.Coefficients{2,6};
    
    y=AllSet_NodalEfficiency_AUC(:,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrixNodalEfficiency(ii,1)=lme.Coefficients{2,4}; %
    PMatrixNodalEfficiency(ii,1)=lme.Coefficients{2,6};

    y=AllSet_Betweenness_AUC(:,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrixBetweenness(ii,1)=lme.Coefficients{2,4}; %
    PMatrixBetweenness(ii,1)=lme.Coefficients{2,6};   
    
    y=AllSet_ClusteringCoefficient_AUC(:,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrixClusteringCoefficient(ii,1)=lme.Coefficients{2,4}; %
    PMatrixClusteringCoefficient(ii,1)=lme.Coefficients{2,6};  
    
    y=AllSet_ClusteringCoefficient_AUC(:,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrixClusteringCoefficient(ii,1)=lme.Coefficients{2,4}; %
    PMatrixClusteringCoefficient(ii,1)=lme.Coefficients{2,6};  
    
    y=AllSet_ParticipantCoefficient_AUC(:,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrixParticipantCoefficient(ii,1)=lme.Coefficients{2,4}; %
    PMatrixParticipantCoefficient(ii,1)=lme.Coefficients{2,6};  
    
    y=AllSet_SubgraphCentrality_AUC(:,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrixSubgraphCentrality(ii,1)=lme.Coefficients{2,4}; %
    PMatrixSubgraphCentrality(ii,1)=lme.Coefficients{2,6};  
    
    y=AllSet_EigenvectorCentrality_AUC(:,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrixEigenvectorCentrality(ii,1)=lme.Coefficients{2,4}; %
    PMatrixEigenvectorCentrality(ii,1)=lme.Coefficients{2,6};  
    
    y=AllSet_PageRankCentrality_AUC(:,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrixPageRankCentrality(ii,1)=lme.Coefficients{2,4}; %
    PMatrixPageRankCentrality(ii,1)=lme.Coefficients{2,6};  
end
close(h);

%%%Plot a figure
addpath /mnt/Data/RfMRILab/Yan/YAN_Program/gretna
load /mnt/Data/RfMRILab/Yan/YAN_Program/Atlas/Dos160_WithName.mat

[pID,pN] = FDR(PMatrixDegree,0.05);
if isempty(pID); pID = 0; end
TMatrixDegree = TMatrixDegree.*(PMatrixDegree<=pID);
NameText = 'Degree';
NodeWeight=TMatrixDegree;

NodeColor=ones(size(NodeWeight));
NodeColor(NodeWeight<0)=3;%should be 2 for <0 and 3 for >0, but here is reversed for FEDN vs. recurrent
NodeColor(NodeWeight>0)=2;
ModularIndex = [1,2,3];
load('/mnt/Data/share/Software/fMRI/DPABI_V2.3_170105/Templates/Dosenbach_Science_160ROIs_Center.mat')
ROILabelUsed=[];
for iL=1:size(NodeWeight,1)
    if NodeWeight(iL)==0
        ROILabelUsed{iL,1}={{'-'}};
    else
        ROILabelUsed{iL,1}={{Dos160_WithName{iL,4}}};
    end
end

addpath /mnt/Data/RfMRILab/Yan/YAN_Program/BrainNetViewer
H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],ROILabelUsed,ModularIndex);
H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],[],ModularIndex);
set(H,'ToolBar','figure')
eval(['print -r300 -djpeg -noui ',figSavePath,'/',NameText,'.jpg;']);




[pID,pN] = FDR(PMatrixNodalEfficiency,0.05);
if isempty(pID); pID = 0; end
TMatrixNodalEfficiency = TMatrixNodalEfficiency.*(PMatrixNodalEfficiency<=pID);
NameText = 'NodalEfficiency';
NodeWeight=TMatrixNodalEfficiency;

NodeColor=ones(size(NodeWeight));
NodeColor(NodeWeight<0)=3;%should be 2 for <0 and 3 for >0, but here is reversed for FEDN vs. recurrent
NodeColor(NodeWeight>0)=2;
ModularIndex = [1,2,3];
load('/mnt/Data3/RfMRILab/ChenX/CX_software/DPABI_V4.2_190919/Templates/Dosenbach_Science_160ROIs_Center.mat')
ROILabelUsed=[];
for iL=1:size(NodeWeight,1)
    if NodeWeight(iL)==0
        ROILabelUsed{iL,1}={{'-'}};
    else
        ROILabelUsed{iL,1}={{Dos160_WithName{iL,4}}};
    end
end

addpath /mnt/Data/RfMRILab/Yan/YAN_Program/BrainNetViewer
H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],ROILabelUsed,ModularIndex);
H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],[],ModularIndex);
set(H,'ToolBar','figure')
eval(['print -r300 -djpeg -noui ',figSavePath,'/',NameText,'.jpg;']);


[pID,pN] = FDR(PMatrixBetweenness,0.05);
if isempty(pID); pID = 0; end
TMatrixBetweenness = TMatrixBetweenness.*(PMatrixBetweenness<=pID);
NameText = 'Betweenness';
NodeWeight=TMatrixBetweenness;

NodeColor=ones(size(NodeWeight));
NodeColor(NodeWeight<0)=3;%should be 2 for <0 and 3 for >0, but here is reversed for FEDN vs. recurrent
NodeColor(NodeWeight>0)=2;
ModularIndex = [1,2,3];
load('/mnt/Data/share/Software/fMRI/DPABI_V2.3_170105/Templates/Dosenbach_Science_160ROIs_Center.mat')
ROILabelUsed=[];
for iL=1:size(NodeWeight,1)
    if NodeWeight(iL)==0
        ROILabelUsed{iL,1}={{'-'}};
    else
        ROILabelUsed{iL,1}={{Dos160_WithName{iL,4}}};
    end
end

addpath /mnt/Data/RfMRILab/Yan/YAN_Program/BrainNetViewer
H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],ROILabelUsed,ModularIndex);
H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],[],ModularIndex);
set(H,'ToolBar','figure')
eval(['print -r300 -djpeg -noui ',figSavePath,'/',NameText,'.jpg;']);



[pID,pN] = FDR(PMatrixClusteringCoefficient,0.05);
if isempty(pID); pID = 0; end
TMatrixClusteringCoefficient = TMatrixClusteringCoefficient.*(PMatrixClusteringCoefficient<=pID);
% NameText = 'ClusteringCoefficient';
% NodeWeight=TMatrixClusteringCoefficient;
% 
% NodeColor=ones(size(NodeWeight));
% NodeColor(NodeWeight<0)=2;
% NodeColor(NodeWeight>0)=3;
% ModularIndex = [1,2,3];
% load('/mnt/Data3/RfMRILab/ChenX/CX_software/DPABI_V4.2_190919/Templates/Dosenbach_Science_160ROIs_Center.mat')
% ROILabelUsed=[];
% for iL=1:size(NodeWeight,1)
%     if NodeWeight(iL)==0
%         ROILabelUsed{iL,1}={{'-'}};
%     else
%         ROILabelUsed{iL,1}={{Dos160_WithName{iL,4}}};
%     end
% end
% 
% addpath /mnt/Data/RfMRILab/Yan/YAN_Program/BrainNetViewer
% % H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],ROILabelUsed,ModularIndex);
% % H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],[],ModularIndex);
% % set(H,'ToolBar','figure')
% % eval(['print -r300 -djpeg -noui ',NameText,'.jpg;']);
% 
% 
% 
[pID,pN] = FDR(PMatrixParticipantCoefficient,0.05);
if isempty(pID); pID = 0; end
TMatrixParticipantCoefficient = TMatrixParticipantCoefficient.*(PMatrixParticipantCoefficient<=pID);
% NameText = 'ParticipantCoefficient';
% NodeWeight=TMatrixParticipantCoefficient;
% 
% NodeColor=ones(size(NodeWeight));
% NodeColor(NodeWeight<0)=2;
% NodeColor(NodeWeight>0)=3;
% ModularIndex = [1,2,3];
% load('/mnt/Data3/RfMRILab/ChenX/CX_software/DPABI_V4.2_190919/Templates/Dosenbach_Science_160ROIs_Center.mat')
% ROILabelUsed=[];
% for iL=1:size(NodeWeight,1)
%     if NodeWeight(iL)==0
%         ROILabelUsed{iL,1}={{'-'}};
%     else
%         ROILabelUsed{iL,1}={{Dos160_WithName{iL,4}}};
%     end
% end
% 
% addpath /mnt/Data/RfMRILab/Yan/YAN_Program/BrainNetViewer
% %H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],ROILabelUsed,ModularIndex);
% % H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],[],ModularIndex);
% % set(H,'ToolBar','figure')
% % eval(['print -r300 -djpeg -noui ',NameText,'.jpg;']);
% % 
% 
% 
[pID,pN] = FDR(PMatrixSubgraphCentrality,0.05);
if isempty(pID); pID = 0; end
TMatrixSubgraphCentrality = TMatrixSubgraphCentrality.*(PMatrixSubgraphCentrality<=pID);
% NameText = 'SubgraphCentrality';
% NodeWeight=TMatrixSubgraphCentrality;
% 
% NodeColor=ones(size(NodeWeight));
% NodeColor(NodeWeight<0)=2;
% NodeColor(NodeWeight>0)=3;
% ModularIndex = [1,2,3];
% load('/mnt/Data3/RfMRILab/ChenX/CX_software/DPABI_V4.2_190919/Templates/Dosenbach_Science_160ROIs_Center.mat')
% ROILabelUsed=[];
% for iL=1:size(NodeWeight,1)
%     if NodeWeight(iL)==0
%         ROILabelUsed{iL,1}={{'-'}};
%     else
%         ROILabelUsed{iL,1}={{Dos160_WithName{iL,4}}};
%     end
% end
% 
% addpath /mnt/Data/RfMRILab/Yan/YAN_Program/BrainNetViewer
% %H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],ROILabelUsed,ModularIndex);
% % H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],[],ModularIndex);
% % set(H,'ToolBar','figure')
% % eval(['print -r300 -djpeg -noui ',NameText,'.jpg;']);
% 
% 
% 
[pID,pN] = FDR(PMatrixEigenvectorCentrality,0.05);
if isempty(pID); pID = 0; end
TMatrixEigenvectorCentrality = TMatrixEigenvectorCentrality.*(PMatrixEigenvectorCentrality<=pID);
% NameText = 'EigenvectorCentrality';
% NodeWeight=TMatrixEigenvectorCentrality;
% 
% NodeColor=ones(size(NodeWeight));
% NodeColor(NodeWeight<0)=2;
% NodeColor(NodeWeight>0)=3;
% ModularIndex = [1,2,3];
% load('/mnt/Data3/RfMRILab/ChenX/CX_software/DPABI_V4.2_190919/Templates/Dosenbach_Science_160ROIs_Center.mat')
% ROILabelUsed=[];
% for iL=1:size(NodeWeight,1)
%     if NodeWeight(iL)==0
%         ROILabelUsed{iL,1}={{'-'}};
%     else
%         ROILabelUsed{iL,1}={{Dos160_WithName{iL,4}}};
%     end
% end
% 
% addpath /mnt/Data/RfMRILab/Yan/YAN_Program/BrainNetViewer
% %H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],ROILabelUsed,ModularIndex);
% % H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],[],ModularIndex);
% % set(H,'ToolBar','figure')
% % eval(['print -r300 -djpeg -noui ',NameText,'.jpg;']);
% 
% 
% 
[pID,pN] = FDR(PMatrixPageRankCentrality,0.05);
if isempty(pID); pID = 0; end
TMatrixPageRankCentrality = TMatrixPageRankCentrality.*(PMatrixPageRankCentrality<=pID);
% NameText = 'PageRankCentrality';
% NodeWeight=TMatrixPageRankCentrality;
% 
% NodeColor=ones(size(NodeWeight));
% NodeColor(NodeWeight<0)=2;
% NodeColor(NodeWeight>0)=3;
% ModularIndex = [1,2,3];
% load('/mnt/Data3/RfMRILab/ChenX/CX_software/DPABI_V4.2_190919/Templates/Dosenbach_Science_160ROIs_Center.mat')
% ROILabelUsed=[];
% for iL=1:size(NodeWeight,1)
%     if NodeWeight(iL)==0
%         ROILabelUsed{iL,1}={{'-'}};
%     else
%         ROILabelUsed{iL,1}={{Dos160_WithName{iL,4}}};
%     end
% end
% 
% addpath /mnt/Data/RfMRILab/Yan/YAN_Program/BrainNetViewer
% %H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],ROILabelUsed,ModularIndex);
% % H = y_CallBrainNetViewer_NodeEdge(Dosenbach_Science_160ROIs_Center,zeros(160,160),0.1,abs(NodeWeight),1,NodeColor,[],[],ModularIndex);
% % set(H,'ToolBar','figure')
% % eval(['print -r300 -djpeg -noui ',NameText,'.jpg;']);



% Output Matrix
SiteIndex = unique(Site);
SubNumSite=zeros(length(SiteIndex), 3);
for i=1:length(SiteIndex)
    DxTemp=Dx(Site==SiteIndex(i));
    SubNumSite(i,:)=[SiteIndex(i),length(find(DxTemp==1)),length(find(DxTemp==-1))];
end

load /mnt/Data/RfMRILab/Yan/YAN_Program/Atlas/Dos160_WithName.mat
MeasureName = {'Assortativity' 'Cp' 'Eglob' 'Eloc' 'Gamma' 'Lambda' 'Lp' 'Modularity' 'Sigma'};
AucOutMat = cell(length(MeasureName), 3);
for mn = 1:length(MeasureName)
    AucOutMat{mn, 1} = MeasureName{mn};
    eval(['tmpT = T_', MeasureName{mn}, ';']);
    AucOutMat{mn, 2} = tmpT;
    eval(['tmpP = P_', MeasureName{mn}, ';']);
    AucOutMat{mn, 3} = tmpP;
end

MeasureNameNode = {'Degree', 'NodalEfficiency', 'Betweenness', 'ClusteringCoefficient', ...
    'ParticipantCoefficient', 'SubgraphCentrality', 'EigenvectorCentrality', 'PageRankCentrality'};
AucOutMatNode = cell(length(MeasureNameNode).*size(AllSet_Degree_AUC,2), 4);
IdxSig=0;
AucOutMatNodeSig = [];
for mn = 1:length(MeasureNameNode)
    for dn = 1:size(AllSet_Degree_AUC,2)
        Idx = (dn-1)*length(MeasureNameNode) + mn;
        AucOutMatNode{Idx, 1} = sprintf('CC%03d', dn);
        AucOutMatNode{Idx, 2} = MeasureNameNode{mn};
        eval(['tmpT = TMatrix', MeasureNameNode{mn}, '(dn);']);
        AucOutMatNode{Idx, 3} = tmpT;
        eval(['tmpP = PMatrix', MeasureNameNode{mn}, '(dn);']);
        AucOutMatNode{Idx, 4} = tmpP;
        if tmpT~=0
            IdxSig = IdxSig+1;
            AucOutMatNodeSig{IdxSig, 1} = MeasureNameNode{mn};
            AucOutMatNodeSig{IdxSig, 2} = sprintf('CC%03d', dn);
            AucOutMatNodeSig{IdxSig, 3} = tmpT;
            AucOutMatNodeSig{IdxSig, 4} = tmpP;
        end
    end
end
save(SavePath, 'AucOutMat', 'AucOutMatNode', 'AucOutMatNodeSig', 'SubNumSite')