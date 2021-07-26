%1. For Figure 1
clear;clc;
%Main analysi, cc200 and scrubbing
%scrubbing
metricspath = '/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SmallWorld/SmallWorld/Dos160SmallWorld/AllSM_Scrubbing/GTA_WieghtedUndirected.mat';
statspath = '/mnt/Data/RfMRILab/ChenX/Small_World/Weighted_BugFree_scrubbing';
figspath = '/mnt/Data/RfMRILab/ChenX/Small_World/Figures_scrubbing';
load /mnt/Data/RfMRILab/ChenX/Small_World/SubID_1625.mat
SubIDAll = SubID_1625;

%Get the IDs.
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

% %%% MDD vs NC
% load([statspath,'/MDD_NC/AucStat_MDD_NC.mat']);
% savepath = [figspath,'/MDD_NC'];
% if ~exist(savepath); mkdir(savepath); end
% % reverse the pallete to follow the tradition 2020/05/22
% g1color = [231/256 20/256 20/256];
% g2color = [95/256 221/256 229/256];
% g1label = 'MDD';
% g2label = 'NC';

% %%% FEDN vs NC
% WantedSubMatrix(find( (Dx==1) .* ((FirstEpisodeScore==1).*((DrugUseScore==-1))==0) ))=0;
% load([statspath,'/FEDN_NC/AucStat_FEDN_NC.mat']); 
% savepath = [figspath,'/FEDN_NC'];
% if ~exist(savepath); mkdir(savepath); end
% g1color = [249/256 216/256 156/256];
% g2color = [95/256 221/256 229/256];
% g1label = 'FEDN';
% g2label = 'NC';

% %%% Recurrent vs NC
% WantedSubMatrix(find( (Dx==1) .* ((FirstEpisodeScore==-1)==0) ))=0;
% load([statspath,'/Recurrent_NC/AucStat_Recurrent_NC.mat']);
% savepath = [figspath,'/Recurrent_NC'];
% if ~exist(savepath); mkdir(savepath); end
% g1color = [243/256 113/256 33/256];
% g2color = [95/256 221/256 229/256];
% g1label = 'Recurrent';
% g2label = 'NC';

%%% FEDN vs Recurrent
WantedSubMatrix(find( ((FirstEpisodeScore==1).*((DrugUseScore==-1))+(FirstEpisodeScore==-1)) ==0 )) = 0;
Dx = FirstEpisodeScore;
load([statspath,'/FEDN_Recurrent/AucStat_FEDN_Recurrent.mat']); 
savepath = [figspath,'/FEDN_Recurrent'];
if ~exist(savepath); mkdir(savepath); end
g1color = [249/256 216/256 156/256];
g2color = [243/256 113/256 33/256];
g1label = 'FEDN';
g2label = 'Recurrent';

%Exclude Site 4
WantedSubMatrix(find(Site==4))=0;

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

SubIDAllPlot=[];

SubIDAllPlot{1}=[];SubIDAllPlot{2}=[];
for i=1:length(SubID)
    if Dx(i) == 1
        SubIDAllPlot{1}=[SubIDAllPlot{1};SubID(i)];
    elseif Dx(i) == -1
        SubIDAllPlot{2}=[SubIDAllPlot{2};SubID(i)];
    end
end

%Get the index
IndexAllPlot{1}=[];IndexAllPlot{2}=[];
%Get the index
for iGroup = 1:length(SubIDAllPlot)
    for i=1:length(SubIDAllPlot{iGroup})
        for j=1:length(SubID)
            if strcmp(SubIDAllPlot{iGroup}{i},SubID{j})
                IndexAllPlot{iGroup}=[IndexAllPlot{iGroup};j];
            end
        end
    end
end


SubID_G1=SubIDAllPlot{1};
SubID_G2=SubIDAllPlot{2};
Index_G1=IndexAllPlot{1};
Index_G2=IndexAllPlot{2};

Age_G1=Age(Index_G1);
Edu_G1=Edu(Index_G1);
Motion_G1=Motion(Index_G1);
Age_G2=Age(Index_G2);
Edu_G2=Edu(Index_G2);
Motion_G2=Motion(Index_G2);

% Diag_BN=cell2mat(BN_Pheno(:,4:15));
% [h,p,ci,stats] = ttest2(Motion_G1,Motion_G2)
% [h,p,ci,stats] = ttest2(Age_G1,Age_G2)
% [h,p,ci,stats] = ttest2(Edu_G1,Edu_G2)

Cov_G1=[Motion_G1,Age_G1,Edu_G1];
Cov_G2=[Motion_G2,Age_G2,Edu_G2];

%Get the mean and sd smallwordness. Call SmallWorldAnalysis code
MAT = load(metricspath);
MAT.SparsityRange=[0.01:0.01:0.5]';
DeltaSparsity=0.01;
AUCRange = [10:34];
Index=[1:size(MAT.CpSet,1)];
NodeNumber = size(MAT.DegreeSet,3);

[~, SubLoc] = ismember(SubID, SubIDAll);
SubLoc(find(SubLoc == 0)) = [];

AllSet_Degree_AUC = zeros(length(SubID),NodeNumber);
AllSet_NodalEfficiency_AUC = zeros(length(SubID),NodeNumber);
AllSet_Betweenness_AUC = zeros(length(SubID),NodeNumber);
AllSet_ClusteringCoefficient_AUC = zeros(length(SubID),NodeNumber);
AllSet_ParticipantCoefficient_AUC = zeros(length(SubID),NodeNumber);
AllSet_SubgraphCentrality_AUC = zeros(length(SubID),NodeNumber);
AllSet_EigenvectorCentrality_AUC = zeros(length(SubID),NodeNumber);
AllSet_PageRankCentrality_AUC = zeros(length(SubID),NodeNumber);

%%%AUC
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

%Plot Eglob and Eloc
figure('DefaultTextFontName','Helvetica','DefaultAxesFontName','Helvetica');
set(gcf,'position',[150, 150, 1000, 500]);

subplot(1,2,1)
MeasureSet=AllSet_Eglob_AUC;
Variable1=MeasureSet(Index_G1);
Variable2=MeasureSet(Index_G2);
CovariateVariable1=Cov_G1;
CovariateVariable2=Cov_G2;
DependentVariable=[Variable1(:);Variable2(:)];

% for violin plot's matrix
outputcsv = [];
outputcsv = DependentVariable;
outputcsv(1:length(Variable1),2) = 1;
outputcsv(length(Variable1)+1:end,2) = 2;
dlmwrite([savepath,'/ForViolinPlot_Eglob.csv'],outputcsv,'delimiter',',');

% CovariateVariable=[CovariateVariable1;CovariateVariable2];
% CovariateVariable=CovariateVariable-repmat(mean(CovariateVariable),size(CovariateVariable,1),1);
% [b,r,SSE_H] = y_regress_ss(DependentVariable,[ones(length(DependentVariable),1),CovariateVariable]);
Fitted=DependentVariable;
Mean_G1=mean(Fitted(1:length(Index_G1),:));
STD_G1=std(Fitted(1:length(Index_G1),:));
Mean_G2=mean(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
STD_G2=std(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
uprim = max(1.2*(Mean_G1+STD_G1),1.2*(Mean_G2+STD_G2));
lowrim = min(0.8*(Mean_G1-STD_G1),0.8*(Mean_G2-STD_G2));
P_Eglob = AucOutMat{3,3};
X = [1:2];
Y = [Mean_G1,Mean_G2];
U = [STD_G1,STD_G2];
L = [0, 0];

bar(X(1),Y(1),'FaceColor', g1color, 'LineWidth',1.5)
hold on
bar(X(2),Y(2),'FaceColor', g2color,'LineWidth',1.5)
errorbar(X,Y,L,U,'k','LineStyle','none','LineWidth',1.5);
ylabel('E_g_l_o_b AUC','FontSize',18,'FontName','Helvetica')
set(gca,'FontSize',15);
set(gca,'FontName','Helvetica');
set(gca,'XTick',[1:2]);
set(gca,'XTickLabel',{});
xlim([0.3 2.7])
ylim([lowrim uprim])
text(1.7,lowrim+0.9*(uprim-lowrim),sprintf('p=%.3f',P_Eglob),'FontSize',18,'FontName','Helvetica')


% %Permutation Test
% Observed = mean(Fitted(1:length(Index_G1),:)) - mean(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
% N=length(Index_G1)+length(Index_G2);
% RandTimes=10000;
% RandDiffSet=[];
% for iRand=1:RandTimes
%     RandPerm=randperm(N);
%     RandDiff=mean(Fitted(RandPerm(1:length(Index_G1)),:)) - mean(Fitted(RandPerm(length(Index_G1)+1:length(Index_G1)+length(Index_G2)),:));
%     RandDiffSet(iRand,1)=RandDiff;
% end
% TwoTailedP=length(find(abs(RandDiffSet)>=abs(Observed)))/RandTimes


subplot(1,2,2)
MeasureSet=AllSet_Eloc_AUC;
Variable1=MeasureSet(Index_G1);
Variable2=MeasureSet(Index_G2);
CovariateVariable1=Cov_G1;
CovariateVariable2=Cov_G2;
DependentVariable=[Variable1(:);Variable2(:)];

% for violin plot's matrix
outputcsv = [];
outputcsv = DependentVariable;
outputcsv(1:length(Variable1),2) = 1;
outputcsv(length(Variable1)+1:end,2) = 2;
dlmwrite([savepath,'/ForViolinPlot_Eloc.csv'],outputcsv,'delimiter',',');

% CovariateVariable=[CovariateVariable1;CovariateVariable2];
% CovariateVariable=CovariateVariable-repmat(mean(CovariateVariable),size(CovariateVariable,1),1);
% [b,r,SSE_H] = y_regress_ss(DependentVariable,[ones(length(DependentVariable),1),CovariateVariable]);
Fitted=DependentVariable;
Mean_G1=mean(Fitted(1:length(Index_G1),:));
STD_G1=std(Fitted(1:length(Index_G1),:));
Mean_G2=mean(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
STD_G2=std(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
uprim = max(1.2*(Mean_G1+STD_G1),1.2*(Mean_G2+STD_G2));
lowrim = min(0.8*(Mean_G1-STD_G1),0.8*(Mean_G2-STD_G2));
P_Eloc = AucOutMat{4,3};
X = [1:2];
Y = [Mean_G1,Mean_G2];
U = [STD_G1,STD_G2];
L = [0, 0];
bar(X(1),Y(1),'FaceColor', g1color,'LineWidth',1.5)
hold on
bar(X(2),Y(2),'FaceColor', g2color,'LineWidth',1.5)
errorbar(X,Y,L,U,'k','LineStyle','none','LineWidth',1.5) ;
ylabel('E_l_o_c AUC','FontSize',18,'FontName','Helvetica')
set(gca,'FontSize',15);
set(gca,'FontName','Helvetica');
set(gca,'XTick',[1:2]);
set(gca,'XTickLabel',{});
xlim([0.3 2.7])
ylim([lowrim uprim])
text(1.7,lowrim+0.9*(uprim-lowrim),sprintf('p=%.3f',P_Eloc),'FontSize',18,'FontName','Helvetica')

% %Permutation Test
% Observed = mean(Fitted(1:length(Index_G1),:)) - mean(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
% N=length(Index_G1)+length(Index_G2);
% RandTimes=10000;
% RandDiffSet=[];
% for iRand=1:RandTimes
%     RandPerm=randperm(N);
%     RandDiff=mean(Fitted(RandPerm(1:length(Index_G1)),:)) - mean(Fitted(RandPerm(length(Index_G1)+1:length(Index_G1)+length(Index_G2)),:));
%     RandDiffSet(iRand,1)=RandDiff;
% end
% TwoTailedP=length(find(abs(RandDiffSet)>=abs(Observed)))/RandTimes

print(gcf,[savepath,'/barplot_EglobEloc.jpg'],'-djpeg','-r300');
close all;

%Plot Lp and Cp
figure('DefaultTextFontName','Helvetica','DefaultAxesFontName','Helvetica');
set(gcf,'position',[150, 150, 1000, 500]);

subplot(1,2,1)
MeasureSet=AllSet_Lp_AUC;
Variable1=MeasureSet(Index_G1);
Variable2=MeasureSet(Index_G2);
CovariateVariable1=Cov_G1;
CovariateVariable2=Cov_G2;
DependentVariable=[Variable1(:);Variable2(:)];

% for violin plot's matrix
outputcsv = [];
outputcsv = DependentVariable;
outputcsv(1:length(Variable1),2) = 1;
outputcsv(length(Variable1)+1:end,2) = 2;
dlmwrite([savepath,'/ForViolinPlot_Lp.csv'],outputcsv,'delimiter',',');

% CovariateVariable=[CovariateVariable1;CovariateVariable2];
% CovariateVariable=CovariateVariable-repmat(mean(CovariateVariable),size(CovariateVariable,1),1);
% [b,r,SSE_H] = y_regress_ss(DependentVariable,[ones(length(DependentVariable),1),CovariateVariable]);
Fitted=DependentVariable;
Mean_G1=mean(Fitted(1:length(Index_G1),:));
STD_G1=std(Fitted(1:length(Index_G1),:));
Mean_G2=mean(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
STD_G2=std(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
uprim = max(1.2*(Mean_G1+STD_G1),1.2*(Mean_G2+STD_G2));
lowrim = min(0.8*(Mean_G1-STD_G1),0.8*(Mean_G2-STD_G2));
P_Lp = AucOutMat{4,3};
X = [1:2];
Y = [Mean_G1,Mean_G2];
U = [STD_G1,STD_G2];
L = [0, 0];

bar(X(1),Y(1),'FaceColor', g1color, 'LineWidth',1.5)
hold on
bar(X(2),Y(2),'FaceColor', g2color,'LineWidth',1.5)
errorbar(X,Y,L,U,'k','LineStyle','none','LineWidth',1.5);
ylabel('L_p AUC','FontSize',18,'FontName','Helvetica')
set(gca,'FontSize',15);
set(gca,'FontName','Helvetica');
set(gca,'XTick',[1:2]);
set(gca,'XTickLabel',{});
xlim([0.3 2.7])
ylim([lowrim uprim])
text(1.7,lowrim+0.9*(uprim-lowrim),sprintf('p=%.3f',P_Lp),'FontSize',18,'FontName','Helvetica')

% %Permutation Test
% Observed = mean(Fitted(1:length(Index_G1),:)) - mean(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
% N=length(Index_G1)+length(Index_G2);
% RandTimes=10000;
% RandDiffSet=[];
% for iRand=1:RandTimes
%     RandPerm=randperm(N);
%     RandDiff=mean(Fitted(RandPerm(1:length(Index_G1)),:)) - mean(Fitted(RandPerm(length(Index_G1)+1:length(Index_G1)+length(Index_G2)),:));
%     RandDiffSet(iRand,1)=RandDiff;
% end
% TwoTailedP=length(find(abs(RandDiffSet)>=abs(Observed)))/RandTimes

subplot(1,2,2)
MeasureSet=AllSet_Cp_AUC;
Variable1=MeasureSet(Index_G1);
Variable2=MeasureSet(Index_G2);
CovariateVariable1=Cov_G1;
CovariateVariable2=Cov_G2;
DependentVariable=[Variable1(:);Variable2(:)];

% for violin plot's matrix
outputcsv = [];
outputcsv = DependentVariable;
outputcsv(1:length(Variable1),2) = 1;
outputcsv(length(Variable1)+1:end,2) = 2;
dlmwrite([savepath,'/ForViolinPlot_Cp.csv'],outputcsv,'delimiter',',');

% CovariateVariable=[CovariateVariable1;CovariateVariable2];
% CovariateVariable=CovariateVariable-repmat(mean(CovariateVariable),size(CovariateVariable,1),1);
% [b,r,SSE_H] = y_regress_ss(DependentVariable,[ones(length(DependentVariable),1),CovariateVariable]);
Fitted=DependentVariable;
Mean_G1=mean(Fitted(1:length(Index_G1),:));
STD_G1=std(Fitted(1:length(Index_G1),:));
Mean_G2=mean(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
STD_G2=std(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
uprim = max(1.2*(Mean_G1+STD_G1),1.2*(Mean_G2+STD_G2));
lowrim = min(0.8*(Mean_G1-STD_G1),0.8*(Mean_G2-STD_G2));
P_Cp = AucOutMat{2,3};
X = [1:2];
Y = [Mean_G1,Mean_G2];
U = [STD_G1,STD_G2];
L = [0, 0];

bar(X(1),Y(1),'FaceColor', g1color, 'LineWidth',1.5)
hold on
bar(X(2),Y(2),'FaceColor', g2color,'LineWidth',1.5)
errorbar(X,Y,L,U,'k','LineStyle','none','LineWidth',1.5);
ylabel('C_p AUC','FontSize',18,'FontName','Helvetica')
set(gca,'FontSize',15);
set(gca,'FontName','Helvetica');
set(gca,'XTick',[1:2]);
set(gca,'XTickLabel',{});
xlim([0.3 2.7])
ylim([lowrim uprim])
text(1.7,lowrim+0.9*(uprim-lowrim),sprintf('p=%.3f',P_Cp),'FontSize',18,'FontName','Helvetica')

% %Permutation Test
% Observed = mean(Fitted(1:length(Index_G1),:)) - mean(Fitted(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
% N=length(Index_G1)+length(Index_G2);
% RandTimes=10000;
% RandDiffSet=[];
% for iRand=1:RandTimes
%     RandPerm=randperm(N);
%     RandDiff=mean(Fitted(RandPerm(1:length(Index_G1)),:)) - mean(Fitted(RandPerm(length(Index_G1)+1:length(Index_G1)+length(Index_G2)),:));
%     RandDiffSet(iRand,1)=RandDiff;
% end
% TwoTailedP=length(find(abs(RandDiffSet)>=abs(Observed)))/RandTimes

print(gcf,[savepath,'/barplot_LpCp.jpg'],'-djpeg','-r300');
close all;


%Get difference at which sparsity
Colormap = [0.5 0.5 0.5;0 0.8 0.8;0 0.8 0;0.8 0.8 0;1 0 1;1 0 0;0 0 0];
Colormap = [0 0.8 0.8;0 0.8 0;0.8 0.8 0;1 0 1;1 0 0;0 0 0];
x = AUCRange/100;

X=[ones(size(Dx)),Dx,Age,Sex,Edu,Motion];
Z={ones(size(Dx)),Dx};
G={Site,Site};

figure;
set(gcf,'position',[150, 150, 1000, 500]);

subplot(1,2,1)
hold on
AUCRange=[10:34];
MeasureSet=[];
Contrast_EglobSet = MAT.EglobSet(SubLoc,:);
for ii=1:length(AUCRange)
    Measure = [Contrast_EglobSet(Index_G1,ii);Contrast_EglobSet(Index_G2,ii)];
%     [T,P]=y_TTest2Cov(Measure(1:length(Index_G1)),Measure(length(Index_G1)+1:length(Index_G1)+length(Index_G2)),Cov_G1,Cov_G2);
%     TMatrix(ii,1)=STATS.tstat;
%     PMatrix(ii,1)=P;
    y = MAT.EglobSet(SubLoc,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrix(ii,1)=lme.Coefficients{2,4}; %
    PMatrix(ii,1)=lme.Coefficients{2,6};
    MeasureSet(:,ii)=Measure;
end
Mean_G1=mean(MeasureSet(1:length(Index_G1),:));
STD_G1=std(MeasureSet(1:length(Index_G1),:));
Mean_G2=mean(MeasureSet(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
STD_G2=std(MeasureSet(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
MeanAll = [Mean_G1,Mean_G2];
STDAll = [STD_G1,STD_G2];
uprim = 1.2*(max(MeanAll)+max(STDAll));
lowrim = 0.8*(min(MeanAll)-max(STDAll));
errorbar(x,Mean_G2(1,:),STD_G2(1,:),'s','Color', g2color,'MarkerSize',6,'LineWidth',1.5)
errorbar(x,Mean_G1(1,:),STD_G1(1,:),'o','Color', g1color,'MarkerSize',6,'LineWidth',1.5)
legend({g2label,g1label},'Location','SouthEast');
for i=1:size(Mean_G1,2)
    if PMatrix(i)<0.05
        h(i) = plot(x(i),lowrim+0.9*(uprim-lowrim),'k*');
        set(get(get(h(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
set(gca,'FontSize',15);
set(gca,'FontName','Helvetica');
xlabel('Density','FontSize',18,'FontName','Helvetica')
xlim([0.08 0.36])
ylim([lowrim uprim])
ylabel('E_g_l_o_b','FontSize',18,'FontName','Helvetica')



subplot(1,2,2)
hold on
AUCRange=[10:34];
MeasureSet=[];
Contrast_ElocSet = MAT.ElocSet(SubLoc,:);
for ii=1:length(AUCRange)
    Measure = [Contrast_ElocSet(Index_G1,ii);Contrast_ElocSet(Index_G2,ii)];
%     [T,P]=y_TTest2Cov(Measure(1:length(Index_G1)),Measure(length(Index_G1)+1:length(Index_G1)+length(Index_G2)),Cov_G1,Cov_G2);
%     TMatrix(ii,1)=STATS.tstat;
%     PMatrix(ii,1)=P;
    y = MAT.ElocSet(SubLoc,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrix(ii,1)=lme.Coefficients{2,4}; %
    PMatrix(ii,1)=lme.Coefficients{2,6};
    MeasureSet(:,ii)=Measure;
end
Mean_G1=mean(MeasureSet(1:length(Index_G1),:));
STD_G1=std(MeasureSet(1:length(Index_G1),:));
Mean_G2=mean(MeasureSet(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
STD_G2=std(MeasureSet(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
MeanAll = [Mean_G1,Mean_G2];
STDAll = [STD_G1,STD_G2];
uprim = 1.2*(max(MeanAll)+max(STDAll));
lowrim = 0.8*(min(MeanAll)-max(STDAll));
errorbar(x,Mean_G2(1,:),STD_G2(1,:),'s','Color', g2color,'MarkerSize',6,'LineWidth',1.5)
errorbar(x,Mean_G1(1,:),STD_G1(1,:),'o','Color', g1color,'MarkerSize',6,'LineWidth',1.5)
legend({g2label,g1label},'Location','SouthEast');
for i=1:size(Mean_G1,2)
    if PMatrix(i)<0.05
        h(i) = plot(x(i),lowrim+0.9*(uprim-lowrim),'k*');
        set(get(get(h(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
set(gca,'FontSize',15);
set(gca,'FontName','Helvetica');
xlabel('Density','FontSize',18,'FontName','Helvetica')
xlim([0.08 0.36])
ylim([lowrim uprim])
ylabel('E_l_o_c','FontSize',18,'FontName','Helvetica')

print(gcf,[savepath,'/densityplot_EglobEloc.jpg'],'-djpeg','-r600');
close all;


figure;
set(gcf,'position',[150, 150, 1000, 500]);

subplot(1,2,1)
hold on
AUCRange=[10:34];
MeasureSet=[];
Contrast_LpSet = MAT.LpSet(SubLoc,:);
for ii=1:length(AUCRange)
    Measure = [Contrast_LpSet(Index_G1,ii);Contrast_LpSet(Index_G2,ii)];
%     [T,P]=y_TTest2Cov(Measure(1:length(Index_G1)),Measure(length(Index_G1)+1:length(Index_G1)+length(Index_G2)),Cov_G1,Cov_G2);
%     TMatrix(ii,1)=STATS.tstat;
%     PMatrix(ii,1)=P;
    y = MAT.LpSet(SubLoc,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrix(ii,1)=lme.Coefficients{2,4}; %
    PMatrix(ii,1)=lme.Coefficients{2,6};
    MeasureSet(:,ii)=Measure;
end
Mean_G1=mean(MeasureSet(1:length(Index_G1),:));
STD_G1=std(MeasureSet(1:length(Index_G1),:));
Mean_G2=mean(MeasureSet(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
STD_G2=std(MeasureSet(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
MeanAll = [Mean_G1,Mean_G2];
STDAll = [STD_G1,STD_G2];
uprim = 1.2*(max(MeanAll)+max(STDAll));
lowrim = 0.8*(min(MeanAll)-max(STDAll));
errorbar(x,Mean_G2(1,:),STD_G2(1,:),'s','Color',g2color,'MarkerSize',6,'LineWidth',1.5)
errorbar(x,Mean_G1(1,:),STD_G1(1,:),'o','Color',g1color,'MarkerSize',6,'LineWidth',1.5)
legend({g2label,g1label},'Location','East');
for i=1:size(Mean_G1,2)
    if PMatrix(i)<0.05
        h(i) = plot(x(i),lowrim+0.9*(uprim-lowrim),'k*');
        set(get(get(h(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
set(gca,'FontSize',15);
set(gca,'FontName','Helvetica');
xlabel('Density','FontSize',18,'FontName','Helvetica')
xlim([0.08 0.36])
ylim([lowrim uprim])
ylabel('L_p','FontSize',18,'FontName','Helvetica')

subplot(1,2,2)
hold on
AUCRange=[10:34];
MeasureSet=[];
Contrast_CpSet = MAT.CpSet(SubLoc,:);
for ii=1:length(AUCRange)
    Measure = [Contrast_CpSet(Index_G1,ii);Contrast_CpSet(Index_G2,ii)];
%     [T,P]=y_TTest2Cov(Measure(1:length(Index_G1)),Measure(length(Index_G1)+1:length(Index_G1)+length(Index_G2)),Cov_G1,Cov_G2);
%     TMatrix(ii,1)=STATS.tstat;
%     PMatrix(ii,1)=P;
    y = MAT.CpSet(SubLoc,ii);
    lme = fitlmematrix(X,y,Z,G);
    TMatrix(ii,1)=lme.Coefficients{2,4}; %
    PMatrix(ii,1)=lme.Coefficients{2,6};
    MeasureSet(:,ii)=Measure;
end
Mean_G1=mean(MeasureSet(1:length(Index_G1),:));
STD_G1=std(MeasureSet(1:length(Index_G1),:));
Mean_G2=mean(MeasureSet(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
STD_G2=std(MeasureSet(length(Index_G1)+1:length(Index_G1)+length(Index_G2),:));
MeanAll = [Mean_G1,Mean_G2];
STDAll = [STD_G1,STD_G2];
uprim = 1.2*(max(MeanAll)+max(STDAll));
lowrim = 0.8*(min(MeanAll)-max(STDAll));
errorbar(x,Mean_G2(1,:),STD_G2(1,:),'s','Color',g2color,'MarkerSize',6,'LineWidth',1.5)
errorbar(x,Mean_G1(1,:),STD_G1(1,:),'o','Color',g1color,'MarkerSize',6,'LineWidth',1.5)
legend({g2label,g1label},'Location','SouthEast');
for i=1:size(Mean_G1,2)
    if PMatrix(i)<0.05
        h(i) = plot(x(i),lowrim+0.9*(uprim-lowrim),'k*');
        set(get(get(h(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
set(gca,'FontSize',15);
set(gca,'FontName','Helvetica');
xlabel('Density','FontSize',18,'FontName','Helvetica')
xlim([0.08 0.36])
ylim([lowrim uprim])
ylabel('C_p','FontSize',18,'FontName','Helvetica')

print(gcf,[savepath,'/densityplot_LpCp.jpg'],'-djpeg','-r600');
close all;