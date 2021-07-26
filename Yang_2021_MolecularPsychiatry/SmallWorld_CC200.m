
addpath('/mnt/Data/RfMRILab/Yan/YAN_Program');
addpath('/mnt/Data/RfMRILab/Yan/YAN_Program/BCT_20120814');


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




DataDir='/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Raw/REST-meta-MDD-Organized/Results/ROISignals_FunImgARCWF';
%Bonds=[1 116; 117 212; 213 228; 229 428; 429 1408; 1409 1568; 1570 1833];
Bonds=[229 428];
j=1;
ROICorrelation_FisherZ_Set=zeros(Bonds(j,2)-Bonds(j,1)+1,Bonds(j,2)-Bonds(j,1)+1,length(SubID));
for i=1:length(SubID)
    ROISignals=load([DataDir,'/','ROISignals_',SubID{i},'.mat']);
    Sig=ROISignals.ROISignals(:,Bonds(j,1):Bonds(j,2));
    ROICorrelation = corrcoef(Sig);
    ROICorrelation_FisherZ = 0.5 * log((1 + ROICorrelation)./(1- ROICorrelation));
    ROICorrelation_FisherZ(find(eye(size(ROICorrelation_FisherZ))))=0;
    ROICorrelation_FisherZ_Set(:,:,i)=ROICorrelation_FisherZ;
end




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

ROICorrelation_FisherZ_Set = ROICorrelation_FisherZ_Set(:,:,WantedSubIndex);



%Perform graph theoretical analysis


ToUpDir = '/mnt/Data4/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SmallWorld/SmallWorld/CC200SmallWorld/AllSM_Weighted';
mkdir(ToUpDir)
OutputTempDir='/mnt/Data4/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SmallWorld/SmallWorld/CC200SmallWorld/Temp';
mkdir(OutputTempDir)

NodeNumber=200;

SubjectNum=length(SubID);

OutputDir = ToUpDir;
mkdir(OutputDir);

%Perform graph theoretical analysis (wieghted undirected)

%SparsityRange = [0.02:0.02:1]';
SparsityRange =[0.01:0.01:0.5]';

nSparsity = length(SparsityRange);

MaxSparsitySet = ones(SubjectNum,1);

CpSet = zeros(SubjectNum,nSparsity);
LpSet = zeros(SubjectNum,nSparsity);
GammaSet = zeros(SubjectNum,nSparsity);
LambdaSet = zeros(SubjectNum,nSparsity);
SigmaSet = zeros(SubjectNum,nSparsity);
ElocSet = zeros(SubjectNum,nSparsity);
EglobSet = zeros(SubjectNum,nSparsity);
AssortativitySet = zeros(SubjectNum,nSparsity);
ModularitySet = zeros(SubjectNum,nSparsity);

DegreeSet = zeros(SubjectNum,nSparsity,NodeNumber);
NodalEfficiencySet = zeros(SubjectNum,nSparsity,NodeNumber);
BetweennessSet = zeros(SubjectNum,nSparsity,NodeNumber);
ClusteringCoefficientSet = zeros(SubjectNum,nSparsity,NodeNumber);
ParticipantCoefficientSet = zeros(SubjectNum,nSparsity,NodeNumber);
SubgraphCentralitySet = zeros(SubjectNum,nSparsity,NodeNumber);
EigenvectorCentralitySet = zeros(SubjectNum,nSparsity,NodeNumber);
PageRankCentralitySet = zeros(SubjectNum,nSparsity,NodeNumber);

%[Cp, Lp, Gamma, Lambda, Sigma, Eloc, Eglob, Assortativity, NO!!Hierarchy, NO!!Synchronization, Modularity]
%[Degree, NodalEfficiency, Betweenness, ClusteringCoefficient, ParticipantCoefficient, NO!!NormalizedParticipantCoefficient, SubgraphCentrality, EigenvectorCentrality, PageRankCentrality]

for i=1:SubjectNum
    FisherZSet = ROICorrelation_FisherZ_Set(:,:,i);
    
    FisherZSet = FisherZSet(1:NodeNumber,1:NodeNumber);
    
    MaxSparsitySet(i) = length(find(FisherZSet>0)) / (NodeNumber*(NodeNumber-1));
    
    MatrixValueThreshold_LastRound = 1;
    
    parfor iSparsity=1:length(SparsityRange)
        
        fprintf('Sparsity: %d for Subject %s',SparsityRange(iSparsity),SubID{i});
        
        if (SparsityRange(iSparsity)<=MaxSparsitySet(i)) % If not reached the maximum sparsity
            
            Matrix = FisherZSet;
            Matrix(find(eye(size(Matrix))))=0;
            
            EdgeNumber=round(SparsityRange(iSparsity)*NodeNumber*(NodeNumber-1));
            MatrixValue=Matrix(:);
            MatrixValueSorted=sort(MatrixValue,'descend');
            MatrixValueThreshold=MatrixValueSorted(EdgeNumber+1);
            
            
            Matrix=Matrix.*(Matrix>MatrixValueThreshold);
            
            [GTA] = y_GraphTheoreticalAnalysis_wu(Matrix,100);
            
            CpSet(i,iSparsity) = GTA.Cp;
            LpSet(i,iSparsity) = GTA.Lp;
            GammaSet(i,iSparsity) = GTA.Gamma;
            LambdaSet(i,iSparsity) = GTA.Lambda;
            SigmaSet(i,iSparsity) = GTA.Sigma;
            ElocSet(i,iSparsity) = GTA.Eloc;
            EglobSet(i,iSparsity) = GTA.Eglob;
            AssortativitySet(i,iSparsity) = GTA.Assortativity;
            ModularitySet(i,iSparsity) = GTA.Modularity;
            
            DegreeSet(i,iSparsity,:) = GTA.Degree;
            NodalEfficiencySet(i,iSparsity,:) = GTA.NodalEfficiency;
            BetweennessSet(i,iSparsity,:) = GTA.Betweenness;
            ClusteringCoefficientSet(i,iSparsity,:) = GTA.ClusteringCoefficient;
            ParticipantCoefficientSet(i,iSparsity,:) = GTA.ParticipantCoefficient;
            SubgraphCentralitySet(i,iSparsity,:) = GTA.SubgraphCentrality;
            EigenvectorCentralitySet(i,iSparsity,:) = GTA.EigenvectorCentrality;
            PageRankCentralitySet(i,iSparsity,:) = GTA.PageRankCentrality;
            
        end
    end
    
    save([OutputTempDir,filesep,'GTA_WieghtedUndirected_',SubID{i},'.mat'],'CpSet','LpSet','GammaSet','LambdaSet','SigmaSet','ElocSet','EglobSet','AssortativitySet','ModularitySet','DegreeSet','NodalEfficiencySet','BetweennessSet','ClusteringCoefficientSet','ParticipantCoefficientSet','SubgraphCentralitySet','EigenvectorCentralitySet','PageRankCentralitySet','MaxSparsitySet','SparsityRange');
    
end

save([OutputDir,filesep,'GTA_WieghtedUndirected.mat'],'SubID','CpSet','LpSet','GammaSet','LambdaSet','SigmaSet','ElocSet','EglobSet','AssortativitySet','ModularitySet','DegreeSet','NodalEfficiencySet','BetweennessSet','ClusteringCoefficientSet','ParticipantCoefficientSet','SubgraphCentralitySet','EigenvectorCentralitySet','PageRankCentralitySet','MaxSparsitySet','SparsityRange');




