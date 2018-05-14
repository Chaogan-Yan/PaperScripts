
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




DataDir='/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Raw/REST-meta-MDD-Organized/Results/ROISignals_FunImgARCWF';
%Bonds=[1 116; 117 212; 213 228; 229 428; 429 1408; 1409 1568; 1570 1833];
Bonds=[429 1408];
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




load /mnt/Data/RfMRILab/Yan/YAN_Program/Atlas/Zalesky980_YeoNetwork.mat
Temp=load('/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_848_794/ROIAnalysis/sgACC_Zalesky/FCMatrix/TMatrix_All.mat');
ROI_Remove=Temp.NaNROI(1:end-1); %We still need the 120th ROI
ROIWanted=1:980;
ROIWanted(ROI_Remove)=[];
YeoNetwork(ROI_Remove)=[];
ROIWanted(find(YeoNetwork==0))=[];
YeoNetwork(find(YeoNetwork==0))=[];

ROICorrelation_FisherZ_Set=ROICorrelation_FisherZ_Set(ROIWanted,ROIWanted,:);



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




