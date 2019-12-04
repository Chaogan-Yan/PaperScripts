clear
tic
%Average FC as between network
MDD = load('/mnt/Data/RfMRILab/Lile/MDD/DrugUse/Dos160_ROICorrelation_MDD.mat');
HC = load('/mnt/Data/RfMRILab/Lile/MDD/DrugUse/Dos160_ROICorrelation_HC.mat');
load /mnt/Data/RfMRILab/Lile/MDD/DrugUse/Dos160_To_Yeo7_Label.mat
load /mnt/Data/RfMRILab/Yan/YAN_Program/Atlas/Dos160_WithName.mat
[SubInfoAll, SubIDAll] = xlsread('/mnt/Data/RfMRILab/Lile/MDD/DrugUse/SubInfo.xlsx', 'k2:n200');

MDDIdx = 2:2:length(MDD.SubName);       %0w: odd; 8w:even
ROICorrSet = cat(3, MDD.ROICorrSet(:,:,MDDIdx), HC.ROICorrSet);
SubName = [MDD.SubName(MDDIdx); HC.SubName];

YeoNetwork = Network142_Yeo;
[~, SortIdx] = sort(YeoNetwork);
ROICorrSet(Network==6, :, :) = [];
ROICorrSet(:, Network==6, :) = [];
Dos160_WithName(Network==6, :) = [];

%read headmotion
load /mnt/Data/RfMRILab/Lile/MDD/DrugUse/data/FDSubName.mat
[~, Loc] = ismember(SubName, SubNameAll);
if min(Loc) > 0
    FD = FDAll(Loc(Loc~=0)); 
else
    FD = FDAll(Loc(Loc~=0));
    SubName = SubNameAll(Loc(Loc~=0));
    ROICorrSet = ROICorrSet(:,:, Loc~=0);
end

%SubInfo
SubID = cell(length(SubName),1);
for sn = 1:length(SubName)
    SubID{sn,1} = SubName{sn}(1:8);
end
[~, Loc] = ismember(SubID, SubIDAll');
SubID = SubID(Loc~=0);
Dx = [ones(length(MDD.SubName)./2, 1); ones(length(HC.SubName), 1).*-1];
Dx = Dx(Loc~=0);
FD = FD(Loc~=0);
ROICorrSet = ROICorrSet(:,:, Loc~=0);

SubInfo = SubInfoAll(Loc(Loc~=0), :);
Sex = SubInfo(:,1);
Sex(Sex==2) = 0;
Age = SubInfo(:,2);
Edu = SubInfo(:,3);

%For Average FC
YeoLabel = YeoNetwork;
LabelIndex = unique(YeoLabel(YeoLabel > 0));
NetworkCorr = zeros(length(LabelIndex),length(LabelIndex));
NetworkCorrSet = zeros(length(LabelIndex),length(LabelIndex),size(ROICorrSet,3));

for iSub=1:size(ROICorrSet,3)
    CorrZ = ROICorrSet(:, :, iSub);
    for j=1:length(LabelIndex)
        for k=1:length(LabelIndex)
            A=double(YeoLabel==LabelIndex(j));
            B=double(YeoLabel==LabelIndex(k));
            Matrix = A*B';
            Matrix = tril(Matrix,-1) + triu(Matrix, 1);
            MatrixIndex = find(Matrix);
            NetworkCorr(j,k) = mean(CorrZ(MatrixIndex), 'omitnan');
        end
    end
    NetworkCorrSet(:,:,iSub)=NetworkCorr;
end

%For test
TMatrix = zeros(size(NetworkCorrSet,2),size(NetworkCorrSet,2));
PMatrix = ones(size(NetworkCorrSet,2),size(NetworkCorrSet,2));
ds = dataset(Dx, Sex, Age, Edu, FD);
for ii=1:size(NetworkCorrSet,2)
    for jj=1:size(NetworkCorrSet,2)
        y=squeeze(NetworkCorrSet(ii,jj,:));
        ds.y = y;
        lm = fitlm(ds, 'y ~ Dx + Sex + Age + FD');
        TMatrix(ii,jj) = lm.Coefficients{2,3};
        PMatrix(ii,jj) = lm.Coefficients{2,4};
        
    end
end

%FDR
TriuMat = triu(ones(size(PMatrix)),0)';
PVector = PMatrix(logical(TriuMat));
addpath /mnt/Data/RfMRILab/Yan/YAN_Program/gretna
[pID,pN] = FDR(PVector,0.05);
if ~isempty(pID)
    PSig = PMatrix<=pID; 
else
    PSig = zeros(size(PMatrix));
end
toc