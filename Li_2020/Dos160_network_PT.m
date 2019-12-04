clear
tic
load /mnt/Data/RfMRILab/Lile/MDD/DrugUse/Dos160_ROICorrelation_MDD.mat
load /mnt/Data/RfMRILab/Lile/MDD/DrugUse/Dos160_To_Yeo7_Label.mat
load /mnt/Data/RfMRILab/Yan/YAN_Program/Atlas/Dos160_WithName.mat

YeoNetwork = Network142_Yeo;
[~, SortIdx] = sort(YeoNetwork);
ROICorrSet(Network==6, :, :) = [];
ROICorrSet(:, Network==6, :) = [];
Dos160_WithName(Network==6, :) = [];

ObsNum = length(SubName);
% 0w and then 8w
SubName = SubName([1:2:ObsNum, 2:2:ObsNum]);        
ROICorrSet = ROICorrSet(:, :, [1:2:ObsNum, 2:2:ObsNum]);

%% read headmotion
load /mnt/Data/RfMRILab/Lile/MDD/DrugUse/data/FDSubName.mat
[~, Loc] = ismember(SubName, SubNameAll);
if min(Loc) > 0
    FD = FDAll(Loc); 
else
    error('FD missing'); 
end
FD = reshape(FD, [], 2);
SubID = cell(length(SubName)/2,1);
for sn = 1:length(SubName)/2
    SubID{sn,1} = SubName{sn}(1:8);
end

%% SubInfo and select sub
[SubInfoAll, SubIDAll] = xlsread('/mnt/Data/RfMRILab/Lile/MDD/DrugUse/SubInfo.xlsx', 1, 'a2:i50');
[~, Loc] = ismember(SubID, SubIDAll');
SubID = SubID(Loc~=0);
ROICorrSet = ROICorrSet(:,:, [Loc;Loc]~=0);
FD = FD(Loc~=0, :);
SubInfo = SubInfoAll(Loc(Loc~=0), :);
HAMD = [SubInfo(:,6), SubInfo(:,7)];

%% For Average FC
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

%% For Ttest
ObsNum = size(ROICorrSet,3);
TMatrix = zeros(size(NetworkCorrSet,2),size(NetworkCorrSet,2));
PMatrix = ones(size(NetworkCorrSet,2),size(NetworkCorrSet,2));
for ii=1:size(NetworkCorrSet,2)
    for jj=1:size(NetworkCorrSet,2)
        y=squeeze(NetworkCorrSet(ii,jj,:));
        [~, p, ~, stats] = ttest(y(ObsNum/2+1:end), y(1:ObsNum/2));
        TMatrix(ii,jj) = stats.tstat;
        PMatrix(ii,jj) = p;
    end
end

%FDR
TrilMat = tril(ones(size(PMatrix)), 0);
PVector = PMatrix(logical(TrilMat));
addpath /mnt/Data/RfMRILab/Yan/YAN_Program/gretna
[pID,~] = FDR(PVector,0.05);
if ~isempty(pID)
    PSig = PMatrix<=pID; 
else
    PSig = zeros(size(PMatrix));
end
toc