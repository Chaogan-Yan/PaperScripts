clear
TplImg = y_Read([fileparts(which('dpabi')), '/Templates/Dosenbach_Science_160ROIs_Radius5_Mask.nii']);
load /mnt/Data/RfMRILab/Lile/MDD/DrugUse/Dos160_To_Yeo7_Label.mat
DataDir = '/mnt/Data/RfMRILab/Lile/MDD/DrugUse/data/MDD/FunImgARCWF';
SaveFile = '/mnt/Data/RfMRILab/Lile/MDD/DrugUse/Dos160_ROICorrelation_MDD.mat';

SubPath = spm_select('FPListRec', DataDir, '^Filtered*');

TplNum = max(TplImg(:));
ROICorrSet = zeros(TplNum, TplNum, size(SubPath,1));
SubNanNum = zeros(size(SubPath,1),1);
TplIdx = cell(TplNum,1);
SubName = cell(size(SubPath,1),1);
tic
for tn = 1:TplNum
    TplIdx{tn,1} = find(TplImg==tn);
end

parfor sn=1:size(SubPath,1)
    Img = y_Read(deblank(SubPath(sn,:)));
    TP = size(Img, 4);
    Img = reshape(Img, [], TP);
    Sig = zeros(TP, TplNum);
    for tn = 1:TplNum
        Sig(:, tn) = mean(Img(TplIdx{tn}, :), 1);
    end
    ROICorrelation = corrcoef(Sig);
    NanNum = find(isnan(ROICorrelation));
    if ~isempty(NanNum)
        NanNum = TplNum-sqrt(TplNum^2 - length(NanNum));
        SubNanNum(sn) = NanNum;
    end
    SubNameTmp = fileparts(deblank(SubPath(sn,:)));
    [~, SubNameTmp] = fileparts(SubNameTmp);
    SubName{sn} = SubNameTmp;
    fprintf('%d NaN region in %s\n', SubNanNum(sn), SubName{sn})
    ROICorr = 0.5 * log((1 + ROICorrelation)./(1- ROICorrelation));
    ROICorrSet(:,:, sn)=ROICorr;
end
toc
save(SaveFile, 'ROICorrSet', 'SubNanNum', 'SubName')