DataDir = '/mnt/Data/RfMRILab/Lile/DFC/data_movie/FunImgRCWFS';
OutDir = '//mnt/Data/RfMRILab/Lile/DFCnew/HBNv2v/ISC/Static_MV_gmd20';
MaskFile = '/mnt/Data/RfMRILab/Lile/DFCnew/HBNv2v/GroupMask_gmd20.nii';

load /mnt/Data/RfMRILab/Lile/DFCnew/HBNv2v/SubID_RP30.mat

[mask, Header] = y_Read(MaskFile);
idx_GM = find(mask(:)>0);
cd(DataDir)
sumData = 0;
for sn = 1:numel(SubID)
    [SubData, ~] = y_ReadAll([DataDir, '/' SubID{sn}]);
    SubData = reshape(SubData, [], size(SubData,4))';
    SubData = SubData(:, idx_GM);
    sumData = sumData + SubData;
end
save('/mnt/Data/RfMRILab/Lile/DFCnew/HBNv2v/ISC/sumData_gmd20.mat', 'sumData')

if ~exist(OutDir, 'dir'); mkdir(OutDir); end
Header.pinfo = [1;0;0];
Header.dt = [16,0];
parfor sn = 1:length(SubID)
    [SubData, ~] = y_ReadAll([DataDir, '/' SubID{sn}]);
    SubData = reshape(SubData, [], size(SubData,4))';
    SubData = SubData(:, idx_GM);
    remData = sumData - SubData;

    SubISC = zeros(length(idx_GM), 1);
    for ii = 1:length(idx_GM)
        SubISC(ii) = corr(SubData(:,ii), remData(:,ii));
    end
    SubISC = 0.5 .* log( (1+SubISC) ./ (1-SubISC) );
    
    hdr = Header;
    OutImg = zeros(size(mask));
    OutImg(idx_GM) = SubISC;
    hdr.fname = [OutDir, '/ISC_', SubID{sn},  '.nii'];
    y_Write(OutImg, hdr, hdr.fname)
end
cd(fileparts(OutDir))
disp('done------done------done------done------done------done');
