%% ASD-EOS neuroanatomic and connectomic profile analysis


%% %%%%%% Organize SubInfo %%%%%%
%% Organize SubInfo - sort the raw subject phenotypical records
clc;clear;
Version = '20230912';
Indir = '/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/ASD_SCZ';
Sites = {'Beijing','JingWei','ShenZhen','XiangYa','HuaShi'};
Folders = {'ForDPABI_20210507','ForDPABI_20210514','ForDPABI_20210511','ForDPABI_20210514','ForDPABI_20210514'};
SubTypes = {'ASD','NC','SCZAcute','SCZProdromal','SCZRisk','ASDRisk'};
ID = {};
Dx = [];
Site = [];
Dir = {};
for iSite = 1:length(Sites)
    SubList0 = dir([Indir,filesep,Sites{iSite},filesep,Folders{iSite},filesep,'FunVoluW',filesep,'sub*']);
    Dx0 = zeros(length(SubList0),1);
    for iSub = 1:length(SubList0)
        Index = cellfun(@(Dx) contains(SubList0(iSub).name(5:5+length(Dx)-1),Dx),SubTypes);
        if length(find(Index))==1
            Dx0(iSub,1) = find(Index);
        elseif length(find(Index))>1 && Index(6) %ASDRisk
            Dx0(iSub,1) = 6;
        else
            Dx0(iSub,1) = nan;
            disp([num2str(iSite),filesep,num2str(iSub)])
        end
    end
    ID = [ID;{SubList0(:).name}'];
    Dir = [Dir;{SubList0(:).folder}'];
    Site = [Site;ones(length(SubList0),1).*iSite];
    Dx = [Dx;Dx0];
end
Dir = cellfun(@(D) D(1:end-8),Dir,'UniformOutput',false); % remove 'FunVoluW'

% Recover the original ID
ID_Raw = {};
for iSub = 1:length(ID)
    List = tdfread([Dir{iSub},filesep,'SubjectID_DPARSFA2BIDS.tsv']);
    Index = cellfun(@(ID_BIDS) strcmp(ID_BIDS,ID{iSub}),cellstr(List.SubjectID_BIDS));
    ID_Raw{iSub,1} = regexprep(List.SubjectID_Original(find(Index),:),' ','');
end

%%% Find site, age and sex
SiteName = {};
Sex = ones(length(ID),1)*99;
Age = zeros(length(ID),1);
Weight = zeros(length(ID),1);
ScanDate = {};
BirthDate = {};

for iSub = 1:length(ID_Raw)
    temp = dir([Dir{iSub},filesep,'T1Raw',filesep,ID_Raw{iSub}]);
    Info = dicominfo([temp(3).folder,filesep,temp(3).name]);
    ScanDate{iSub,1} = Info.StudyDate;
    if strcmp(Info.PatientSex,'F')
        Sex(iSub,1) = 0;
    elseif strcmp(Info.PatientSex,'M')
        Sex(iSub,1) = 1;
    end
    Age(iSub,1) = str2num(Info.PatientAge(1:3));
    Weight(iSub,1) = Info.PatientWeight;
    BirthDate{iSub,1} = Info.PatientBirthDate;
    SiteName{iSub,1} = Info.InstitutionName;
    disp(num2str(iSub))
end

SiteFromData = zeros(length(ID_Raw),1);
for i = 1:length(ID_Raw)
    switch SiteName{i}
        case 'PKU 6th HP'
            SiteFromData(i) = 1;
        case {' The 2nd XiangYa Hospital of CSU Hospital',' The 2nd XiangYa Hospital of CSU/F82F55/Hospital'}
            SiteFromData(i) = 4;    
        case 'SH Mental Health Center'
            SiteFromData(i) = 2;   
        case 'SZKANGNING'
            SiteFromData(i) = 3;
        case {'EAST CHINA NORMAL UNIVERSITY','East China Normal University'}
            SiteFromData(i) = 5;
    end
end
Site = SiteFromData; %% Some raw site indices from Shanghai are wrong, therefore use the records in the DICOM header



%% Organize SubInfo - deal with human error in the raw records
% Wrong Dx
DxError_ID = {'ASDchenmin20180508','ASDlikunming20180801','ASDyangsiqi20180807','SCZAcutewangsijia20180705','ASDzhangyuran20200111','ASDzhengwenwen20200111'};
DxCorrect = [3,4,5,4,6,6];
for iSub = 1:length(DxError_ID)
    Index = find(contains(ID,DxError_ID{iSub}));disp(Dx(Index))
    Dx(Index) = DxCorrect(iSub);
end

% Wrong age
AgeError = {'ASDliuyiran20180305','ASDSiYuanWang20180816'};
AgeCorrect = [10,4.7];
for iSub = 1:length(AgeError)
    Index = find(contains(ID,AgeError{iSub}));disp(Age(Index))
    Age(Index) = AgeCorrect(iSub);
end

% Wrong Sex
SexError = {'NCWangziyan20170119'};
SexCorrect = [0];
for iSub = 1:length(SexError)
    Index = find(contains(ID,SexError{iSub}));disp(Sex(Index))
    Sex(Index) = SexCorrect(iSub);
end

% Repetitive subjects (e.g. typo in name)
RepeatError = {'ASDZhangyanxiang20170118','ASDmaminbo20190722'};
for iSub = 1:length(RepeatError)
    Index = find(contains(ID,RepeatError{iSub}));disp(Index)
    ID(Index) = [];
    ID_Raw(Index) = [];
    Sex(Index) = [];
    BirthDate(Index) = [];
    ScanDate(Index) = [];
    Weight(Index) = [];
    Age(Index) = [];
    SiteName(Index) = [];
    Dir(Index) = [];
    Site(Index) = [];
    Dx(Index) = [];
end

[C,IA,IC] = unique(ID);
SubInfo.ID = ID(IA);
SubInfo.ID_Raw = ID_Raw(IA);
SubInfo.Dir = Dir(IA);
SubInfo.Site = Site(IA);
SubInfo.Sex = Sex(IA);
SubInfo.Dx = Dx(IA);
SubInfo.SiteName = SiteName(IA);
SubInfo.Age = Age(IA);
SubInfo.Weight = Weight(IA);
SubInfo.ScanDate = ScanDate(IA);
SubInfo.BirthDate = BirthDate(IA);
save(['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_',Version,'.mat'],'SubInfo');



%% Organize SubInfo - add functional MRI info: 1) pseudo series flags 2) head motion
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat');
RestingFlags = ones(length(SubInfo.ID),1);
MotionFD = zeros(length(SubInfo.ID),1);
Motion3mmDegFlags = ones(length(SubInfo.ID),1);
for iSub = 1:length(SubInfo.ID)
    temp = dir([SubInfo.Dir{iSub},filesep,'*DPABI_Format*.csv']);
    ConvertInfo = readtable([temp.folder,filesep,temp.name]);
    for j = 1:size(ConvertInfo,1)
        if strcmp(SubInfo.ID_Raw{iSub},ConvertInfo{j,1}) && strcmp('PseudoSeries',ConvertInfo{j,3})
            RestingFlags(iSub) = 0;
        end
    end
    MotionInfo = readtable([SubInfo.Dir{iSub},filesep,'RealignParameter',filesep,'HeadMotion.tsv'],'FileType','text');
    for j = 1:size(MotionInfo,1)
        if strcmp(SubInfo.ID{iSub},MotionInfo{j,1})
            MotionFD(iSub,1) = MotionInfo{j,21};
            if ~isempty(find(MotionInfo{j,2:7}>3))
                Motion3mmDegFlags(iSub,1) = 0;
            end
        end
    end
    disp(num2str(iSub))
end
SubInfo.MotionFD = MotionFD;
SubInfo.Motion3mmDegFlags = Motion3mmDegFlags;
SubInfo.RestFlag = RestingFlags;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat','SubInfo');



%% Organize SubInfo - add fmriprep report results for QC
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/QCtempFinal.mat');
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat');
QCFlags = ones(length(SubInfo.ID),1);

[C,IA,IB] = intersect(SubInfo.ID,BadList);
QCFlags(IA) = 0;

SubInfo.QCFlag = QCFlags;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat','SubInfo');


%% Organize SubInfo - Create sublist for functional analysi
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat')
Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);
Fields = fieldnames(SubInfo);
for i = 1:length(Fields)
    SubInfo.(Fields{i}) = SubInfo.(Fields{i})(Index);
end
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','SubInfo');



%% %%%%%% ANCOVA %%%%%%
%% ANCOVA - Subcortical volume - ComBat and statistic
clc;clear;
addpath(genpath('/mnt/Data6/RfMRILab/Lubin/Software/ComBatHarmonization-master'));
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat')
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubcorticalData_Version1_20210521.mat');
ValueAll = [Value_Beijing;Value_XiangYa;Value_Jingwei;Value_Shenzhen;Value_Huashi];
IDAll = [ID_Beijing;ID_XiangYa;ID_Jingwei;ID_Shenzhen;ID_Huashi];
for iSub = 1:length(SubInfo.ID)
    for j = 1:length(IDAll)
        if strcmp(SubInfo.ID{iSub},IDAll{j})
            AllVolume(iSub,:) = ValueAll(j,:);
        end
    end
end

SubInfo.ICV = AllVolume(:,end);

AreaIndex = [1,5,6,7,8,12,13,15,19,23,24,25,26,27,28,29]; % discard 11-brainstem,16/30-ventralDC,40 Optic-Chiasm

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 );
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
Sex = SubInfo.Sex(Index);
Site = SubInfo.Site(Index);
Volume = AllVolume(Index,AreaIndex);
ICV = SubInfo.ICV(Index);

AllCov_comBat = [Dx,Age,Sex];
Volume_Combat = combat(Volume', Site', AllCov_comBat, 1);
Volume_Combat = Volume_Combat';

DxCov=[];
DxIndex = [1,2,3];
for i=1:length(DxIndex)-1
    DxCov(:,i) = Dx==DxIndex(i);
end

AllCov = [DxCov,ones(length(Dx),1),Age,Sex,ICV];

Contrast=zeros(1,size(AllCov,2));
Contrast(1:size(DxCov,2))=1;

for i = 1:size(Volume_Combat,2)
    [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(Volume_Combat(:,i),AllCov,Contrast,'F');
    Result_F(i,1) = TF_ForContrast;
    Result_P(i,1) = 1-fcdf(TF_ForContrast,size(DxCov,2),size(DxCov,	1)-size(AllCov,2));
end

for i  = 1:length(AreaIndex)
    disp([num2str(i),filesep,num2str(AreaIndex(i)),filesep,Title{AreaIndex(i)}])
end

Result_P_Corrected = mafdr(Result_P,'BHFDR','true');
SigFlags = find(Result_P_Corrected<0.05); 
Result_P_Corrected(find(Result_P_Corrected<0.05))
SigArea = Title(AreaIndex(SigFlags))
%  "Left-Lateral-Vent…"  "Left-Thalamus-Proper"  "Left-Hippocampus"    "Left-Amygdala"    "Right-Lateral-Ven…"    "Right-Hippocampus"    "Right-Amygdala"

% Save ROI values
Volume_Combat_All = zeros(length(SubInfo.ID),length(SigFlags));
Volume_Combat_All(Index,:) = Volume_Combat(:,SigFlags);

ROIValue.SubcortVol_F = Result_F(SigFlags);
ROIValue.SubcortVol_P = Result_P_Corrected(SigFlags);
ROIValue.SubcortVol_Sig = Volume_Combat_All;
ROIValue.SubcortVol_SigFlag = SigArea;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat','SubInfo','ROIValue')



%% ANCOVA - Subcortical volume - (resample subcortical ROI to fs5)
TemplateList = {'LatVentricle_L_mask.nii','Thalamus_L_mask.nii','Hippocampus_L_mask.nii','Amygdala_L_mask.nii','LatVentricle_R_mask.nii','Hippocampus_R_mask.nii','Amygdala_R_mask.nii'};
DataAll = [];
mkdir('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/Anat');

for iTemp = 1:length(TemplateList)
    FValue = Result_F(SigFlags(iTemp));
    [Data, VoxelSize, FileList, Header] = y_ReadAll(['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Templates',filesep,TemplateList{iTemp}]);
    Data(find(Data)) = FValue;
    Header.dt = 16;
    if iTemp == 1
        DataAll = Data;
    else
        DataAll = DataAll+Data;
    end
    name =  regexprep(TemplateList{iTemp},'_mask','');
    y_Write(Data,Header,['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/Anat',...
        filesep,'Volume_',name]);
end
y_Write(DataAll,Header,['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/Anat',...
    filesep,'Volume_SubcorticalAll.nii']);





%% ANCOVA - Cortical structural metrics - ComBat 
clc;clear;
Date = '20230912_ICVCovariate';
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat')
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults',filesep,Date,...
    filesep,'Combat_Metrics',filesep,'Cortical_Structure'];
mkdir(OutDir);
Metrics  = {'Thickness','Area','Curv','Sulc','Volume'};
Hemispheres = {'AnatSurfLH','AnatSurfRH'};
Masks = {'/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V5.1_201201/DPABISurf/SurfTemplates/fsaverage_lh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V5.1_201201/DPABISurf/SurfTemplates/fsaverage_rh_cortex.label.gii'};

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 );
ID = SubInfo.ID(Index);
Dir = SubInfo.Dir(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
Sex = SubInfo.Sex(Index);
Site = SubInfo.Site(Index);

AllCov_comBat = [Dx,Age,Sex];

for iHemi = 1:length(Hemispheres)
    MaskData = y_ReadAll(Masks{iHemi});
    MaskIndex = find(MaskData);
    for iMetric = 1:length(Metrics)
        OutputDir = [OutDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric}];
        mkdir(OutputDir);
        ImageList = cellfun(@(Dir,ID) dir([Dir,filesep,'ResultsS',filesep,Hemispheres{iHemi},filesep,...
            Metrics{iMetric},filesep,'fsaverage',filesep,'s',ID,'*']),...
            Dir,ID);
        ImageList = cellfun(@(Dir,ID) [Dir,filesep,ID],{ImageList(:).folder},{ImageList(:).name},'UniformOutput',false);
        [Data, VoxelSize, FileList, Header] = y_ReadAll(ImageList);
        Data = Data(MaskIndex,:);
        Volume_Combat = combat(Data, Site', AllCov_comBat, 1);
        for iSub = 1:size(Volume_Combat,2)
            Brain = zeros(size(MaskData));
            Brain(MaskIndex) = Volume_Combat(:,iSub);
            y_Write(Brain,Header,[OutputDir,filesep,ID{iSub}]);
            disp([Hemispheres{iHemi},filesep,Metrics{iMetric},filesep,'Sub',num2str(iSub)])
        end
    end
end
        
       
%% ANCOVA - Cortical structural metrics - Statistic 
clc;clear;
Date = '20230912_ICVCovariate';
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat')
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/Cortical_Structure';
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults',filesep,Date,filesep,'Cortical_Metrics'];
mkdir(OutDir);
Metrics  = {'Thickness','Area','Curv','Sulc','Volume'};
Hemispheres = {'AnatSurfLH','AnatSurfRH'};
Masks = {'/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V5.1_201201/DPABISurf/SurfTemplates/fsaverage_lh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V5.1_201201/DPABISurf/SurfTemplates/fsaverage_rh_cortex.label.gii'};

for iHemi = 1:length(Hemispheres)
    parfor iMetric = 1:length(Metrics)
        OutputDir = [OutDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric},filesep];
        mkdir(OutputDir);
        
        ASDIndex = find(SubInfo.Dx==1 & SubInfo.QCFlag>0 & SubInfo.Age<=18);
        NCIndex = find(SubInfo.Dx==2 & SubInfo.QCFlag>0 & SubInfo.Age<=18);
        SCZIndex = find(SubInfo.Dx==3 & SubInfo.QCFlag>0 & SubInfo.Age<=18);
        
        ASDList = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric},...
            filesep,ID,'.gii'],SubInfo.ID(ASDIndex),'UniformOutput',false);
        NCList = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric},...
            filesep,ID,'.gii'],SubInfo.ID(NCIndex),'UniformOutput',false);
        SCZList = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric},...
            filesep,ID,'.gii'],SubInfo.ID(SCZIndex),'UniformOutput',false);
  
        switch iMetric % Kiho IM, et al, Cerebral Cortex, 2008; {'Area','Curv','Sulc','Thickness','Volume'};
            case 1 % 'Area'
                ScalingFactor = 2/3;
            case 2 % 'Curv'
                ScalingFactor = -1/3;
            case 3 % 'Sulc'
                ScalingFactor = 1/3;
            case 4 % 'Thickness'
                ScalingFactor = 1/3;
            case 5 % 'Volume'
                ScalingFactor = 1;
        end                
                
        CovariateMatrix = {[SubInfo.Age(ASDIndex),SubInfo.Sex(ASDIndex),SubInfo.ICV(ASDIndex).^(ScalingFactor)];...
            [SubInfo.Age(NCIndex),SubInfo.Sex(NCIndex),SubInfo.ICV(NCIndex).^(ScalingFactor)];...
            [SubInfo.Age(SCZIndex),SubInfo.Sex(SCZIndex),SubInfo.ICV(SCZIndex).^(ScalingFactor)]};

        y_ANCOVA1_Image({ASDList,NCList,SCZList},[OutputDir,filesep,'Ancova'],Masks{iHemi},{},CovariateMatrix);
        disp(['Have down ',Hemispheres{iHemi},' ',Metrics{iMetric}])
    end
end


%% ANCOVA - Cortical structural metrics - (resample cortical ROI to fs5)
% Deal with header problem
clc;clear;
InputDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/Anat';
InName = {'Thickness_L_1_Mask.gii',...
    'Thickness_R_1_Mask.gii',...
    'Thickness_R_2_Mask.gii'};
[Data,~,~,Header] = y_ReadAll('/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/ASD_SCZ/XiangYa/ForDPABI_20210514/ResultsS/AnatSurfLH/Thickness/fsaverage/ssub-SCZRiskdengyihan20200613_space-fsaverage_hemi-L.thickness.gii');
for iSub = 1:length(InName)
    Data = y_ReadAll([InputDir,filesep,InName{iSub}]);
    y_Write(Data,Header,[InputDir,filesep,InName{iSub}(1:end-4),'_NewHeader']);
end

[DPABIPath, fileN, extn] = fileparts(which('DPABI.m'));
CommandInit=sprintf('docker run -ti --rm -v %s:/opt/freesurfer/license.txt -v %s:/data ', fullfile(DPABIPath, 'DPABISurf', 'FreeSurferLicense', 'license.txt'), InputDir); %YAN Chao-Gan, 181214. Remove -t because there is a tty issue in windows

InName = {'Thickness_L_1_Mask_NewHeader.gii',...
    'Thickness_R_1_Mask_NewHeader.gii',...
    'Thickness_R_2_Mask_NewHeader.gii'};
OutName = {'Thickness_L_1_Mask_fs5.gii',...
    'Thickness_R_1_Mask_fs5.gii',...
    'Thickness_R_2_Mask_fs5.gii'};
Hemis = {'lh','rh','rh'};

for iROI = 1:3
    Command=sprintf('%s cgyan/dpabi mri_surf2surf --mapmethod nnf --srcsubject fsaverage --trgsubject fsaverage5 --hemi %s --sval %s --tval %s', CommandInit,Hemis{iROI},['/data/',InName{iROI}],['/data/',OutName{iROI}]);
    system(Command);
end


%% ANCOVA - Seed-based FC - generate FC maps
% Seed-based FC using brain areas showed significant anatomical difference
clc;clear;
Date = '20230912_ICVCovariate';
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat');
InDir = '/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/ASD_SCZ/';
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults',filesep,Date,filesep,'SeedFC'];
mkdir(OutDir);
load('/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/ASD_SCZ/Beijing/ForDPABI_20210507/DPABISurf_AutoSave_2021_5_8_10_24.mat');

ID = SubInfo.ID;
Dir = SubInfo.Dir;

% Index = [find(Dx==1);find(Dx==2);find(Dx==3)]; 
Index = [1:length(ID)];%% calculate SCA for all subjects

Cfg.WorkingDir=InDir;
Cfg.OutDir=OutDir;
Cfg.StartingDirName='FunSurfWCF';
Cfg.StartingDirName_Volume='FunVoluWCF';
Cfg.FunctionalSessionNumber=1;
Cfg.SubjectNum=length(ID);

FunSessionPrefixSet={''}; %The first session doesn't need a prefix. From the second session, need a prefix such as 'S2_';
%     for iFunSession=2:Cfg.FunctionalSessionNumber
%         FunSessionPrefixSet=[FunSessionPrefixSet;{['S',num2str(iFunSession),'_']}];
%     end

%Make Subject Specific ROIs
for iFunSession=1:Cfg.FunctionalSessionNumber
    mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles']);
    for i=1:Cfg.SubjectNum
        SourceFile=[Dir{Index(i)},filesep,'fmriprep',filesep,ID{Index(i)},filesep,'ses-',num2str(iFunSession),filesep,'func',filesep,ID{Index(i)},'_ses-',num2str(iFunSession),'_task-rest_space-MNI152NLin2009cAsym_res-2_desc-aseg_dseg.nii.gz'];
        [Data Head]=y_Read(SourceFile);
        y_Write(Data==10,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LTha.nii']);
        y_Write(Data==17,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LHippo.nii']);
        y_Write(Data==53,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_RHippo.nii']);
        y_Write(Data==18,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LAmyg.nii']);
        y_Write(Data==54,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_RAmyg.nii']);
        disp(num2str(i))
    end
end


%Extract ROI Signals and Functional Connectivity Analysis
ROIDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/Anat/';
Cfg.CalFC.ROIDefSurfLH={[ROIDir,'Thickness_L_1_Mask_fs5.gii']};
Cfg.CalFC.ROIDefSurfRH={[ROIDir,'Thickness_R_1_Mask_fs5.gii'],...
    [ROIDir,'Thickness_R_2_Mask_fs5.gii']};

for iFunSession=1:Cfg.FunctionalSessionNumber
    mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'ROISignals_',Cfg.StartingDirName]);
    mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'ROISignals_',Cfg.StartingDirName]);
    if (Cfg.IsProcessVolumeSpace==1)
        mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'ROISignals_',Cfg.StartingDirName_Volume]);
    end
    mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName]);
    
    %Extract the ROI time courses
    parfor i=1:Cfg.SubjectNum
        ROISignalsSurfLH=[];
        ROISignalsSurfRH=[];
        ROISignalsVolu=[];
        % Left Hemi
        if ~isempty(Cfg.CalFC.ROIDefSurfLH)
            DirName=dir(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},'*fsaverage5_hemi-L_bold.func.gii'));
            for iFile=1:length(DirName)
                FileName=DirName(iFile).name;
                [ROISignalsSurfLH] = y_ExtractROISignal_Surf(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},FileName), ...
                    Cfg.CalFC.ROIDefSurfLH, ...
                    [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'ROISignals_',Cfg.StartingDirName,filesep,ID{Index(i)}], ...
                    '', ... % Will not restrict into the brain mask in extracting ROI signals
                    Cfg.CalFC.IsMultipleLabel);
            end
        end
        
        % Right Hemi
        if ~isempty(Cfg.CalFC.ROIDefSurfRH)
            DirName=dir(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},'*fsaverage5_hemi-R_bold.func.gii'));
            for iFile=1:length(DirName)
                FileName=DirName(iFile).name;
                [ROISignalsSurfRH] = y_ExtractROISignal_Surf(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},FileName), ...
                    Cfg.CalFC.ROIDefSurfRH, ...
                    [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'ROISignals_',Cfg.StartingDirName,filesep,ID{Index(i)}], ...
                    '', ... % Will not restrict into the brain mask in extracting ROI signals
                    Cfg.CalFC.IsMultipleLabel);
            end
        end
        
        % Volume
        
        ROIDefVolu={[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LTha.nii'];...
            [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LAmyg.nii'];...
            [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LHippo.nii'];...
            [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_RAmyg.nii'];...
            [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_RHippo.nii']};
%             [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_RTha.nii']};
%             [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_Optic.nii']};

            
        if ~isempty(ROIDefVolu)  % YAN Chao-Gan, 190708: if (Cfg.IsProcessVolumeSpace==1) && (~isempty(Cfg.CalFC.ROIDefVolu))
            [ROISignalsVolu] = y_ExtractROISignal([Dir{Index(i)},filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName_Volume,filesep,ID{Index(i)}], ...
                ROIDefVolu, ...
                [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'ROISignals_',Cfg.StartingDirName_Volume,filesep,ID{Index(i)}], ...
                '', ... % Will not restrict into the brain mask in extracting ROI signals
                Cfg.CalFC.IsMultipleLabel);
        end
        
        ROISignals = [ROISignalsSurfLH, ROISignalsSurfRH, ROISignalsVolu];
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROISignals_',ID{Index(i)},'.mat'], ROISignals, '');
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROISignals_',ID{Index(i)},'.txt'], ROISignals, ' ''-ASCII'', ''-DOUBLE'',''-TABS''');
        ROICorrelation = corrcoef(ROISignals);
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_',ID{Index(i)},'.mat'], ROICorrelation, '');
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_',ID{Index(i)},'.txt'], ROICorrelation, ' ''-ASCII'', ''-DOUBLE'',''-TABS''');
        ROICorrelation_FisherZ = 0.5 * log((1 + ROICorrelation)./(1- ROICorrelation));
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_FisherZ_',ID{Index(i)},'.mat'], ROICorrelation_FisherZ, '');
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_FisherZ_',ID{Index(i)},'.txt'], ROICorrelation_FisherZ, ' ''-ASCII'', ''-DOUBLE'',''-TABS''');
    end
end


%Calculate Seed Based Functional Connectivity
Cfg.MaskFileSurfLH='/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_lh_cortex.label.gii';
Cfg.MaskFileSurfRH='/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_rh_cortex.label.gii';
Cfg.MaskFileVolu='/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/Reslice_freesurfer_subcortical_mask_1mm.nii';
for iFunSession=1:Cfg.FunctionalSessionNumber
    mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'FC_SeedSurfLHSurfRHVolu_',Cfg.StartingDirName]);
    mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'FC_SeedSurfLHSurfRHVolu_',Cfg.StartingDirName]);
    if (Cfg.IsProcessVolumeSpace==1)
        mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'FC_SeedSurfLHSurfRHVolu_',Cfg.StartingDirName_Volume]);
    end
    
    parfor i=1:Cfg.SubjectNum
        % Calculate Functional Connectivity by Seed based Correlation Anlyasis
        % Left Hemi
        DirName=dir(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},'*fsaverage5_hemi-L_bold.func.gii'));
        for iFile=1:length(DirName)
            FileName=DirName(iFile).name;
            y_SCA_Surf(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},FileName), ...
                {[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROISignals_',ID{Index(i)},'.txt']}, ... %This is the ROI Signals extracted by the previous step
                [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'FC_SeedSurfLHSurfRHVolu_',Cfg.StartingDirName,filesep,'FC_',ID{Index(i)},'.func.gii'], ...
                Cfg.MaskFileSurfLH, ...
                1);
        end
        
        % Right Hemi
        DirName=dir(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},'*fsaverage5_hemi-R_bold.func.gii'));
        for iFile=1:length(DirName)
            FileName=DirName(iFile).name;
            y_SCA_Surf(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},FileName), ...
                {[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROISignals_',ID{Index(i)},'.txt']}, ...
                [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'FC_SeedSurfLHSurfRHVolu_',Cfg.StartingDirName,filesep,'FC_',ID{Index(i)},'.func.gii'], ...
                Cfg.MaskFileSurfRH, ...
                1);
        end
        
        % Volume
        if (Cfg.IsProcessVolumeSpace==1)
            % Calculate Functional Connectivity by Seed based Correlation Anlyasis
            y_SCA([Dir{Index(i)},filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName_Volume,filesep,ID{Index(i)}], ...
                {[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROISignals_',ID{Index(i)},'.txt']}, ...
                [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'FC_SeedSurfLHSurfRHVolu_',Cfg.StartingDirName_Volume,filesep,'FC_',ID{Index(i)},'.nii'], ...
                Cfg.MaskFileVolu, ...
                1);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%***Smooth the results***%%%%%%%%%%%%%%%%

Cfg.StartingDirName = 'Results';
Cfg.WorkingDir=Cfg.OutDir;
%Smooth on Results

%Get ready for later freesurfer usage.
[DPABIPath, fileN, extn] = fileparts(which('DPABI.m'));
ProgramPath=fullfile(DPABIPath, 'DPARSF');
TemplatePath=fullfile(DPABIPath, 'Templates');
if ispc
    CommandInit=sprintf('docker run -i --rm -v %s:/opt/freesurfer/license.txt -v %s:/data ', fullfile(DPABIPath, 'DPABISurf', 'FreeSurferLicense', 'license.txt'), Cfg.WorkingDir); %YAN Chao-Gan, 181214. Remove -t because there is a tty issue in windows
else
    CommandInit=sprintf('docker run -ti --rm -v %s:/opt/freesurfer/license.txt -v %s:/data ', fullfile(DPABIPath, 'DPABISurf', 'FreeSurferLicense', 'license.txt'), Cfg.WorkingDir); %YAN Chao-Gan, 181214. Remove -t because there is a tty issue in windows
end
if isdeployed % If running within docker with compiled version
    CommandInit=sprintf('export SUBJECTS_DIR=%s/freesurfer && ', Cfg.WorkingDir);
else
    CommandInit=sprintf('%s -e SUBJECTS_DIR=/data/freesurfer cgyan/dpabi', CommandInit);
end

mkdir([OutDir,filesep,'freesurfer/fsaverage5/']);
copyfile('/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/ASD_SCZ/Beijing/ForDPABI_20210507/freesurfer/fsaverage5/',[OutDir,filesep,'freesurfer/fsaverage5/']);

fprintf(['Smoothing the resutls...\n']);
FailList = []; %% very strange, random subject would fail smoothing while using parfor
for iFunSession=1:Cfg.FunctionalSessionNumber
    for i=1:Cfg.SubjectNum
        try
        %Check the DSpaces need to be normalized
        DirDSpace = dir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName]);
        if strcmpi(DirDSpace(3).name,'.DS_Store')  %110908 YAN Chao-Gan, for MAC OS compatablie
            StartIndex=4;
        else
            StartIndex=3;
        end
        DSpaceSet=[];
        for iDir=StartIndex:length(DirDSpace)
            if DirDSpace(iDir).isdir
                if ~((length(DirDSpace(iDir).name)>=28 && strcmpi(DirDSpace(iDir).name(1:28),'ROISignals_SurfLHSurfRHVolu_')))
                    DSpaceSet = [DSpaceSet;{DirDSpace(iDir).name}];
                end
            end
            
        end
        
        for iDSpace=1:length(DSpaceSet)
            switch DSpaceSet{iDSpace}
                case {'AnatSurfLH','AnatSurfRH'}
                    DirLevel1 = dir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace}]);
                    if strcmpi(DirLevel1(3).name,'.DS_Store')  %110908 YAN Chao-Gan, for MAC OS compatablie
                        StartIndexLevel1=4;
                    else
                        StartIndexLevel1=3;
                    end
                    for iDirLevel1=StartIndexLevel1:length(DirLevel1)
                        if DirLevel1(iDirLevel1).isdir
                            DirLevel2 = dir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace},filesep,DirLevel1(iDirLevel1).name]);
                            if strcmpi(DirLevel2(3).name,'.DS_Store')  %110908 YAN Chao-Gan, for MAC OS compatablie
                                StartIndexLevel2=4;
                            else
                                StartIndexLevel2=3;
                            end
                            for iDirLevel2=StartIndexLevel2:length(DirLevel2)
                                if DirLevel2(iDirLevel2).isdir
                                    SpaceName=DirLevel2(iDirLevel2).name;
                                    [tmp1 tmp2]=mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,'S',filesep,DSpaceSet{iDSpace},filesep,DirLevel1(iDirLevel1).name,filesep,SpaceName]);
                                    FileList=[];
                                    DirImg=dir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace},filesep,DirLevel1(iDirLevel1).name,filesep,SpaceName,filesep,'*',ID{Index(i)},'*.gii']);
                                    for j=1:length(DirImg)
                                        FileList=[FileList,' ',DirImg(j).name];
                                    end
                                    if strcmpi(DSpaceSet{iDSpace},'AnatSurfLH')
                                        Command = sprintf('%s parallel -j %g mri_surf2surf --s %s --hemi lh --sval /data/%s/AnatSurfLH/%s/%s/{1}  --fwhm %g --cortex --tval /data/%sS/AnatSurfLH/%s/%s/s{1} ::: %s', ...
                                            CommandInit, Cfg.ParallelWorkersNumber, SpaceName, [FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],DirLevel1(iDirLevel1).name, SpaceName, Cfg.Smooth.FWHMSurf, [FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],DirLevel1(iDirLevel1).name, SpaceName,FileList);
                                        system(Command);
                                        
                                    elseif strcmpi(DSpaceSet{iDSpace},'AnatSurfRH')
                                        Command = sprintf('%s parallel -j %g mri_surf2surf --s %s --hemi rh --sval /data/%s/AnatSurfRH/%s/%s/{1}  --fwhm %g --cortex --tval /data/%sS/AnatSurfRH/%s/%s/s{1} ::: %s', ...
                                            CommandInit, Cfg.ParallelWorkersNumber, SpaceName, [FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],DirLevel1(iDirLevel1).name, SpaceName, Cfg.Smooth.FWHMSurf, [FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],DirLevel1(iDirLevel1).name, SpaceName,FileList);
                                        system(Command);
                                    end
                                end
                            end
                        end
                    end
                    
                case {'FunSurfLH','FunSurfRH'}
                    SpaceName='fsaverage5';
                    %Check the measures need to be normalized
                    DirMeasure = dir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace}]);
                    if strcmpi(DirMeasure(3).name,'.DS_Store')  %110908 YAN Chao-Gan, for MAC OS compatablie
                        StartIndex=4;
                    else
                        StartIndex=3;
                    end
                    MeasureSet=[];
                    for iDir=StartIndex:length(DirMeasure)
                        if DirMeasure(iDir).isdir
                            if ~((length(DirMeasure(iDir).name)>=10 && strcmpi(DirMeasure(iDir).name(1:10),'ROISignals'))) %~((length(DirMeasure(iDir).name)>10 && strcmpi(DirMeasure(iDir).name(end-10:end),'_ROISignals')))
                                MeasureSet = [MeasureSet;{DirMeasure(iDir).name}];
                            end
                        end
                    end
                    for iMeasure=1:length(MeasureSet)
                        [tmp1 tmp2]=mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,'S',filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure}]);
                        FileList=[];
                        DirImg=dir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure},filesep,'*',ID{Index(i)},'*.gii']);
                        for j=1:length(DirImg)
                            FileList=[FileList,' ',DirImg(j).name];
                        end
                        if strcmpi(DSpaceSet{iDSpace},'FunSurfLH')
                            Command = sprintf('%s parallel -j %g mri_surf2surf --s %s --hemi lh --sval /data/%s/FunSurfLH/%s/{1}  --fwhm %g --cortex --tval /data/%sS/FunSurfLH/%s/s{1} ::: %s', ...
                                CommandInit, Cfg.ParallelWorkersNumber, SpaceName, [FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],MeasureSet{iMeasure}, Cfg.Smooth.FWHMSurf, [FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],MeasureSet{iMeasure},FileList);
                            system(Command);
                        elseif strcmpi(DSpaceSet{iDSpace},'FunSurfRH')
                            Command = sprintf('%s parallel -j %g mri_surf2surf --s %s --hemi rh --sval /data/%s/FunSurfRH/%s/{1}  --fwhm %g --cortex --tval /data/%sS/FunSurfRH/%s/s{1} ::: %s', ...
                                CommandInit, Cfg.ParallelWorkersNumber, SpaceName, [FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],MeasureSet{iMeasure}, Cfg.Smooth.FWHMSurf, [FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],MeasureSet{iMeasure},FileList);
                            system(Command);
                        end
                    end
                    
                case {'FunVolu','AnatVolu'}
                    %Check the measures need to be normalized
                    DirMeasure = dir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace}]);
                    if strcmpi(DirMeasure(3).name,'.DS_Store')  %110908 YAN Chao-Gan, for MAC OS compatablie
                        StartIndex=4;
                    else
                        StartIndex=3;
                    end
                    MeasureSet=[];
                    for iDir=StartIndex:length(DirMeasure)
                        if DirMeasure(iDir).isdir
                            if ~((length(DirMeasure(iDir).name)>=10 && strcmpi(DirMeasure(iDir).name(1:10),'ROISignals'))) %~((length(DirMeasure(iDir).name)>10 && strcmpi(DirMeasure(iDir).name(end-10:end),'_ROISignals')))
                                MeasureSet = [MeasureSet;{DirMeasure(iDir).name}];
                            end
                        end
                    end
                    for iMeasure=1:length(MeasureSet)
                        [tmp1 tmp2]=mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,'S',filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure}]);
                        FileList=[];
                        DirImg=dir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure},filesep,'*',ID{Index(i)},'*.nii.gz']);
                        if ~isempty(DirImg)
                            for j=1:length(DirImg)
                                gunzip([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure},filesep,DirImg(j).name]);
                                delete([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure},filesep,DirImg(j).name]);
                            end
                        end
                        DirImg=dir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure},filesep,'*',ID{Index(i)},'*.nii']);
                        if ~isempty(DirImg)
                            for j=1:length(DirImg)
                                FileList=[FileList;{[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure},filesep,DirImg(j).name]}];
                            end
                            SPMJOB = load([ProgramPath,filesep,'Jobmats',filesep,'Smooth.mat']);
                            SPMJOB.matlabbatch{1,1}.spm.spatial.smooth.data = FileList;
                            SPMJOB.matlabbatch{1,1}.spm.spatial.smooth.fwhm = Cfg.Smooth.FWHMVolu;
                            spm_jobman('run',SPMJOB.matlabbatch);
                            
                            DirTemp=dir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure},filesep,'ss*']);
                            if isempty(DirTemp)
                                movefile([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure},filesep,'s*'],[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,'S',filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure}]);
                            else
                                movefile([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure},filesep,'ss*'],[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName,'S',filesep,DSpaceSet{iDSpace},filesep,MeasureSet{iMeasure}]);
                            end
                        end
                    end
            end
        end
        catch
            FailList = [FailList,i];
        end
        disp(num2str(i))
    end
end
% tempList = zeros(size(SubInfo.ID));
% tempList(FailList) = 1;
% SubInfo.FailFCSmooth = tempList;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat','SubInfo','ROIValue')



%% ANCOVA - Seed-based FC - Combat
clc;clear;
Date = '20230912_ICVCovariate';
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat')
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SeedFC/ResultsS';
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults',filesep,Date,...
    filesep,'Combat_Metrics',filesep,'SeedFC'];
mkdir(OutDir);
Metrics  = {'FC_SeedSurfLHSurfRHVolu_FunSurfWCF','FC_SeedSurfLHSurfRHVolu_FunSurfWCF','FC_SeedSurfLHSurfRHVolu_FunVoluWCF'};
Hemispheres = {'FunSurfLH','FunSurfRH','FunVolu'};
Masks = {'/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_lh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_rh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/Reslice_freesurfer_subcortical_mask_1mm.nii'};

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Dir = SubInfo.Dir(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
Sex = SubInfo.Sex(Index);
Site = SubInfo.Site(Index);

AllCov_comBat = [Dx,Age,Sex];
nROI = 8; %%%%%%%%% Of note! &&&&&&&&
for iHemi = 1:length(Hemispheres)
    MaskData = y_ReadAll(Masks{iHemi});
    MaskDataOne = reshape(MaskData,[],1);
    MaskIndex = find(MaskDataOne);
    for iROI = 1:nROI
        OutputDir = [OutDir,filesep,Hemispheres{iHemi},filesep,'ROI',num2str(iROI)];
        mkdir(OutputDir);
        ImageList = cellfun(@(ID) dir([InDir,filesep,Hemispheres{iHemi},filesep,...
            Metrics{iHemi},filesep,'szROI',num2str(iROI),'FC_',ID,'*']),ID);
        ImageList = cellfun(@(Dir,ID) [Dir,filesep,ID],{ImageList(:).folder},{ImageList(:).name},'UniformOutput',false);
        if contains(Hemispheres{iHemi},'Surf')
            [Data, VoxelSize, FileList, Header] = y_ReadAll(ImageList);
            Data = Data(MaskIndex,:);
            Volume_Combat = combat(Data, Site', AllCov_comBat, 1);
            for iSub = 1:size(Volume_Combat,2)
                Brain = zeros(size(MaskData));
                Brain(MaskIndex) = Volume_Combat(:,iSub);
                y_Write(Brain,Header,[OutputDir,filesep,ID{iSub}]);
                disp([Hemispheres{iHemi},filesep,num2str(iROI),filesep,'Sub',num2str(iSub)])
            end
        elseif contains(Hemispheres{iHemi},'Volu')
            [Data, VoxelSize, FileList, Header] = y_ReadAll(ImageList);
            Data = reshape(Data,[],size(Data,4));
            Data = Data(MaskIndex,:);
            Volume_Combat = combat(Data, Site', AllCov_comBat, 1);
            for iSub = 1:size(Volume_Combat,2)
                BrainOne = zeros(size(MaskDataOne));
                BrainOne(MaskIndex) = Volume_Combat(:,iSub);
                Brain = reshape(BrainOne,size(MaskData,1),size(MaskData,2),size(MaskData,3));
                y_Write(Brain,Header,[OutputDir,filesep,ID{iSub}]);
                disp([Hemispheres{iHemi},filesep,num2str(iROI),filesep,'Sub',num2str(iSub)])
            end            
        end
    end
end



%% ANCOVA - Seed-based FC - Statistic 
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat');
% load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Templates/PALMSettings.mat');
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Templates/MonteCarloSettings.mat');
InDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/SeedFC'];
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SeedFC/StatisticResult'];
mkdir(OutDir);

Hemispheres = {'FunSurfLH','FunSurfRH','FunVolu'};
Suffix = {'.gii','.gii','.nii'};
nROI = 8; %%%%%%%%% change %%%%%%%%%%%
Masks = {'/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_lh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_rh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Templates/RemoveCerebellum_Reslice_freesurfer_subcortical_mask_1mm.nii'};
%     '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/Reslice_freesurfer_subcortical_mask_1mm.nii'};
% PALMSettings = {PALMSettingSurfL,PALMSettingSurfR,PALMSettingVolu};
MonteCarloSettings = {MonteCarloSurfL,MonteCarloSurfR,''};

for iHemi = 1:length(Hemispheres)
    parfor iROI = 1:nROI
        OutputDir = [OutDir,filesep,Hemispheres{iHemi},filesep,'ROI',num2str(iROI)];
        mkdir(OutputDir);
        
        ASDIndex = find(SubInfo.Dx==1 & SubInfo.QCFlag>0 & SubInfo.Age<=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);
        NCIndex = find(SubInfo.Dx==2 & SubInfo.QCFlag>0 & SubInfo.Age<=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);
        SCZIndex = find(SubInfo.Dx==3 & SubInfo.QCFlag>0 & SubInfo.Age<=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);
        
        ASDList = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,'ROI',num2str(iROI),...
            filesep,ID,Suffix{iHemi}],SubInfo.ID(ASDIndex),'UniformOutput',false);
        NCList = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,'ROI',num2str(iROI),...
            filesep,ID,Suffix{iHemi}],SubInfo.ID(NCIndex),'UniformOutput',false);
        SCZList = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,'ROI',num2str(iROI),...
            filesep,ID,Suffix{iHemi}],SubInfo.ID(SCZIndex),'UniformOutput',false);
        
        CovariateMatrix = {[SubInfo.Age(ASDIndex),SubInfo.Sex(ASDIndex),SubInfo.MotionFD(ASDIndex)];...
            [SubInfo.Age(NCIndex),SubInfo.Sex(NCIndex),SubInfo.MotionFD(NCIndex)];...
            [SubInfo.Age(SCZIndex),SubInfo.Sex(SCZIndex),SubInfo.MotionFD(SCZIndex)]};
        
        %         y_ANCOVA1_Image({ASDList,NCList,SCZList},[OutputDir,filesep,'ROI',num2str(iROI),'_AncovaPT'],Masks{iHemi},{},CovariateMatrix,PALMSettings{iHemi});
        y_ANCOVA1_Image({ASDList,NCList,SCZList},[OutputDir,filesep,'ROI',num2str(iROI),'_Ancova_MonteCarlo'],Masks{iHemi},{},CovariateMatrix);
        
%         % calculate Monte Carlo simulation table
%         if iHemi<3
%             [~,~,~,Header]=y_ReadAll([OutputDir,filesep,'ROI',num2str(iROI),'_Ancova_MonteCarlo.gii']);
%             HeaderInfo = Header.private.metadata(4).value;
%             FWHMIndex = strfind(HeaderInfo,'FWHM');
%             FWHM = str2num(HeaderInfo(FWHMIndex+5:FWHMIndex+9));
%             
%             SimReport=w_MonteCarlo_Surf_Bin(...
%                 {MonteCarloSettings{iHemi}.SurfPath},...
%                 FWHM,...
%                 0.001,... %%VertexP
%                 0.01667,... %%Alpha (ClusterP), 0.05/3
%                 1000,... % number of iteration
%                 OutputDir,...
%                 {MonteCarloSettings{iHemi}.MskFile},...
%                 {MonteCarloSettings{iHemi}.AreaFile});
%         end
        disp(['Have down ',Hemispheres{iHemi},' ROI',num2str(iROI)])
    end
end

% Record the peak F-value for drawing figure
% for iROI = 1:8
%     [BrainL,~,~,~]= y_ReadAll([OutDir,filesep,'FunSurfLH',filesep,'ROI',num2str(iROI),filesep,'Ancova.gii']);
%     [BrainR,~,~,~]= y_ReadAll([OutDir,filesep,'FunSurfRH',filesep,'ROI',num2str(iROI),filesep,'Ancova.gii']);
%     [BrainV,~,~,~]= y_ReadAll([OutDir,filesep,'FunVolu',filesep,'ROI',num2str(iROI),filesep,'Ancova.nii']);
%     MaxF(iROI,1) = max([max(BrainL),max(BrainR),max(reshape(BrainV,1,[]))]);
%     disp(['Max F value for ROI',num2str(iROI),' is ',num2str(MaxF(iROI,1))])
% end
% ROIValue.SeedFC.MaxFValue = MaxF;
% save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat','SubInfo','ROIValue')



%% ANCOVA - Seed-based FC - Sum up significant seed-based FC brain area
clc;clear;
InputList = {'/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/SeedFC/Corrected',...
    '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/SeedFC/Uncorrected'};
YeoL = '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V7.0_230110/DPABISurf/SurfTemplates/fsaverage5_lh_Yeo2011_7Networks_N1000.label.gii';
YeoR = '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V7.0_230110/DPABISurf/SurfTemplates/fsaverage5_rh_Yeo2011_7Networks_N1000.label.gii';


for i = 1:2
    InputDir = InputList{i};
    
    % Surf_L
    List = dir([InputDir,filesep,'FC_L_*.gii']);
    names = {List(:).name};
    List = cellfun(@(File) [InputDir,filesep,File],names,'UniformOutput',false)';
    [Data, VoxelSize, FileList, Header] = y_ReadAll(List);
    AllVolume = zeros(size(Data));
    Index = find(Data);
    AllVolume(Index) = 1;
    AllBrain = sum(AllVolume,2);
    y_Write(AllBrain,Header,[InputDir,filesep,'FC_Sum_L']);
    % add Yeo7 network
    [YeoData, ~, ~, YeoHeader] = y_ReadAll(YeoL);
    MaskIndex = find(AllBrain);
    AllBrain(MaskIndex) = YeoData(MaskIndex);
    y_Write(AllBrain,Header,[InputDir,filesep,'FC_Sum_L_Yeo7']);
    
    
    % Surf_R
    List = dir([InputDir,filesep,'FC_R_*.gii']);
    names = {List(:).name};
    List = cellfun(@(File) [InputDir,filesep,File],names,'UniformOutput',false)';
    [Data, VoxelSize, FileList, Header] = y_ReadAll(List);
    AllVolume = zeros(size(Data));
    Index = find(Data);
    AllVolume(Index) = 1;
    AllBrain = sum(AllVolume,2);
    y_Write(AllBrain,Header,[InputDir,filesep,'FC_Sum_R']);
    % add Yeo7 network
    [YeoData, ~, ~, YeoHeader] = y_ReadAll(YeoR);
    MaskIndex = find(AllBrain);
    AllBrain(MaskIndex) = YeoData(MaskIndex);
    y_Write(AllBrain,Header,[InputDir,filesep,'FC_Sum_R_Yeo7']);
end



%% ANCOVA - Regional metrics - Combat 
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat')
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/Regional_Metrics'];
mkdir(OutDir);
Metrics  = {'ALFF_FunSurfWC','fALFF_FunSurfWC','ReHo_FunSurfWCF','DegreeCentrality_FunSurfWCF'};
Prefixs = {'szALFF_','szfALFF_','szReHo_','szDegreeCentrality_Bilateral_PositiveWeightedSumBrain_'};
Hemispheres = {'FunSurfLH','FunSurfRH','FunVolu'};
Masks = {'/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_lh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_rh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/Reslice_freesurfer_subcortical_mask_1mm.nii'};

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Dir = SubInfo.Dir(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
Sex = SubInfo.Sex(Index);
Site = SubInfo.Site(Index);

AllCov_comBat = [Dx,Age,Sex];
for iHemi = 1:length(Hemispheres)
    MaskData = y_ReadAll(Masks{iHemi});
    MaskDataOne = reshape(MaskData,[],1);
    MaskIndex = find(MaskDataOne);
    for iMetric = 1:length(Metrics)
        OutputDir = [OutDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric}];
        mkdir(OutputDir);
        if contains(Hemispheres{iHemi},'Surf')
            ImageList = cellfun(@(Dir,ID) dir([Dir,filesep,'ResultsS',filesep,Hemispheres{iHemi},filesep,...
                Metrics{iMetric},filesep,Prefixs{iMetric},ID,'*']),Dir,ID);
            ImageList = cellfun(@(Dir,ID) [Dir,filesep,ID],{ImageList(:).folder},{ImageList(:).name},'UniformOutput',false);
            [Data, VoxelSize, FileList, Header] = y_ReadAll(ImageList);
            Data = Data(MaskIndex,:);
            Volume_Combat = combat(Data, Site', AllCov_comBat, 1);
            for iSub = 1:size(Volume_Combat,2)
                Brain = zeros(size(MaskData));
                Brain(MaskIndex) = Volume_Combat(:,iSub);
                y_Write(Brain,Header,[OutputDir,filesep,ID{iSub}]);
                disp([Hemispheres{iHemi},filesep,Metrics{iMetric},filesep,'Sub',num2str(iSub)])
            end
        elseif contains(Hemispheres{iHemi},'Volu')
            ImageList = cellfun(@(Dir,ID) dir([Dir,filesep,'ResultsS',filesep,Hemispheres{iHemi},filesep,...
                regexprep(Metrics{iMetric},'Surf','Volu'),filesep,regexprep(Prefixs{iMetric},'_Bilateral_','_'),ID,'*']),Dir,ID);
            ImageList = cellfun(@(Dir,ID) [Dir,filesep,ID],{ImageList(:).folder},{ImageList(:).name},'UniformOutput',false);
            [Data, VoxelSize, FileList, Header] = y_ReadAll(ImageList);
            Data = reshape(Data,[],size(Data,4));
            Data = Data(MaskIndex,:);
            Volume_Combat = combat(Data, Site', AllCov_comBat, 1);
            for iSub = 1:size(Volume_Combat,2)
                BrainOne = zeros(size(MaskDataOne));
                BrainOne(MaskIndex) = Volume_Combat(:,iSub);
                Brain = reshape(BrainOne,size(MaskData,1),size(MaskData,2),size(MaskData,3));
                y_Write(Brain,Header,[OutputDir,filesep,ID{iSub}]);
                disp([Hemispheres{iHemi},filesep,Metrics{iMetric},filesep,'Sub',num2str(iSub)])
            end            
        end
    end
end



%% ANCOVA - Regional metrics - Statistic
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat');
% load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Templates/PALMSettings.mat');
% load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Templates/MonteCarloSettings.mat');
InDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/Regional_Metrics'];
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Regional_Metrics'];
mkdir(OutDir);

Hemispheres = {'FunSurfLH','FunSurfRH','FunVolu'};
Metrics  = {'ALFF_FunSurfWC','fALFF_FunSurfWC','ReHo_FunSurfWCF','DegreeCentrality_FunSurfWCF'};
Suffix = {'.gii','.gii','.nii'};
Masks = {'/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_lh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_rh_cortex.label.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Templates/RemoveCerebellum_Reslice_freesurfer_subcortical_mask_1mm.nii'};
%     '/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/Reslice_freesurfer_subcortical_mask_1mm.nii'};
% PALMSettings = {PALMSettingSurfL,PALMSettingSurfR,PALMSettingVolu};
% MonteCarloSettings = {MonteCarloSurfL,MonteCarloSurfR,''};

for iHemi = 1:length(Hemispheres)
    parfor iMetric = 1:length(Metrics)
        OutputDir = [OutDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric}];
        mkdir(OutputDir);
        
        ASDIndex = find(SubInfo.Dx==1 & SubInfo.QCFlag>0 & SubInfo.Age<=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);
        NCIndex = find(SubInfo.Dx==2 & SubInfo.QCFlag>0 & SubInfo.Age<=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);
        SCZIndex = find(SubInfo.Dx==3 & SubInfo.QCFlag>0 & SubInfo.Age<=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);
        
        ASDList = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric},...
            filesep,ID,Suffix{iHemi}],SubInfo.ID(ASDIndex),'UniformOutput',false);
        NCList = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric},...
            filesep,ID,Suffix{iHemi}],SubInfo.ID(NCIndex),'UniformOutput',false);
        SCZList = cellfun(@(ID) [InDir,filesep,Hemispheres{iHemi},filesep,Metrics{iMetric},...
            filesep,ID,Suffix{iHemi}],SubInfo.ID(SCZIndex),'UniformOutput',false);
        
        CovariateMatrix = {[SubInfo.Age(ASDIndex),SubInfo.Sex(ASDIndex),SubInfo.MotionFD(ASDIndex)];...
            [SubInfo.Age(NCIndex),SubInfo.Sex(NCIndex),SubInfo.MotionFD(NCIndex)];...
            [SubInfo.Age(SCZIndex),SubInfo.Sex(SCZIndex),SubInfo.MotionFD(SCZIndex)]};
        
        %         y_ANCOVA1_Image({ASDList,NCList,SCZList},[OutputDir,filesep,'ROI',num2str(iMetric),'_AncovaPT'],Masks{iHemi},{},CovariateMatrix,PALMSettings{iHemi});
        y_ANCOVA1_Image({ASDList,NCList,SCZList},[OutputDir,filesep,Metrics{iMetric},'_Ancova'],Masks{iHemi},{},CovariateMatrix);
        
%         % calculate Monte Carlo simulation table
%         if iHemi<3
%             [~,~,~,Header]=y_ReadAll([OutputDir,filesep,'ROI',num2str(iMetric),'_Ancova_MonteCarlo.gii']);
%             HeaderInfo = Header.private.metadata(4).value;
%             FWHMIndex = strfind(HeaderInfo,'FWHM');
%             FWHM = str2num(HeaderInfo(FWHMIndex+5:FWHMIndex+9));
%             
%             SimReport=w_MonteCarlo_Surf_Bin(...
%                 {MonteCarloSettings{iHemi}.SurfPath},...
%                 FWHM,...
%                 0.001,... %%VertexP
%                 0.01667,... %%Alpha (ClusterP), 0.05/3
%                 1000,... % number of iteration
%                 OutputDir,...
%                 {MonteCarloSettings{iHemi}.MskFile},...
%                 {MonteCarloSettings{iHemi}.AreaFile});
%         end
        disp(['Have down ',Hemispheres{iHemi},' metric',num2str(iMetric)])
    end
end



%% ANCOVA - ROI-ROI FC - generate FC matrix
% Extract anatomical-ROI （e.g. surface thickness, subcortical volume) time series and calculate ROI-ROI FC
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat');
InDir = '/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/ASD_SCZ/';
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_ROI2ROI'];
mkdir(OutDir);
load('/mnt/Data6/RfMRILab/Lubin/DataPreprocessing/ASD_SCZ/Beijing/ForDPABI_20210507/DPABISurf_AutoSave_2021_5_8_10_24.mat');

ID = SubInfo.ID;
Dir = SubInfo.Dir;
Site = SubInfo.Site;
Sex = SubInfo.Sex;
Dx = SubInfo.Dx;
Age = SubInfo.Age;

% Index = [find(Dx==1);find(Dx==2);find(Dx==3)];
Index = [1:length(Dx)]'; %% calculate SCA for all subjects

Cfg.WorkingDir=InDir;
Cfg.OutDir=OutDir;

Cfg.StartingDirName='FunSurfWCFS';
Cfg.StartingDirName_Volume='FunVoluWCFS';

Cfg.FunctionalSessionNumber=1;
Cfg.SubjectNum=length(Index);

FunSessionPrefixSet={''}; %The first session doesn't need a prefix. From the second session, need a prefix such as 'S2_';
%     for iFunSession=2:Cfg.FunctionalSessionNumber
%         FunSessionPrefixSet=[FunSessionPrefixSet;{['S',num2str(iFunSession),'_']}];
%     end

%Make Subject Specific ROIs
% 
for iFunSession=1:Cfg.FunctionalSessionNumber
    mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles']);
    for i=1:Cfg.SubjectNum
        SourceFile=[Dir{Index(i)},filesep,'fmriprep',filesep,ID{Index(i)},filesep,'ses-',num2str(iFunSession),filesep,'func',filesep,ID{Index(i)},'_ses-',num2str(iFunSession),'_task-rest_space-MNI152NLin2009cAsym_res-2_desc-aseg_dseg.nii.gz'];
        [Data Head]=y_Read(SourceFile);
        y_Write(Data==10,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LTha.nii']);
%         y_Write(Data==49,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_RTha.nii']);
        y_Write(Data==17,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LHippo.nii']);
        y_Write(Data==53,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_RHippo.nii']);
        y_Write(Data==18,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LAmyg.nii']);
        y_Write(Data==54,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_RAmyg.nii']);
%         y_Write(Data==85,Head,[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_Optic.nii']);
        disp(num2str(i))
    end
end


%Extract ROI Signals and Functional Connectivity Analysis
ROIDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/Anat/';
Cfg.CalFC.ROIDefSurfLH={[ROIDir,'Thickness_L_1_Mask_fs5.gii']};
Cfg.CalFC.ROIDefSurfRH={[ROIDir,'Thickness_R_1_Mask_fs5.gii'],...
    [ROIDir,'Thickness_R_2_Mask_fs5.gii']};

for iFunSession=1:Cfg.FunctionalSessionNumber
    mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'ROISignals_',Cfg.StartingDirName]);
    mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'ROISignals_',Cfg.StartingDirName]);
    if (Cfg.IsProcessVolumeSpace==1)
        mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'ROISignals_',Cfg.StartingDirName_Volume]);
    end
    mkdir([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName]);
    
    %Extract the ROI time courses
    parfor i=1:Cfg.SubjectNum
        ROISignalsSurfLH=[];
        ROISignalsSurfRH=[];
        ROISignalsVolu=[];
        % Left Hemi
        if ~isempty(Cfg.CalFC.ROIDefSurfLH)
            DirName=dir(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},'*fsaverage5_hemi-L_bold.func.gii'));
            for iFile=1:length(DirName)
                FileName=DirName(iFile).name;
                [ROISignalsSurfLH] = y_ExtractROISignal_Surf(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},FileName), ...
                    Cfg.CalFC.ROIDefSurfLH, ...
                    [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'ROISignals_',Cfg.StartingDirName,filesep,ID{Index(i)}], ...
                    '', ... % Will not restrict into the brain mask in extracting ROI signals
                    Cfg.CalFC.IsMultipleLabel);
            end
        end
        
        % Right Hemi
        if ~isempty(Cfg.CalFC.ROIDefSurfRH)
            DirName=dir(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},'*fsaverage5_hemi-R_bold.func.gii'));
            for iFile=1:length(DirName)
                FileName=DirName(iFile).name;
                [ROISignalsSurfRH] = y_ExtractROISignal_Surf(fullfile(Dir{Index(i)},[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],ID{Index(i)},FileName), ...
                    Cfg.CalFC.ROIDefSurfRH, ...
                    [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'ROISignals_',Cfg.StartingDirName,filesep,ID{Index(i)}], ...
                    '', ... % Will not restrict into the brain mask in extracting ROI signals
                    Cfg.CalFC.IsMultipleLabel);
            end
        end
        
        % Volume
        
        ROIDefVolu={[Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LTha.nii'];...
            [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LHippo.nii'];...
            [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_LAmyg.nii'];...
            [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_RHippo.nii'];...
            [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROIFiles',filesep,ID{Index(i)},'_RAmyg.nii']};
        
        if ~isempty(ROIDefVolu)  % YAN Chao-Gan, 190708: if (Cfg.IsProcessVolumeSpace==1) && (~isempty(Cfg.CalFC.ROIDefVolu))
            [ROISignalsVolu] = y_ExtractROISignal([Dir{Index(i)},filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName_Volume,filesep,ID{Index(i)}], ...
                ROIDefVolu, ...
                [Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'ROISignals_',Cfg.StartingDirName_Volume,filesep,ID{Index(i)}], ...
                '', ... % Will not restrict into the brain mask in extracting ROI signals
                Cfg.CalFC.IsMultipleLabel);
        end
        
        ROISignals = [ROISignalsSurfLH, ROISignalsSurfRH, ROISignalsVolu];
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROISignals_',ID{Index(i)},'.mat'], ROISignals, '');
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROISignals_',ID{Index(i)},'.txt'], ROISignals, ' ''-ASCII'', ''-DOUBLE'',''-TABS''');
        ROICorrelation = corrcoef(ROISignals);
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_',ID{Index(i)},'.mat'], ROICorrelation, '');
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_',ID{Index(i)},'.txt'], ROICorrelation, ' ''-ASCII'', ''-DOUBLE'',''-TABS''');
        ROICorrelation_FisherZ = 0.5 * log((1 + ROICorrelation)./(1- ROICorrelation));
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_FisherZ_',ID{Index(i)},'.mat'], ROICorrelation_FisherZ, '');
        y_CallSave([Cfg.OutDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_FisherZ_',ID{Index(i)},'.txt'], ROICorrelation_FisherZ, ' ''-ASCII'', ''-DOUBLE'',''-TABS''');
    end
end



%% ANCOVA - ROI-ROI FC - Combat 
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
InDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_ROI2ROI/Results/ROISignals_SurfLHSurfRHVolu_FunSurfWCFS/'];
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/FC_ROI2ROI'];
mkdir(OutDir);

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Dir = SubInfo.Dir(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
Sex = SubInfo.Sex(Index);
Site = SubInfo.Site(Index);

nROI = 8;
Mask = tril(ones(nROI,nROI),-1);
MaskIndex = find(Mask);

AllCov_comBat = [Dx,Age,Sex];
for iSub = 1:length(ID)
    load([InDir,filesep,'ROICorrelation_FisherZ_',ID{iSub},'.mat']);
    DataOneDim = reshape(Data(MaskIndex),1,[]);
    FCAll(:,iSub) = DataOneDim;
end

Volume_Combat = combat(FCAll, Site', AllCov_comBat, 1);

for iSub = 1:size(Volume_Combat,2)
    CorrMatrix = zeros(nROI,nROI);
    CorrMatrix(MaskIndex) = Volume_Combat(:,iSub);
    save([OutDir,filesep,'ROICorrelation_FisherZ_',ID{iSub},'.mat'],'CorrMatrix');
end
       


%% ANCOVA - ROI-ROI FC - Statistic 
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat');
InDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/FC_ROI2ROI'];

ID = SubInfo.ID;
Dir = SubInfo.Dir;
Site = SubInfo.Site;
Sex = SubInfo.Sex;
Dx = SubInfo.Dx;
Age = SubInfo.Age;
MotionFD = SubInfo.MotionFD;

nROI = 8; %%%%%%
F_FC_ROI2ROI = zeros(nROI,nROI);
P_FC_ROI2ROI = zeros(nROI,nROI);


DxCov=[];
DxIndex = [1,2,3];
for i=1:length(DxIndex)-1
    DxCov(:,i) = Dx==DxIndex(i);
end

FCAll = zeros(length(ID),nROI,nROI);
for i = 1:length(ID)
    load([InDir,filesep,'ROICorrelation_FisherZ_',ID{i},'.mat']);
    FCAll(i,:,:) = CorrMatrix;
end

for iROI = 1:nROI
    for jROI = 1:nROI
        AllCov = [DxCov,ones(length(Dx),1),Age,Sex,MotionFD];
        Contrast=zeros(1,size(AllCov,2));
        Contrast(1:size(DxCov,2))=1;
        if iROI~=jROI
            [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(FCAll(:,iROI,jROI),AllCov,Contrast,'F');
            F_FC_ROI2ROI(iROI,jROI) = TF_ForContrast;
            P_FC_ROI2ROI(iROI,jROI) = 1-fcdf(TF_ForContrast,size(DxCov,2),size(DxCov,1)-size(AllCov,2));
        end
    end
end

P_Corrected = mafdr(P_FC_ROI2ROI(find(tril(P_FC_ROI2ROI))),'BHFDR','true');

P_FC_ROI2ROI_Corrected = zeros(size(P_FC_ROI2ROI));
Index = find(tril(P_FC_ROI2ROI));
P_FC_ROI2ROI_Corrected(Index) = P_Corrected;
[SigIndex1,SigIndex2] = find(P_FC_ROI2ROI_Corrected<0.05 & P_FC_ROI2ROI_Corrected>0);
ROILabels = {'IPS_L','IPS_R','TPJ_R','Tha_L','Hippo_L','Amyg_L','Hippo_R','Amyg_R'}; %%%%%
ROILabel_Sig = {};
for i =1:length(SigIndex1)
    ROILabel_Sig{i} = [num2str(i),': ',ROILabels{SigIndex2(i)},' - ',ROILabels{SigIndex1(i)}, ', F: ',num2str(F_FC_ROI2ROI(SigIndex1(i),SigIndex2(i))), ', p: ',num2str(P_FC_ROI2ROI_Corrected(SigIndex1(i),SigIndex2(i)))];
    disp(ROILabel_Sig{i})
end

ROIFC_Sig = zeros(length(ID),length(SigIndex1));
for i = 1:length(SigIndex1)
    ROIFC_Sig(:,i) = FCAll(:,SigIndex1(i),SigIndex2(i));
end

ROIValue.ROI2ROIFC.AncovaF = F_FC_ROI2ROI;
ROIValue.ROI2ROIFC.AncovaP = P_FC_ROI2ROI;
ROIValue.ROI2ROIFC.AncovaP_Corrected = P_FC_ROI2ROI_Corrected;
ROIValue.ROI2ROIFC.FC_Sig = ROIFC_Sig;
ROIValue.ROI2ROIFC.Label_Sig = ROILabel_Sig;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','ROIValue','SubInfo');



%% ANCOVA - Graph theoretical metrics - Combat preparing
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20211119_ICVCovariate/NetworkAnalysis_GraphTheory_Revised';
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/GraphTheory'];
mkdir(OutDir);

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Dir = SubInfo.Dir(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
Sex = SubInfo.Sex(Index);
Site = SubInfo.Site(Index);

FieldList1 = {'Lp_AUC','Cp_AUC','Eloc_AUC','Eglob_AUC',...
    'Assortativity_AUC','Modularity_AUC','ClusteringCoefficient_AUC'}; % One value for each subject 
FieldList2 = {'Degree_AUC','NodalEfficiency_AUC','Betweenness_AUC','ParticipantCoefficient_AUC',...
    'SubgraphCentrality_AUC','EigenvectorCentrality_AUC','PageRankCentrality_AUC','ElocSet','EglobSet'}; % One value for each node/step (multiple values for each subject)

AllValue = struct();
GTAFlag = ones(length(ID),1);
for iSub = 1:length(ID)
    try
        GTA_Value=load([InDir,filesep,'GTA_',ID{iSub},'.mat']);
    catch
        GTAFlag(iSub) =0;
    end
    for iField = 1:length(FieldList1)
        if ~isfield(AllValue,FieldList1{iField})
            AllValue.(FieldList1{iField}) = GTA_Value.(FieldList1{iField});
        else
            AllValue.(FieldList1{iField}) = [AllValue.(FieldList1{iField}),GTA_Value.(FieldList1{iField})];
        end
    end
    
    for iField = 1:length(FieldList2)
        GTA_Value.(FieldList2{iField}) = reshape(GTA_Value.(FieldList2{iField}),[],1);
        if ~isfield(AllValue,FieldList2{iField})
            AllValue.(FieldList2{iField}) = GTA_Value.(FieldList2{iField});
        else
            AllValue.(FieldList2{iField}) = [AllValue.(FieldList2{iField}),GTA_Value.(FieldList2{iField})];
        end
    end
    disp(num2str(iSub))
end

%% ANCOVA - Graph theoretical metrics - (Visual check) 
for iMetric = 1:length(FieldList1)
    figure,hist(AllValue.(FieldList1{iMetric}))
    title(FieldList1{iMetric});
end


for iMetric = 1:length(FieldList2)
    figure,hist(mean(AllValue.(FieldList2{iMetric}),1))
    title([FieldList2{iMetric},filesep,'Across ROI']);
    figure,hist(mean(AllValue.(FieldList2{iMetric}),2))
    title([FieldList2{iMetric},filesep,'Across Subject']);
end

%% ANCOVA - Graph theoretical metrics - (Deal with outliers)
% SubgraphCentrality_AUC - across ROI mean
Mean_Sub = mean(AllValue.SubgraphCentrality_AUC,1);
Index = find(((Mean_Sub-mean(Mean_Sub))./std(Mean_Sub))>3);
GTAFlag(Index)=0;

%% ANCOVA - Graph theoretical metrics - Combat 
AllCov_comBat = [Dx,Age,Sex];
GTAIndex = find(GTAFlag);
GTA_Combat = struct();
Value_Raw = []; 

% FieldList1
for iField = 1:length(FieldList1)
    Value_Raw = [Value_Raw;AllValue.(FieldList1{iField})];
end
% Value_Raw = (Value_Raw-repmat(mean(Value_Raw,2),1,size(Value_Raw,2)))./repmat(std(Value_Raw,0,2),1,size(Value_Raw,2));
Value_Combat = zeros(size(Value_Raw));
Value_Combat(:,GTAIndex) = combat(Value_Raw(:,GTAIndex), Site(GTAIndex)', AllCov_comBat(GTAIndex,:), 1);
for iField = 1:length(FieldList1)
    GTA_Combat.(FieldList1{iField}) = Value_Combat(iField,:);
end

% FieldList2
for iField = 1:length(FieldList2)
    Value_Raw = AllValue.(FieldList2{iField});
    Value_Combat = zeros(size(Value_Raw));
    Value_Combat(:,GTAIndex) = combat(Value_Raw(:,GTAIndex), Site(GTAIndex)', AllCov_comBat(GTAIndex,:), 1);
    GTA_Combat.(FieldList2{iField}) = Value_Combat;
end

ROIValue.Graph.ValueRaw = AllValue;
ROIValue.Graph.ValueCombat = GTA_Combat;
ROIValue.Graph.FieldList1 = FieldList1;
ROIValue.Graph.FieldList2 = FieldList2;
SubInfo.GTAFlag = GTAFlag;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','ROIValue','SubInfo');



%% ANCOVA - Graph theoretical metrics - Statistic 
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat');

%%% Chekc GTAFlag! %%%
Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2 & SubInfo.GTAFlag);

ID = SubInfo.ID(Index);
Dx = SubInfo.Dx(Index);
Site = SubInfo.Site(Index);
Age = SubInfo.Age(Index);
Sex = SubInfo.Sex(Index);
MotionFD = SubInfo.MotionFD(Index);

DxCov=[];
DxIndex = [1,2,3];
for i=1:length(DxIndex)-1
    DxCov(:,i) = Dx==DxIndex(i);
end  

SiteCov=[];
SiteIndex = unique(Site);
for i=1:length(SiteIndex)-1
    SiteCov(:,i) = Site==SiteIndex(i);
end

FieldList1 = ROIValue.Graph.FieldList1;
FieldList2 = ROIValue.Graph.FieldList2;

%%% Single values should not use combat results %%%
AllCov = [DxCov,ones(length(Dx),1),Age,Sex,MotionFD,SiteCov];
Contrast=zeros(1,size(AllCov,2));
Contrast(1:size(DxCov,2))=1;

for iMetric = 1:length(FieldList1)
    GTA_Value = ROIValue.Graph.ValueRaw.(FieldList1{iMetric});
    GTA_Value = GTA_Value(:,Index)';
    RegressedMetrics.(FieldList1{iMetric}) = zeros(1,length(ROIValue.Graph.ValueRaw.(FieldList1{iMetric})));
    
    [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(GTA_Value,AllCov,Contrast,'F');
    Stat_GTA_ANCOVA.F.(FieldList1{iMetric}) = TF_ForContrast;
    Stat_GTA_ANCOVA.P.(FieldList1{iMetric}) = 1-fcdf(TF_ForContrast,size(DxCov,2),length(Age)-size(AllCov,2));
    % save metrics after regressing out covariates
    [b,r,SSE,SSR] = y_regress_ss(GTA_Value,AllCov(:,2:end));
    RegressedMetrics.(FieldList1{iMetric})(Index)=b(1).*ones(length(Age),1)+r;
end

%%% Dependent values should not use combat results - EglobSet/ElocSet %%%
for iMetric = 8:9 % EglobSet/ElocSet
    GTA_Value = ROIValue.Graph.ValueRaw.(FieldList2{iMetric});
    GTA_Value = GTA_Value(:,Index)';
    RegressedMetrics.(FieldList2{iMetric}) = zeros(size(ROIValue.Graph.ValueRaw.(FieldList2{iMetric})));

    for i = 1:size(GTA_Value,2)
        [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(GTA_Value(:,i),AllCov,Contrast,'F');
        Stat_GTA_ANCOVA.F.(FieldList2{iMetric})(i) = TF_ForContrast;
        Stat_GTA_ANCOVA.P.(FieldList2{iMetric})(i) = 1-fcdf(TF_ForContrast,size(DxCov,2),length(Age)-size(AllCov,2));
        % save metrics after regressing out covariates
        [b,r,SSE,SSR] = y_regress_ss(GTA_Value(:,i),AllCov(:,2:end));
        RegressedMetrics.(FieldList2{iMetric})(i,Index) = b(1).*ones(length(Age),1)+r;
    end
    Stat_GTA_ANCOVA.P_Corrected.(FieldList2{iMetric}) = mafdr(Stat_GTA_ANCOVA.P.(FieldList2{iMetric}),'BHFDR','true');
end

%%% Multiple-ROI values should use combat results %%%
AllCov = [DxCov,ones(length(Dx),1),Age,Sex,MotionFD];
Contrast=zeros(1,size(AllCov,2));
Contrast(1:size(DxCov,2))=1;

for iMetric = 1:length(FieldList2)-2 % no EglobSet/ElocSet
    GTA_Value = ROIValue.Graph.ValueCombat.(FieldList2{iMetric});
    GTA_Value = GTA_Value(:,Index)';
    RegressedMetrics.(FieldList2{iMetric}) = zeros(size(ROIValue.Graph.ValueRaw.(FieldList2{iMetric})));

    for i = 1:size(GTA_Value,2)
        [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(GTA_Value(:,i),AllCov,Contrast,'F');
        Stat_GTA_ANCOVA.F.(FieldList2{iMetric})(i) = TF_ForContrast;
        Stat_GTA_ANCOVA.P.(FieldList2{iMetric})(i) = 1-fcdf(TF_ForContrast,size(DxCov,2),length(Age)-size(AllCov,2));
        % save metrics after regressing out covariates
        [b,r,SSE,SSR] = y_regress_ss(GTA_Value(:,i),AllCov(:,2:end));
        RegressedMetrics.(FieldList2{iMetric})(i,Index) = b(1).*ones(length(Age),1)+r;
    end
    Stat_GTA_ANCOVA.P_Corrected.(FieldList2{iMetric}) = mafdr(Stat_GTA_ANCOVA.P.(FieldList2{iMetric}),'BHFDR','true');
end

for iMetric = 2:2 % NodalEfficiency
    GTA_Value = ROIValue.Graph.ValueCombat.(FieldList2{iMetric});
    GTA_Value = GTA_Value(:,Index)';
    RegressedMetrics.(FieldList2{iMetric}) = zeros(size(ROIValue.Graph.ValueRaw.(FieldList2{iMetric})));

    for i = 1:size(GTA_Value,2)
        [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(GTA_Value(:,i),AllCov,Contrast,'F');
        Stat_GTA_ANCOVA.F.(FieldList2{iMetric})(i) = TF_ForContrast;
        Stat_GTA_ANCOVA.P.(FieldList2{iMetric})(i) = 1-fcdf(TF_ForContrast,size(DxCov,2),length(Age)-size(AllCov,2));
        % save metrics after regressing out covariates
        [b,r,SSE,SSR] = y_regress_ss(GTA_Value(:,i),AllCov(:,2:end));
        RegressedMetrics.(FieldList2{iMetric})(i,Index) = b(1).*ones(length(Age),1)+r;
    end
    
    GradIndex1 = find(ROIValue.Gradient.Ancova.P_Grad_ANCOVA_Corrected(1,:)<0.05 | ROIValue.Gradient.Ancova.P_Grad_ANCOVA_Corrected(2,:)<0.05);
    Stat_GTA_ANCOVA.P_Corrected.NodalEffi_Gradient(GradIndex1) = mafdr(Stat_GTA_ANCOVA.P.(FieldList2{iMetric})(GradIndex1),'BHFDR','true');
    GradIndex2 = find(ROIValue.Gradient.Ancova.P_Grad_ANCOVA_Corrected(1,:)>0.05 & ROIValue.Gradient.Ancova.P_Grad_ANCOVA_Corrected(2,:)>0.05);
    Stat_GTA_ANCOVA.P_Corrected.NodalEffi_Gradient(GradIndex2) = mafdr(Stat_GTA_ANCOVA.P.(FieldList2{iMetric})(GradIndex2),'BHFDR','true');
end


ROIValue.Graph.Ancova = Stat_GTA_ANCOVA;
ROIValue.Graph.ValueCovRegressed = RegressedMetrics;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','ROIValue','SubInfo');

%% ANCOVA - Graph theoretical metrics - (Inspect SigROI)
clc;clear; 
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat');
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/DPABISurf_Schaefer2018_400_Tian2020_54_Info.mat');

for i = 1:length(ROIValue.Graph.FieldList2)
    SigIndex = find(ROIValue.Graph.Ancova.P_Corrected.(ROIValue.Graph.FieldList2{i})<0.05);
    disp([ROIValue.Graph.FieldList2{i},': ',num2str(length(SigIndex))]);
    if i <= length(ROIValue.Graph.FieldList2)-2
        [DPABISurf_Schaefer2018_400_Tian2020_54_NodeName{SigIndex}]
    end
end



%% ANCOVA - Connetome - Combat 
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20211119_ICVCovariate/NetworkAnalysis_FCMatrix_UseGUI';
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/FC_Network'];
mkdir(OutDir);

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Dir = SubInfo.Dir(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
Sex = SubInfo.Sex(Index);
Site = SubInfo.Site(Index);

nROI = 454;
Mask = tril(ones(nROI,nROI),-1);
MaskIndex = find(Mask);

AllCov_comBat = [Dx,Age,Sex];
for iSub = 1:length(ID)
    load([InDir,filesep,'NetworkMatrix_',ID{iSub},'.mat']);
    DataOneDim = reshape(NetworkMatrix(MaskIndex),1,[]);
    FCAll(:,iSub) = DataOneDim;
end

Volume_Combat = combat(FCAll, Site', AllCov_comBat, 1);

for iSub = 1:size(Volume_Combat,2)
    CorrMatrix = zeros(nROI,nROI);
    CorrMatrix(MaskIndex) = Volume_Combat(:,iSub);
    CorrMatrix = CorrMatrix + CorrMatrix';
    save([OutDir,filesep,'NetworkMatrix_',ID{iSub},'.mat'],'CorrMatrix');
    disp(num2str(iSub))
end
       


%% ANCOVA - Connetome - statistic
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat');
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/FC_Network';
OutDir= '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_Network';
mkdir(OutDir);

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Sex = SubInfo.Sex(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
MotionFD = SubInfo.MotionFD(Index);

ASDIndex = find(Dx==1); NCIndex = find(Dx==2); SCZIndex = find(Dx==3);

ASDList = cellfun(@(ID) [InDir,filesep,'NetworkMatrix_',ID,'.mat'],ID(ASDIndex),'UniformOutput',false);
NCList = cellfun(@(ID) [InDir,filesep,'NetworkMatrix_',ID,'.mat'],ID(NCIndex),'UniformOutput',false);
SCZList = cellfun(@(ID) [InDir,filesep,'NetworkMatrix_',ID,'.mat'],ID(SCZIndex),'UniformOutput',false);

CovariateMatrix = {[Age(ASDIndex),Sex(ASDIndex),MotionFD(ASDIndex)];...
    [Age(NCIndex),Sex(NCIndex),MotionFD(NCIndex)];...
    [Age(SCZIndex),Sex(SCZIndex),MotionFD(SCZIndex)]};

y_ANCOVA1_Image({ASDList,NCList,SCZList},[OutDir,filesep,'FC_Ancova'],'',{},CovariateMatrix);


%%%%%%%% Plot in DPABINet_VIEW %%%%%%%%
%%%%%%%% Save FMatrix_Ancova_FDR_0001.mat %%%%%%%%


%% ANCOVA - Connetome - (Revision template file for DPABINet_VIEW)
clc;clear;
Raw = load('/mnt/Data4/RfMRILab/Lubin/Project/RMD_Cross_Diagnosis/StatisticalResults/NetworkAnalysis_20220803/Short_DPABISurf_Schaefer2018_400_Tian2020_54_Info.mat');
OutDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_Network';
Prefix = 'Template_forDPABINetVIEW_';

Thresholds = {'005','01'};
for iThre = 1:length(Thresholds)
    load(['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_Network/EdgeMatrix_Ancova_FDR_Dot',Thresholds{iThre},'.mat']);
    Index = find(any(EdgeMatrix));
    NodeLabel = Raw.DPABISurf_Schaefer2018_400_Tian2020_54_NodeName(Index);
    NodeNetIndex = Raw.DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork(Index);
    NodeCenter = Raw.ROICenter_Schafer_Tian_454(Index,:);
    YeoIndex = unique(NodeNetIndex);
    NetColorMap = Raw.YeoSCNetwork_ColorMap(YeoIndex,:);
    NetLabel = Raw.YeoSCNetwork_Label(YeoIndex,:);
    save([OutDir,filesep,Prefix,'ANCOVA_Dot',Thresholds{iThre}],'NodeLabel','NodeCenter','NodeNetIndex','NetColorMap','NetLabel');
    NetworkMatrix = EdgeMatrix(Index,Index);
    save([OutDir,filesep,'Ancova_Tresholded_Dot',Thresholds{iThre},'.mat'],'NetworkMatrix');
end



%% ANCOVA - Network connectivity - Combat 
% High order averaged FC matrix
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_HighOrderNet/FC_Matrix';
OutDir = ['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/FC_HighOrderNet'];
mkdir(OutDir);

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Dir = SubInfo.Dir(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
Sex = SubInfo.Sex(Index);
Site = SubInfo.Site(Index);

nNet = 8;
Mask = tril(ones(nNet,nNet));
% Mask = tril(ones(nNet,nNet),-1);
MaskIndex = find(Mask);

AllCov_comBat = [Dx,Age,Sex];
for iSub = 1:length(ID)
    load([InDir,filesep,'NetworkMatrix_',ID{iSub},'.mat']);
    DataOneDim = reshape(NetworkMatrix(MaskIndex),1,[]);
    FCAll(:,iSub) = DataOneDim;
end

Volume_Combat = combat(FCAll, Site', AllCov_comBat, 1);

for iSub = 1:size(Volume_Combat,2)
    CorrMatrix = zeros(nNet,nNet);
    CorrMatrix(MaskIndex) = Volume_Combat(:,iSub);
    CorrMatrix = CorrMatrix + tril(CorrMatrix,-1)';
    save([OutDir,filesep,'NetworkMatrix_',ID{iSub},'.mat'],'CorrMatrix');
    disp(num2str(iSub))
end



%% ANCOVA - Network connectivity - statistic
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat');
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/FC_HighOrderNet';
OutDir= '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_HighOrderNet';
mkdir(OutDir);

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Sex = SubInfo.Sex(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
MotionFD = SubInfo.MotionFD(Index);

ASDIndex = find(Dx==1); NCIndex = find(Dx==2); SCZIndex = find(Dx==3);

ASDList = cellfun(@(ID) [InDir,filesep,'NetworkMatrix_',ID,'.mat'],ID(ASDIndex),'UniformOutput',false);
NCList = cellfun(@(ID) [InDir,filesep,'NetworkMatrix_',ID,'.mat'],ID(NCIndex),'UniformOutput',false);
SCZList = cellfun(@(ID) [InDir,filesep,'NetworkMatrix_',ID,'.mat'],ID(SCZIndex),'UniformOutput',false);

CovariateMatrix = {[Age(ASDIndex),Sex(ASDIndex),MotionFD(ASDIndex)];...
    [Age(NCIndex),Sex(NCIndex),MotionFD(NCIndex)];...
    [Age(SCZIndex),Sex(SCZIndex),MotionFD(SCZIndex)]};

y_ANCOVA1_Image({ASDList,NCList,SCZList},[OutDir,filesep,'FC_Ancova'],'',{},CovariateMatrix);



%% ANCOVA - Functional gradient - generate gradient values
% clc;clear;
% load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
% InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20211119_ICVCovariate/NetworkAnalysis_FCMatrix_UseGUI';
% 
% SubList = cellfun(@(ID) dir([InDir,filesep,'NetworkMatrix_',ID,'.mat']),SubInfo.ID,'UniformOutput',false);
% SubList = cellfun(@(File) [File.folder,filesep,File.name],SubList,'UniformOutput',false);
% AllData = {};
% for iSub = 1:length(SubList)
% 	load(SubList{iSub});
%     if length(find(isnan(NetworkMatrix)))>10
%         disp(iSub); %sub-ASDQingHongCao20201111
%         AllData{iSub} = NetworkMatrix0;
%     else
%         NetworkMatrix(find(isnan(NetworkMatrix)))=0;
%         NetworkMatrix(find(isinf(NetworkMatrix)))=1;
%         AllData{iSub} = NetworkMatrix;
%     end
%     NetworkMatrix0 = NetworkMatrix;
% end
% 
% Gp = GradientMaps('kernel','normalized angle','approach','dm','alignment','pa','n_components',6);
% Gp = Gp.fit(AllData);



%% ANCOVA - Functional gradient - Combat and statistic
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat');
Exist = load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20211119.mat');
[C,IA,IB] = intersect(Exist.SubInfo.ID,SubInfo.ID);
ROIValue.Gradient.GradValue = Exist.ROIValue.Gradient.Gp.aligned(IA)';

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Dx = SubInfo.Dx(Index);
Site = SubInfo.Site(Index);
Age = SubInfo.Age(Index);
Sex = SubInfo.Sex(Index);
MotionFD = SubInfo.MotionFD(Index);

% Sex(find(Sex==0)) = -1;
SiteCov=[];
SiteIndex = unique(Site);
for i=1:length(SiteIndex)-1
    SiteCov(:,i) = Site==SiteIndex(i);
end

DxCov=[];
DxIndex = [1,2,3];
for i=1:length(DxIndex)-1
    DxCov(:,i) = Dx==DxIndex(i);
end

AllCov_ComBat = [Dx,Age,Sex];
for iGrad = 1:size(ROIValue.Gradient.GradValue{1},2)
    for iSub = 1:length(ID)
        GradAll(iGrad,:,iSub) = ROIValue.Gradient.GradValue{iSub}(:,iGrad);
    end
    Grad_Combat(iGrad,:,:) = combat(squeeze(GradAll(iGrad,:,:)), Site', AllCov_ComBat, 1);
end

for iGrad = 1:size(Grad_Combat,1)
    AllCov = [DxCov,ones(length(Dx),1),Age,Sex,MotionFD];
    Contrast=zeros(1,size(AllCov,2));
    Contrast(1:size(DxCov,2))=1;
    
    % Grad value for each ROI
    for iROI = 1:size(Grad_Combat,2)
        [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(squeeze(Grad_Combat(iGrad,iROI,:)),AllCov,Contrast,'F');
        F_Grad_ANCOVA(iGrad,iROI) = TF_ForContrast;
        P_Grad_ANCOVA(iGrad,iROI) = 1-fcdf(TF_ForContrast,size(DxCov,2),size(DxCov,1)-size(AllCov,2));
    end
    P_Grad_ANCOVA_Corrected(iGrad,:) =  mafdr(P_Grad_ANCOVA(iGrad,:),'BHFDR','true');
    
    % Grad range
    GradRange(:,iGrad) = (max(squeeze(Grad_Combat(iGrad,:,:)))-min(squeeze(Grad_Combat(iGrad,:,:))))';
    Index = find(GradRange(:,iGrad)<mean(GradRange(:,iGrad))+3*std(GradRange(:,iGrad)));
    [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(GradRange(:,iGrad),AllCov,Contrast,'F');
    F_GradRange_ANCOVA(iGrad) = TF_ForContrast;
    P_GradRange_ANCOVA(iGrad) = 1-fcdf(TF_ForContrast,size(DxCov,2),size(DxCov,1)-size(AllCov,2));
end

% Calculate overall dispersion
Dispersion = zeros(length(Dx),1);
for iSub = 1:size(Grad_Combat,3)
    Centriod = mean(Grad_Combat(1:2,:,iSub)');
    Distance = pdist2(Grad_Combat(1:2,:,iSub)',Centriod);
    Dispersion(iSub) = sum(Distance.^2);
end
Index = find(abs(Dispersion)<mean(Dispersion)+3*std(Dispersion));
[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(Dispersion(Index),AllCov(Index,:),Contrast,'F');
F_GradDispersion_ANCOVA = TF_ForContrast;
P_GradDispersion_ANCOVA = 1-fcdf(TF_ForContrast,size(DxCov,2),size(DxCov,1)-size(AllCov,2));

% Calculate dispersion in each network
Template = load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/DPABISurf_Schaefer2018_400_Tian2020_54_Info.mat');
Dispersion_All = zeros(length(Dx),8);
for iNet = 1:8
    NetIndex = find(Template.DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork==iNet);
    for iSub = 1:size(Grad_Combat,3)
        Centriod = mean(Grad_Combat(1:2,NetIndex,iSub)');
        Distance = pdist2(Grad_Combat(1:2,NetIndex,iSub)',Centriod);
        Dispersion_All(iSub,iNet) = sum(Distance.^2);
    end
    [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(Dispersion_All(:,iNet),AllCov,Contrast,'F');
    F_GradDispersionAll_ANCOVA(iNet) = TF_ForContrast;
    P_GradDispersionAll_ANCOVA(iNet) = 1-fcdf(TF_ForContrast,size(DxCov,2),size(DxCov,1)-size(AllCov,2));
end
P_GradDispersionAll_ANCOVA_Corrected =  mafdr(P_GradDispersionAll_ANCOVA,'BHFDR','true');

ROIValue.Gradient.GradValueCombat = Grad_Combat;
ROIValue.Gradient.Ancova.F_Grad_ANCOVA = F_Grad_ANCOVA;
ROIValue.Gradient.Ancova.P_Grad_ANCOVA = P_Grad_ANCOVA;
ROIValue.Gradient.Ancova.P_Grad_ANCOVA_Corrected = P_Grad_ANCOVA_Corrected;
ROIValue.Gradient.Ancova.F_GradRange_ANCOVA = F_GradRange_ANCOVA;
ROIValue.Gradient.Ancova.P_GradRange_ANCOVA = P_GradRange_ANCOVA;
ROIValue.Gradient.Ancova.F_GradDispersion_ANCOVA = F_GradDispersion_ANCOVA;
ROIValue.Gradient.Ancova.P_GradDispersion_ANCOVA = P_GradDispersion_ANCOVA;
ROIValue.Gradient.Ancova.F_GradDispersionAll_ANCOVA = F_GradDispersionAll_ANCOVA;
ROIValue.Gradient.Ancova.P_GradDispersionAll_ANCOVA = P_GradDispersionAll_ANCOVA;
ROIValue.Gradient.Ancova.P_GradDispersionAll_ANCOVA_Corrected = P_GradDispersionAll_ANCOVA_Corrected;
ROIValue.Gradient.Dispersion = Dispersion;
ROIValue.Gradient.Dispersion_All = Dispersion_All;
ROIValue.Gradient.GradRange = GradRange;

save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','SubInfo','ROIValue');



%% ANCOVA - Functional gradient - (Write gradient values on brains)
%%% Write the gradient values on the brains %%%
OutDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Gradient';
mkdir(OutDir);
[LeftTemplate, ~, ~, ~] = y_ReadAll('/mnt/Data4/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_lh_Schaefer2018_400Parcels_7Networks_order.label.gii');
[RightTemplate, ~, ~, ~] = y_ReadAll('/mnt/Data4/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_rh_Schaefer2018_400Parcels_7Networks_order.label.gii');
[~, ~, ~, HeaderL] = y_ReadAll('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/SeedFC/Corrected/FC_L_ROI2.gii');
[~, ~, ~, HeaderR] = y_ReadAll('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/SeedFC/Corrected/FC_R_ROI1.gii');
[VBrain,HeaderV] = y_Read('/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/Templates/Tian2020_Subcortex_Atlas/Tian_Subcortex_S4_3T_2009cAsym.nii');

%%% Write Ancova results brain %%%
OutputDir = [OutDir,filesep,'Ancova'];
mkdir(OutputDir);
for iGrad = 1:2
    GradBrainL = zeros(length(LeftTemplate),1);
    GradBrainR = zeros(length(RightTemplate),1);
    GradBrainV = zeros(size(VBrain));
    for iROI = 1:200
        if P_Grad_ANCOVA_Corrected(iGrad,iROI)<0.05
            GradBrainL(find(LeftTemplate==iROI)) = F_Grad_ANCOVA(iGrad,iROI);
        end
    end
    for iROI = 201:400
        if P_Grad_ANCOVA_Corrected(iGrad,iROI)<0.05
            GradBrainR(find(RightTemplate==iROI-200)) = F_Grad_ANCOVA(iGrad,iROI);
        end
    end
    for iROI = 401:454
        if P_Grad_ANCOVA_Corrected(iGrad,iROI)<0.05
            GradBrainV(find(VBrain==iROI-400)) = F_Grad_ANCOVA(iGrad,iROI);
        end
    end
    y_Write(GradBrainL,HeaderL,[OutputDir,filesep,'Ancova_Gradient',num2str(iGrad),'_Left']);
    y_Write(GradBrainR,HeaderR,[OutputDir,filesep,'Ancova_Gradient',num2str(iGrad),'_Right']);
    y_Write(GradBrainV,HeaderV,[OutputDir,filesep,'Ancova_Gradient',num2str(iGrad),'_Volume']);
end 

%%% Write brain areas showed sig diff of gradient, color it using Yeo7 %%%
GradBrainL = zeros(length(LeftTemplate),1);
GradBrainR = zeros(length(RightTemplate),1);
GradBrainV = zeros(size(VBrain));
for iROI = 1:200
    if P_Grad_ANCOVA_Corrected(1,iROI)<0.05 || P_Grad_ANCOVA_Corrected(2,iROI)<0.05
        GradBrainL(find(LeftTemplate==iROI)) = Template.DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork(iROI);
    end
end
for iROI = 201:400
    if P_Grad_ANCOVA_Corrected(1,iROI)<0.05 || P_Grad_ANCOVA_Corrected(2,iROI)<0.05
        GradBrainR(find(RightTemplate==iROI-200)) = Template.DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork(iROI);  
    end
end
for iROI = 401:454
    if P_Grad_ANCOVA_Corrected(1,iROI)<0.05 || P_Grad_ANCOVA_Corrected(2,iROI)<0.05
        GradBrainV(find(VBrain==iROI-400)) = Template.DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork(iROI);
    end
end
y_Write(GradBrainL,HeaderL,[OutputDir,filesep,'Ancova_Gradient12inYeo_Left']);
y_Write(GradBrainR,HeaderR,[OutputDir,filesep,'Ancova_Gradient12inYeo_Right']);
y_Write(GradBrainV,HeaderV,[OutputDir,filesep,'Ancova_Gradient12inYeo_Volume']);


%%% Write mean gradient brain %%%
OutputDir = [OutDir,filesep,'MeanGradient'];
mkdir(OutputDir);
DxList = {'ASD','NC','SCZ'};
ROIValue.Gradient.GradValueMean = zeros(2,454,3); % (grad,Dx,ROI)
for iDx = 1:3
    for iGrad = 1:2
        Index = find(Dx == iDx);
        GradMean = mean(squeeze(Grad_Combat(iGrad,:,Index)),2);
        GradBrainL = zeros(length(LeftTemplate),1);
        GradBrainR = zeros(length(RightTemplate),1);
        GradBrainV = zeros(size(VBrain));
        for iROI = 1:200
            GradBrainL(find(LeftTemplate==iROI)) = GradMean(iROI);
        end
        for iROI = 201:400
            GradBrainR(find(RightTemplate==iROI-200)) = GradMean(iROI);
        end
        for iROI = 401:454
            GradBrainV(find(VBrain==iROI-400)) = GradMean(iROI);
        end
        y_Write(GradBrainL,HeaderL,[OutputDir,filesep,'MeanMap_Gradient',num2str(iGrad),'_',DxList{iDx},'_Left']);
        y_Write(GradBrainR,HeaderR,[OutputDir,filesep,'MeanMap_Gradient',num2str(iGrad),'_',DxList{iDx},'_Right']);
        y_Write(GradBrainV,HeaderV,[OutputDir,filesep,'MeanMap_Gradient',num2str(iGrad),'_',DxList{iDx},'_Volume']);
        ROIValue.Gradient.GradValueMean(iGrad,:,iDx) = GradMean;
    end
end

save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','SubInfo','ROIValue');


% gradient_in_euclidean([squeeze(GradMean(1,1,:)),squeeze(GradMean(1,2,:))]);%ASD
% gradient_in_euclidean([squeeze(GradMean(2,1,:)),squeeze(GradMean(2,2,:))]);%NC
% gradient_in_euclidean([squeeze(GradMean(3,1,:)),squeeze(GradMean(3,2,:))]);%SCZ














%% %%%%%% Extract ROI Values %%%%%%

%% Extract subcortical volume values 
% Subcortical volume data has been extracted during ANCOVA analysis

%% Extract cortical anatomical metrics
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat');
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/Cortical_Structure';
DataType = {['AnatSurfLH',filesep,'Thickness'],...
    ['AnatSurfRH',filesep,'Thickness'],...
    ['AnatSurfRH',filesep,'Thickness']};
ROIDir = {'/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/Anat/Thickness_L_1_Mask.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/Anat/Thickness_R_1_Mask.gii',...
    '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/Anat/Thickness_R_2_Mask.gii'};

CorticalMetrics = zeros(length(SubInfo.ID),length(DataType));
for i = 1:length(DataType)
    SubList = cellfun(@(ID) dir([InDir,filesep,DataType{i},filesep,ID,'*']),SubInfo.ID,'UniformOutput',false);
    SubList = cellfun(@(File) [File.folder,filesep,File.name],SubList,'UniformOutput',false);
    Flags = cellfun(@(Dir) ~strcmp(Dir,filesep),SubList);
    Index = find(Flags);
    SubList = SubList(Index);
    [AllVolume, VoxelSize, FileList, Header] = y_ReadAll(SubList);
    [MaskData, VoxelSize1, FileList1, Header1] = y_ReadAll(ROIDir{i});
    MaskIndex = find(MaskData);
    CorticalMetrics(Index,i) = mean(AllVolume(MaskIndex,:),1)';
    disp(num2str(i))
end

ROIValue.CorticalMetrics_Sig = CorticalMetrics;
ROIValue.CorticalMetrics_SigFlag = DataType;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat','SubInfo','ROIValue');

%% Extract seed-based FC values
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat');
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/SeedFC';
ROIDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/SeedFC/Corrected';
ROIName = {         'FC_R_ROI1.gii',...
    'FC_L_ROI2.gii',                ...
    'FC_L_ROI3.gii',                ...
    'FC_L_ROI5.gii','FC_R_ROI5.gii',...
    'FC_L_ROI6.gii','FC_R_ROI6.gii',...
    'FC_L_ROI7.gii','FC_R_ROI7.gii',...
    'FC_L_ROI8.gii','FC_R_ROI8.gii'};

for i = 1:length(ROIName)
    temp = strsplit(ROIName{i},{'_','.'});
    ROIFlags(i) = str2num(temp{3}(4:end));
    if strcmp(temp{2}(1),'L')
        Hemi = 'FunSurfLH';
    elseif strcmp(temp{2}(1),'R')
        Hemi = 'FunSurfRH';
    end
    SubList = cellfun(@(ID) dir([InDir,filesep,Hemi,filesep,'ROI',num2str(ROIFlags(i)),filesep,ID,'*']),SubInfo.ID,'UniformOutput',false);
    SubList = cellfun(@(File) [File.folder,filesep,File.name],SubList,'UniformOutput',false);
    [AllVolume, VoxelSize, FileList, Header] = y_ReadAll(SubList);
    [MaskData, VoxelSize1, FileList1, Header1] = y_ReadAll([ROIDir,filesep,ROIName{i}]);

    MaskIndex = find(MaskData);
    SeedFC(:,i) = mean(AllVolume(MaskIndex,:),1)';
    disp(num2str(i))
end

ROIValue.SeedFC.FC = SeedFC;
ROIValue.SeedFC.ROIName = ROIName;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','SubInfo','ROIValue');

%% Extract connectome values
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/FC_Network';

% FDR005
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_Network/EdgeMatrix_Ancova_FDR_Dot005.mat');
SubList = cellfun(@(ID) dir([InDir,filesep,'NetworkMatrix_',ID,'.mat']),SubInfo.ID,'UniformOutput',false);
SubList = cellfun(@(File) [File.folder,filesep,File.name],SubList,'UniformOutput',false);

Index = find(triu(EdgeMatrix));
for iSub = 1:length(SubList)
	load(SubList{iSub});
    AllData(iSub,:) = CorrMatrix(Index);
    disp(num2str(iSub))
end

ROIValue.Connetcome.SigFC_FDR005 = AllData;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','SubInfo','ROIValue');

%% Extract network connectivity values
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_HighOrderNet/FC_Matrix';

% FDR05
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_HighOrderNet/EdgeMatrix_Ancova_FDR_Dot05.mat');
SubList = cellfun(@(ID) dir([InDir,filesep,'NetworkMatrix_',ID,'.mat']),SubInfo.ID,'UniformOutput',false);
SubList = cellfun(@(File) [File.folder,filesep,File.name],SubList,'UniformOutput',false);

Index = find(triu(EdgeMatrix,-1));
for iSub = 1:length(SubList)
	load(SubList{iSub});
    AllData(iSub,:) = NetworkMatrix(Index);
    disp(num2str(iSub))
end

ROIValue.NetworkFC.SigFC_FDR05 = AllData;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','SubInfo','ROIValue');

%% Extract ROI2ROI-FC values
% ROI2ROI-FC data has been extracted during ANCOVA analysis

%% Extract gradient values
% gradient data has been extracted during ANCOVA analysis

%% Extract graph theoretical matrics
% graph theoretical matrics data has been extracted during ANCOVA analysis

%% Extract regional metrics
% clc;clear;
% load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
% InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/Regional_Metrics';
% Metrics  = {'ALFF_FunSurfWC','fALFF_FunSurfWC','ReHo_FunSurfWCF','DegreeCentrality_FunSurfWCF'};
% Hemispheres = {'FunSurfLH','FunSurfRH','FunVolu'};
% HemiLabels = {'L','R','V'};
% HemiSurfix = {'.gii','.gii','.nii'};
% ROINames = {'ALFF','fALFF','ReHo','DC'};
% nROIs = {[4,4,2],[4,3,0],[1,2,0],[2,1,2]}; % number of clusters in left, right and subcortical brain
% 
% 
% for iMetric = 1:length(ROINames)
%     for iHemi = 1:length(Hemispheres)
%         SubList = cellfun(@(ID) dir([InDir,filesep,Hemispheres{iHemi},filesep,...
%             Metrics{iMetric},filesep,ID,'*']),SubInfo.ID);
%         SubList = cellfun(@(Dir,ID) [Dir,filesep,ID],{SubList(:).folder},{SubList(:).name},'UniformOutput',false);
%         [AllVolume, VoxelSize, FileList, Header] = y_ReadAll(SubList);
%         if iHemi <3
%             if nROIs{iMetric}(iHemi) >0
%                 for iROI = 1:nROIs{iMetric}(iHemi)
%                     [MaskData, VoxelSize1, FileList1, Header1] = y_ReadAll(['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/RegionalMetric',...
%                         filesep,ROINames{iMetric},'_',HemiLabels{iHemi},'_',num2str(iROI),HemiSurfix{iHemi}]);
%                     MaskIndex = find(MaskData);
%                     RegionalMetric(:,iMetric,iHemi,iROI) = mean(AllVolume(MaskIndex,:),1)';
%                     MetricLabel{:,iMetric,iHemi,iROI} = [ROINames{iMetric},'_',HemiLabels{iHemi},'_',num2str(iROI)];
%                     disp([ROINames{iMetric},'_',HemiLabels{iHemi},'_',num2str(iROI)])
%                 end
%             end
%         else
%             AllVolume = reshape(AllVolume,[],size(AllVolume,4));
%             if nROIs{iMetric}(iHemi) >0
%                 for iROI = 1:nROIs{iMetric}(iHemi)
%                     [MaskData, VoxelSize1, FileList1, Header1] = y_ReadAll(['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/RegionalMetric',...
%                         filesep,ROINames{iMetric},'_',HemiLabels{iHemi},'_',num2str(iROI),HemiSurfix{iHemi}]);
%                     MaskData = reshape(MaskData,[],1);
%                     MaskIndex = find(MaskData);
%                     RegionalMetric(:,iMetric,iHemi,iROI) = mean(AllVolume(MaskIndex,:),1)';
%                     MetricLabel{:,iMetric,iHemi,iROI} = [ROINames{iMetric},'_',HemiLabels{iHemi},'_',num2str(iROI)];
%                     disp([ROINames{iMetric},'_',HemiLabels{iHemi},'_',num2str(iROI)])
%                 end
%             end
%         end
%     end
% end
% 
% ROIValue.RegionalMetric.Value = RegionalMetric;
% ROIValue.RegionalMetric.Label = MetricLabel;
% save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','SubInfo','ROIValue');



%% %%%%%% Post-hoc analysis %%%%%%

%% Post-hoc - Anatomical metrics
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat')
Pairs = [1,2;2,3;1,3]; % ASD-NC, NC-SCZ, ASD-SCZ
AnatResults = [ROIValue.CorticalMetrics_Sig,ROIValue.SubcortVol_Sig];
Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 );

AnatResults = AnatResults(Index,:);
Sex = SubInfo.Sex(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
ICV = SubInfo.ICV(Index);

for iPair = 1:size(Pairs,1)
    GroupIndex1 = find(Dx ==Pairs(iPair,1));
    GroupIndex2 = find(Dx ==Pairs(iPair,2));
    
    DxCov = [ones(length(GroupIndex1),1);ones(length(GroupIndex2),1).*-1];
    AgeCov = [Age(GroupIndex1);Age(GroupIndex2)];
    SexCov = [Sex(GroupIndex1);Sex(GroupIndex2)];
    ICVCov = [ICV(GroupIndex1);ICV(GroupIndex2)];
    AnatResults1 = [AnatResults(GroupIndex1,:);AnatResults(GroupIndex2,:)];

    for iROI = 1:size(AnatResults1,2)
        if iROI < 4 %%%%%%%% cortical thickness
            ScalingFactor = 1/3;
        else
            ScalingFactor = 1;
        end
        ICVCov = ICVCov.^ScalingFactor;
        AllCov = [DxCov,ones(length(DxCov),1),AgeCov,SexCov,ICVCov];
        Contrast=zeros(1,size(AllCov,2));
        Contrast(1)=1;
        AnatValue = AnatResults1(:,iROI);
        [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(AnatValue,AllCov,Contrast,'T');
        PostHocT_Anat(iROI,iPair) = TF_ForContrast
        if TF_ForContrast>=0
            PostHocP_Anat(iROI,iPair) = 1-tcdf(TF_ForContrast,length(DxCov)-2)
        else
            PostHocP_Anat(iROI,iPair) = tcdf(TF_ForContrast,length(DxCov)-2)
        end
    end
end
 
ROIValue.PostHoc.T = PostHocT_Anat;
ROIValue.PostHoc.P = PostHocP_Anat;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat','ROIValue','SubInfo');



%% Post-hoc - Seed-based FC
% average significant FCs in left brain and right brain 
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

Pairs = [1,2;2,3;1,3]; % ASD-NC, NC-SCZ, ASD-SCZ
Sex = SubInfo.Sex(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
Motion = SubInfo.MotionFD(Index);

nSeed = 7;

FCResults = zeros(length(Dx),nSeed);
ROIIndex = {[1];...
    [2];...
    [3];...
    [4,5];...
    [6,7];...
    [8,9];...
    [10,11]};

for i = 1:nSeed
    FCResults(:,i) = mean(ROIValue.SeedFC.FC(:,ROIIndex{i}),2);
end
    
for iPair = 1:size(Pairs,1)
    GroupIndex1 = find(Dx ==Pairs(iPair,1));
    GroupIndex2 = find(Dx ==Pairs(iPair,2));
    
    DxCov = [ones(length(GroupIndex1),1);ones(length(GroupIndex2),1).*-1];
    AgeCov = [Age(GroupIndex1);Age(GroupIndex2)];
    SexCov = [Sex(GroupIndex1);Sex(GroupIndex2)];
    MotionCov = [Motion(GroupIndex1);Motion(GroupIndex2)];
    FCResults1 = [FCResults(GroupIndex1,:);FCResults(GroupIndex2,:)];    

    AllCov = [DxCov,ones(length(DxCov),1),AgeCov,SexCov,MotionCov];
    Contrast=zeros(1,size(AllCov,2));
    Contrast(1)=1;
    
    for iROI = 1:size(FCResults1,2)
        FCValue = FCResults1(:,iROI);
        [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(FCValue,AllCov,Contrast,'T');
        PostHoc_T_SeedFC(iROI,iPair) = TF_ForContrast;
        if TF_ForContrast>=0
            PostHoc_P_SeedFC(iROI,iPair) = 1-tcdf(TF_ForContrast,length(DxCov)-2);
        else
            PostHoc_P_SeedFC(iROI,iPair) = tcdf(TF_ForContrast,length(DxCov)-2);
        end
    end
end

ROIValue.SeedFC.PostHoc.T = PostHoc_T_SeedFC;
ROIValue.SeedFC.PostHoc.P = PostHoc_P_SeedFC;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','ROIValue','SubInfo');



%% Post-hoc - ROI-ROI FC
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
Pairs = [1,2;2,3;1,3]; % ASD-NC, NC-SCZ, ASD-SCZ

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

Sex = SubInfo.Sex(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
Motion = SubInfo.MotionFD(Index);

FCResults = ROIValue.ROI2ROIFC.FC_Sig(Index,:);

for iPair = 1:size(Pairs,1)
    GroupIndex1 = find(Dx ==Pairs(iPair,1));
    GroupIndex2 = find(Dx ==Pairs(iPair,2));
    
    DxCov = [ones(length(GroupIndex1),1);ones(length(GroupIndex2),1).*-1];
    AgeCov = [Age(GroupIndex1);Age(GroupIndex2)];
    SexCov = [Sex(GroupIndex1);Sex(GroupIndex2)];
    MotionCov = [Motion(GroupIndex1);Motion(GroupIndex2)];
    FCResults1 = [FCResults(GroupIndex1,:);FCResults(GroupIndex2,:)];    
    
    AllCov = [DxCov,ones(length(DxCov),1),AgeCov,SexCov,MotionCov];
    Contrast=zeros(1,size(AllCov,2));
    Contrast(1)=1;
    
    for iROI = 1:size(FCResults1,2)
        FCValue = FCResults1(:,iROI);
        [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(FCValue,AllCov,Contrast,'T');
        PostHoc_T_ROI2ROIFC(iROI,iPair) = TF_ForContrast;
        if TF_ForContrast>=0
            PostHoc_P_ROI2ROIFC(iROI,iPair) = 1-tcdf(TF_ForContrast,length(DxCov)-2);
        else
            PostHoc_P_ROI2ROIFC(iROI,iPair) = tcdf(TF_ForContrast,length(DxCov)-2);
        end
    end
end

ROIValue.ROI2ROIFC.PostHoc.T = PostHoc_T_ROI2ROIFC;
ROIValue.ROI2ROIFC.PostHoc.P = PostHoc_P_ROI2ROIFC;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','ROIValue','SubInfo');



%% Post-hoc - Connectome
% 454*454
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/FC_Network';
OutDir= '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_Network';
mkdir(OutDir);

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Sex = SubInfo.Sex(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
MotionFD = SubInfo.MotionFD(Index);

nROI = 454;
AllMatrix = zeros(nROI,nROI,length(ID));
for iSub = 1:length(ID)
    load([InDir,filesep,'NetworkMatrix_',ID{iSub},'.mat']);
    AllMatrix(:,:,iSub) = CorrMatrix;
end

Thresholds = {'001','005'};
for iThre = 1:length(Thresholds)
    load(['/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_Network/EdgeMatrix_Ancova_FDR_Dot',Thresholds{iThre},'.mat']);
    Mask = EdgeMatrix;
    
    PostHoc_T_NetworkFC = zeros(nROI,nROI,3);
    PostHoc_P_NetworkFC = zeros(nROI,nROI,3);
    
    Pairs = [1,2;2,3;1,3]; % ASD-NC, NC-SCZ, ASD-SCZ
    for iPair = 1:size(Pairs,1)
        GroupIndex1 = find(Dx ==Pairs(iPair,1));
        GroupIndex2 = find(Dx ==Pairs(iPair,2));
        
        DxCov = [ones(length(GroupIndex1),1);ones(length(GroupIndex2),1).*-1];
        AgeCov = [Age(GroupIndex1);Age(GroupIndex2)];
        SexCov = [Sex(GroupIndex1);Sex(GroupIndex2)];
        MotionCov = [MotionFD(GroupIndex1);MotionFD(GroupIndex2)];
        
        AllCov = [DxCov,ones(length(DxCov),1),AgeCov,SexCov,MotionCov];
        Contrast=zeros(1,size(AllCov,2));
        Contrast(1)=1;
        df(iPair) = length(AgeCov)-2;
        
        for iROI = 1:size(AllMatrix,1)
            for jROI = 1:size(AllMatrix,2)
                if Mask(iROI,jROI)~=0
                    FCValue = [squeeze(AllMatrix(iROI,jROI,GroupIndex1));squeeze(AllMatrix(iROI,jROI,GroupIndex2))];
                    [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(FCValue,AllCov,Contrast,'T');
                    PostHoc_T_NetworkFC(iROI,jROI,iPair) = TF_ForContrast;
                    if TF_ForContrast>=0
                        PostHoc_P_NetworkFC(iROI,jROI,iPair) = 1-tcdf(TF_ForContrast,length(DxCov)-2);
                    else
                        PostHoc_P_NetworkFC(iROI,jROI,iPair) = tcdf(TF_ForContrast,length(DxCov)-2);
                    end
                end
                disp([num2str(iROI),filesep,num2str(jROI)])
            end
        end
    end
    PostHoc_T_NetworkFC_ASD_NC = PostHoc_T_NetworkFC(:,:,1);
    PostHoc_T_NetworkFC_SCZ_NC = -PostHoc_T_NetworkFC(:,:,2); %%%%minus
    PostHoc_T_NetworkFC_ASD_SCZ = PostHoc_T_NetworkFC(:,:,3);
    PostHoc_P_NetworkFC_ASD_NC = PostHoc_P_NetworkFC(:,:,1);
    PostHoc_P_NetworkFC_SCZ_NC = PostHoc_P_NetworkFC(:,:,2);
    PostHoc_P_NetworkFC_ASD_SCZ = PostHoc_P_NetworkFC(:,:,3);
    
    StatOpt.TestFlag = 'T';
    StatOpt.Df = df(1);
    save([OutDir,filesep,'TMatrix_PostHoc_FDR_',Thresholds{iThre},'_ASD_NC.mat'],...
        'PostHoc_T_NetworkFC_ASD_NC','PostHoc_P_NetworkFC_ASD_NC','StatOpt');
    StatOpt.Df = df(2);
    save([OutDir,filesep,'TMatrix_PostHoc_FDR_',Thresholds{iThre},'_SCZ_NC.mat'],...
        'PostHoc_T_NetworkFC_SCZ_NC','PostHoc_P_NetworkFC_SCZ_NC','StatOpt');
    StatOpt.Df = df(3);
    save([OutDir,filesep,'TMatrix_PostHoc_FDR_',Thresholds{iThre},'_ASD_SCZ.mat'],...
        'PostHoc_T_NetworkFC_ASD_SCZ','PostHoc_P_NetworkFC_ASD_SCZ','StatOpt');
end



%% Post-hoc - Connectome - (Revision template file for DPABINet_VIEW)
clc;clear;
Raw = load('/mnt/Data4/RfMRILab/Lubin/Project/RMD_Cross_Diagnosis/StatisticalResults/NetworkAnalysis_20220803/Short_DPABISurf_Schaefer2018_400_Tian2020_54_Info.mat');
OutDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_Network';
Prefix = 'Template_forDPABINetVIEW_';

Pairs = {'ASD_NC','SCZ_NC','ASD_SCZ'}; 
Correction = '005';
for iPair = 1:length(Pairs)
    Stat = load([OutDir,filesep,'TMatrix_PostHoc_FDR_',Correction,'_',Pairs{iPair},'.mat']);
    SigMatrix = Stat.(['PostHoc_P_NetworkFC_',Pairs{iPair}])<0.0167./2 & Stat.(['PostHoc_P_NetworkFC_',Pairs{iPair}])>0; %%%% Critical note! 0.0167./2 for two tail!!%%%% %%%
    Index = find(any(SigMatrix));
    NodeLabel = Raw.DPABISurf_Schaefer2018_400_Tian2020_54_NodeName(Index);
    NodeNetIndex = Raw.DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork(Index);
    NodeCenter = Raw.ROICenter_Schafer_Tian_454(Index,:);
    YeoIndex = unique(NodeNetIndex);
    NetColorMap = Raw.YeoSCNetwork_ColorMap(YeoIndex,:);
    NetLabel = Raw.YeoSCNetwork_Label(YeoIndex,:);
    save([OutDir,filesep,Prefix,'PostHoc_FDR_',Correction,'_',Pairs{iPair}],'NodeLabel','NodeCenter','NodeNetIndex','NetColorMap','NetLabel');
    Stat.(['PostHoc_T_NetworkFC_',Pairs{iPair}])(find(SigMatrix==0)) = 0;
    AfterEffect_T = Stat.(['PostHoc_T_NetworkFC_',Pairs{iPair}])(Index,Index);
    StatOpt = Stat.StatOpt;
    save([OutDir,filesep,'TMatrix_PostHoc_FDR_',Correction,'_',Pairs{iPair},'_Tresholded'],'AfterEffect_T');
end



%% Post-hoc - Network connectivity 
% 8*8
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
InDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Combat_Metrics/FC_HighOrderNet';
OutDir= '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/FC_HighOrderNet';
mkdir(OutDir);

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Sex = SubInfo.Sex(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
MotionFD = SubInfo.MotionFD(Index);

nROI = 8;
AllMatrix = zeros(nROI,nROI,length(ID));
for iSub = 1:length(ID)
    load([InDir,filesep,'NetworkMatrix_',ID{iSub},'.mat']);
    AllMatrix(:,:,iSub) = CorrMatrix;
end

Thresholds = {'05'};
for iThre = 1:length(Thresholds)
    load([OutDir,filesep,'EdgeMatrix_Ancova_FDR_Dot',Thresholds{iThre},'.mat']);
    Mask = EdgeMatrix;
    
    PostHoc_T_NetworkFC = zeros(nROI,nROI,3);
    PostHoc_P_NetworkFC = zeros(nROI,nROI,3);
    
    Pairs = [1,2;2,3;1,3]; % ASD-NC, NC-SCZ, ASD-SCZ
    for iPair = 1:size(Pairs,1)
        GroupIndex1 = find(Dx ==Pairs(iPair,1));
        GroupIndex2 = find(Dx ==Pairs(iPair,2));
        
        DxCov = [ones(length(GroupIndex1),1);ones(length(GroupIndex2),1).*-1];
        AgeCov = [Age(GroupIndex1);Age(GroupIndex2)];
        SexCov = [Sex(GroupIndex1);Sex(GroupIndex2)];
        MotionCov = [MotionFD(GroupIndex1);MotionFD(GroupIndex2)];
        
        AllCov = [DxCov,ones(length(DxCov),1),AgeCov,SexCov,MotionCov];
        Contrast=zeros(1,size(AllCov,2));
        Contrast(1)=1;
        df(iPair) = length(AgeCov)-2;
        
        for iROI = 1:size(AllMatrix,1)
            for jROI = 1:size(AllMatrix,2)
                if Mask(iROI,jROI)~=0
                    FCValue = [squeeze(AllMatrix(iROI,jROI,GroupIndex1));squeeze(AllMatrix(iROI,jROI,GroupIndex2))];
                    [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(FCValue,AllCov,Contrast,'T');
                    PostHoc_T_NetworkFC(iROI,jROI,iPair) = TF_ForContrast;
                    if TF_ForContrast>=0
                        PostHoc_P_NetworkFC(iROI,jROI,iPair) = 1-tcdf(TF_ForContrast,length(DxCov)-2);
                    else
                        PostHoc_P_NetworkFC(iROI,jROI,iPair) = tcdf(TF_ForContrast,length(DxCov)-2);
                    end
                end
                disp([num2str(iROI),filesep,num2str(jROI)])
            end
        end
    end
    PostHoc_T_NetworkFC_ASD_NC = PostHoc_T_NetworkFC(:,:,1);
    PostHoc_T_NetworkFC_SCZ_NC = -PostHoc_T_NetworkFC(:,:,2); %%%%minus
    PostHoc_T_NetworkFC_ASD_SCZ = PostHoc_T_NetworkFC(:,:,3);
    PostHoc_P_NetworkFC_ASD_NC = PostHoc_P_NetworkFC(:,:,1);
    PostHoc_P_NetworkFC_SCZ_NC = PostHoc_P_NetworkFC(:,:,2);
    PostHoc_P_NetworkFC_ASD_SCZ = PostHoc_P_NetworkFC(:,:,3);
    
    StatOpt.TestFlag = 'T';
    StatOpt.Df = df(1);
    save([OutDir,filesep,'TMatrix_PostHoc_FDR_',Thresholds{iThre},'_ASD_NC.mat'],...
        'PostHoc_T_NetworkFC_ASD_NC','PostHoc_P_NetworkFC_ASD_NC','StatOpt');
    StatOpt.Df = df(2);
    save([OutDir,filesep,'TMatrix_PostHoc_FDR_',Thresholds{iThre},'_SCZ_NC.mat'],...
        'PostHoc_T_NetworkFC_SCZ_NC','PostHoc_P_NetworkFC_SCZ_NC','StatOpt');
    StatOpt.Df = df(3);
    save([OutDir,filesep,'TMatrix_PostHoc_FDR_',Thresholds{iThre},'_ASD_SCZ.mat'],...
        'PostHoc_T_NetworkFC_ASD_SCZ','PostHoc_P_NetworkFC_ASD_SCZ','StatOpt');
end



%% Post-hoc - Regional metrics
% clc;clear;
% load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
% Pairs = [1,2;2,3;1,3]; % ASD-NC, NC-SCZ, ASD-SCZ
% 
% % Site = SubInfo.Site;
% Sex = SubInfo.Sex;
% Dx = SubInfo.Dx;
% Age = SubInfo.Age;
% MotionFD = SubInfo.MotionFD;
% 
% ExtractedValue = ROIValue.RegionalMetric.Value;
% 
% for iPair = 1:size(Pairs,1)
%     DxCov = [ones(length(find(Dx==Pairs(iPair,1))),1);ones(length(find(Dx==Pairs(iPair,2))),1).*-1];
%     AgeCov = [Age(find(Dx==Pairs(iPair,1)),:);Age(find(Dx==Pairs(iPair,2)),:)];
%     SexCov = [Sex(find(Dx==Pairs(iPair,1)),:);Sex(find(Dx==Pairs(iPair,2)),:)];
% %     SiteNew = [Site(find(Dx==Pairs(iPair,1)),:);Site(find(Dx==Pairs(iPair,2)),:)];
% %     SiteCov=[];
% %     SiteIndex = unique(SiteNew);
% %     for i=1:length(SiteIndex)-1
% %         SiteCov(:,i) = SiteNew==SiteIndex(i);
% %     end
%     MotionCov = [MotionFD(find(Dx==Pairs(iPair,1)),:);MotionFD(find(Dx==Pairs(iPair,2)),:)];
%     AllCov = [DxCov,ones(length(DxCov),1),AgeCov,SexCov,MotionCov]; %%%%%%     no need for sitecov after combat
%     Contrast=zeros(1,size(AllCov,2));
%     Contrast(1)=1;
%     
%     for iMetric = 1:size(ExtractedValue,2)
%         for iHemi = 1:size(ExtractedValue,3)
%             for iCluster = 1:size(ExtractedValue,4)
%                 MetricValue = [ExtractedValue(find(Dx==Pairs(iPair,1)),iMetric,iHemi,iCluster);ExtractedValue(find(Dx==Pairs(iPair,2)),iMetric,iHemi,iCluster)];
%                 [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(MetricValue,AllCov,Contrast,'T');
%                 AfterEffect_T_RegionalMetric(iMetric,iHemi,iCluster,iPair) = TF_ForContrast
%                 if TF_ForContrast>=0
%                     AfterEffect_P_RegionalMetric(iMetric,iHemi,iCluster,iPair) = 1-tcdf(TF_ForContrast,length(DxCov)-2)
%                 else
%                     AfterEffect_P_RegionalMetric(iMetric,iHemi,iCluster,iPair) = tcdf(TF_ForContrast,length(DxCov)-2)
%                 end
%             end
%         end
%     end
% end
% 
% save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/RegionalMetric/Multiple_Compare_RegionalMetric.mat',...
%     'AfterEffect_T_RegionalMetric','AfterEffect_P_RegionalMetric');



%% Post-hoc - functional gradient
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat')
OutDir= '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Gradient/PostHoc';
mkdir(OutDir);

Index = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2);

ID = SubInfo.ID(Index);
Sex = SubInfo.Sex(Index);
Dx = SubInfo.Dx(Index);
Age = SubInfo.Age(Index);
MotionFD = SubInfo.MotionFD(Index);
GradValueCombat = ROIValue.Gradient.GradValueCombat(:,:,Index);


nROI = 454;
Pairs = [1,2;2,3;1,3]; % ASD-NC, NC-SCZ, ASD-SCZ
PostHoc_T_Gradient = zeros(2,size(GradValueCombat,2),size(Pairs,1)); %% only do gradient 1-2
PostHoc_P_Gradient = ones(2,size(GradValueCombat,2),size(Pairs,1));
PostHoc_T_GradRange = zeros(2,size(Pairs,1)); %% only do gradient 1-2
PostHoc_P_GradRange = ones(2,size(Pairs,1)); 
PostHoc_T_Dispersion = zeros(size(Pairs,1),1); 
PostHoc_P_Dispersion = ones(size(Pairs,1),1); 
for iPair = 1:size(Pairs,1)
    GroupIndex1 = find(Dx ==Pairs(iPair,1));
    GroupIndex2 = find(Dx ==Pairs(iPair,2));
    
    DxCov = [ones(length(GroupIndex1),1);ones(length(GroupIndex2),1).*-1];
    AgeCov = [Age(GroupIndex1);Age(GroupIndex2)];
    SexCov = [Sex(GroupIndex1);Sex(GroupIndex2)];
    MotionCov = [MotionFD(GroupIndex1);MotionFD(GroupIndex2)];
    
    AllCov = [DxCov,ones(length(DxCov),1),AgeCov,SexCov,MotionCov];
    Contrast=zeros(1,size(AllCov,2));
    Contrast(1)=1;
    df(iPair) = length(AgeCov)-2;
    
    for iGrad = 1:2
        for iROI = 1:nROI
            if ROIValue.Gradient.Ancova.P_Grad_ANCOVA_Corrected(iGrad,iROI)<0.05
                GradientValue = [squeeze(GradValueCombat(iGrad,iROI,GroupIndex1));squeeze(GradValueCombat(iGrad,iROI,GroupIndex2))];
                [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(GradientValue,AllCov,Contrast,'T');
                PostHoc_T_Gradient(iGrad,iROI,iPair) = TF_ForContrast;
                if TF_ForContrast>=0
                    PostHoc_P_Gradient(iGrad,iROI,iPair) = 1-tcdf(TF_ForContrast,length(DxCov)-2);
                else
                    PostHoc_P_Gradient(iGrad,iROI,iPair) = tcdf(TF_ForContrast,length(DxCov)-2);
                end
            end
        end
        
        % Grad Range
        Range = ROIValue.Gradient.GradRange([GroupIndex1;GroupIndex2],iGrad);
        [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(Range,AllCov,Contrast,'T');
        PostHoc_T_GradRange(iGrad,iPair) = TF_ForContrast;
        if TF_ForContrast>=0
            PostHoc_P_GradRange(iGrad,iPair) = 1-tcdf(TF_ForContrast,length(DxCov)-2);
        else
            PostHoc_P_GradRange(iGrad,iPair) = tcdf(TF_ForContrast,length(DxCov)-2);
        end
    end
    
    % Grad Dispersion
    Dispersion = ROIValue.Gradient.Dispersion([GroupIndex1;GroupIndex2]);
    [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(Dispersion,AllCov,Contrast,'T');
    PostHoc_T_Dispersion(iPair) = TF_ForContrast;
    if TF_ForContrast>=0
        PostHoc_P_Dispersion(iPair) = 1-tcdf(TF_ForContrast,length(DxCov)-2);
    else
        PostHoc_P_Dispersion(iPair) = tcdf(TF_ForContrast,length(DxCov)-2);
    end
end

[mean(ROIValue.Gradient.GradRange(find(Dx==1),1)),mean(ROIValue.Gradient.GradRange(find(Dx==2),1)),mean(ROIValue.Gradient.GradRange(find(Dx==3),1))]
[mean(ROIValue.Gradient.GradRange(find(Dx==1),2)),mean(ROIValue.Gradient.GradRange(find(Dx==2),2)),mean(ROIValue.Gradient.GradRange(find(Dx==3),2))]
[mean(ROIValue.Gradient.Dispersion(find(Dx==1))),mean(ROIValue.Gradient.Dispersion(find(Dx==2))),mean(ROIValue.Gradient.Dispersion(find(Dx==3)))]


ROIValue.Gradient.PostHoc.T_Gradient = PostHoc_T_Gradient;
ROIValue.Gradient.PostHoc.P_Gradient = PostHoc_P_Gradient;
ROIValue.Gradient.PostHoc.T_GradRange = PostHoc_T_GradRange;
ROIValue.Gradient.PostHoc.P_GradRange = PostHoc_P_GradRange;
ROIValue.Gradient.PostHoc.T_Dispersion = PostHoc_T_Dispersion;
ROIValue.Gradient.PostHoc.P_Dispersion = PostHoc_P_Dispersion;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','ROIValue','SubInfo');


%%% Write Post hoc results brain
[LeftTemplate, ~, ~, ~] = y_ReadAll('/mnt/Data4/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_lh_Schaefer2018_400Parcels_7Networks_order.label.gii');
[RightTemplate, ~, ~, ~] = y_ReadAll('/mnt/Data4/RfMRILab/Lubin/Software/DPABI_V6.0_210501/DPABISurf/SurfTemplates/fsaverage5_rh_Schaefer2018_400Parcels_7Networks_order.label.gii');
[~, ~, ~, HeaderL] = y_ReadAll('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/SeedFC/Corrected/FC_L_ROI2.gii');
[~, ~, ~, HeaderR] = y_ReadAll('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/SignificantArea/SeedFC/Corrected/FC_R_ROI1.gii');
[VBrain,HeaderV] = y_Read('/mnt/Data6/RfMRILab/Lubin/Software/DPABI_V6.0_210501/Templates/Tian2020_Subcortex_Atlas/Tian_Subcortex_S4_3T_2009cAsym.nii');

OutputDir = '/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Gradient/PostHoc';
mkdir(OutputDir);
Pairs = {'ASD_NC','SCZ_NC','ASD_SCZ'};
for iGrad = 1:2
    for iPair = 1:length(Pairs)
        if iPair ~= 2
            TValues = PostHoc_T_Gradient(iGrad,:,iPair);
        else
            TValues = -PostHoc_T_Gradient(iGrad,:,iPair); % reverse to SCZ-NC
        end
        GradBrainL = zeros(length(LeftTemplate),1);
        GradBrainR = zeros(length(RightTemplate),1);
        GradBrainV = zeros(size(VBrain));
        for iROI = 1:200
            if PostHoc_P_Gradient(iGrad,iROI,iPair)<0.05/6 %%% multiple compare!
                GradBrainL(find(LeftTemplate==iROI)) = TValues(iROI);
            end
        end
        for iROI = 201:400
            if PostHoc_P_Gradient(iGrad,iROI,iPair)<0.05/6
                GradBrainR(find(RightTemplate==iROI-200)) = TValues(iROI);
            end
        end
        for iROI = 401:454
            if PostHoc_P_Gradient(iGrad,iROI,iPair)<0.05/6
                GradBrainV(find(VBrain==iROI-400)) = TValues(iROI);
            end
        end
        y_Write(GradBrainL,HeaderL,[OutputDir,filesep,'PostHoc_Gradient',num2str(iGrad),'_',Pairs{iPair},'_Left']);
        y_Write(GradBrainR,HeaderR,[OutputDir,filesep,'PostHoc_Gradient',num2str(iGrad),'_',Pairs{iPair},'_Right']);
        y_Write(GradBrainV,HeaderV,[OutputDir,filesep,'PostHoc_Gradient',num2str(iGrad),'_',Pairs{iPair},'_Volume']);
    end
end



%% Post-hoc - Graph theory
clc;clear;
load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat');

%%% Chekc GTAFlag! %%%
GTAIndex = find(SubInfo.Dx<=3 & SubInfo.QCFlag & SubInfo.Age <=18 & SubInfo.RestFlag & SubInfo.MotionFD<0.2 & SubInfo.GTAFlag);

ID = SubInfo.ID(GTAIndex);
Dx = SubInfo.Dx(GTAIndex);
Site = SubInfo.Site(GTAIndex);
Age = SubInfo.Age(GTAIndex);
Sex = SubInfo.Sex(GTAIndex);
MotionFD = SubInfo.MotionFD(GTAIndex);

FieldList1 = ROIValue.Graph.FieldList1;
FieldList2 = ROIValue.Graph.FieldList2;
Pairs = [1,2;2,3;1,3]; % ASD-NC, NC-SCZ, ASD-SCZ
Labels = {'ASD_NC', 'NC_SCZ', 'ASD_SCZ'};
for iPair = 1:size(Pairs,1)
    GroupIndex1 = find(Dx ==Pairs(iPair,1));
    GroupIndex2 = find(Dx ==Pairs(iPair,2));
    
    DxCov = [ones(length(GroupIndex1),1);ones(length(GroupIndex2),1).*-1];
    AgeCov = [Age(GroupIndex1);Age(GroupIndex2)];
    SexCov = [Sex(GroupIndex1);Sex(GroupIndex2)];
    MotionCov = [MotionFD(GroupIndex1);MotionFD(GroupIndex2)];
    Site1 = [Site(GroupIndex1);Site(GroupIndex2)];
    SiteCov=[];
    SiteIndex = unique(Site1);
    for i=1:length(SiteIndex)-1
        SiteCov(:,i) = Site1==SiteIndex(i);
    end
    
    %%% Single values should not use combat results %%%
    AllCov = [DxCov,ones(length(DxCov),1),AgeCov,SexCov,MotionCov,SiteCov];
    Contrast=zeros(1,size(AllCov,2));
    Contrast(1)=1;
    
    for iMetric = 1:length(FieldList1)
            GTA_Value = ROIValue.Graph.ValueRaw.(FieldList1{iMetric})(GTAIndex);
            GTA_Value = GTA_Value(:,[GroupIndex1;GroupIndex2])';
            if ROIValue.Graph.Ancova.P.(FieldList1{iMetric}) <0.05
                [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(GTA_Value,AllCov,Contrast,'T');
                Stat_GTA_PostHoc.(Labels{iPair}).T.(FieldList1{iMetric}) = TF_ForContrast;
                if TF_ForContrast>=0
                    Stat_GTA_PostHoc.(Labels{iPair}).P.(FieldList1{iMetric}) = 1-tcdf(TF_ForContrast,length(DxCov)-2);
                else
                    Stat_GTA_PostHoc.(Labels{iPair}).P.(FieldList1{iMetric}) = tcdf(TF_ForContrast,length(DxCov)-2);
                end
            else
                Stat_GTA_PostHoc.(Labels{iPair}).T.(FieldList1{iMetric}) = 0;
                Stat_GTA_PostHoc.(Labels{iPair}).P.(FieldList1{iMetric}) = 1;
            end
    end
    
    %%% Dependent values should not use combat results - EglobSet/ElocSet %%%
    for iMetric = 8:9 % EglobSet/ElocSet
        GTA_Value = ROIValue.Graph.ValueRaw.(FieldList2{iMetric})(:,GTAIndex);
        GTA_Value = GTA_Value(:,[GroupIndex1;GroupIndex2])';
        for i = 1:size(GTA_Value,2)
            if ROIValue.Graph.Ancova.P_Corrected.(FieldList2{iMetric})(i) <0.05
                [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(GTA_Value(:,i),AllCov,Contrast,'T');
                Stat_GTA_PostHoc.(Labels{iPair}).T.(FieldList2{iMetric})(i) = TF_ForContrast;
                if TF_ForContrast>=0
                    Stat_GTA_PostHoc.(Labels{iPair}).P.(FieldList2{iMetric})(i) = 1-tcdf(TF_ForContrast,length(DxCov)-2);
                else
                    Stat_GTA_PostHoc.(Labels{iPair}).P.(FieldList2{iMetric})(i) = tcdf(TF_ForContrast,length(DxCov)-2);
                end
            else
                Stat_GTA_PostHoc.(Labels{iPair}).T.(FieldList2{iMetric})(i) = 0;
                Stat_GTA_PostHoc.(Labels{iPair}).P.(FieldList2{iMetric})(i) = 1;
            end
        end
    end
    
    %%% Multiple-ROI values should use combat results %%%
    AllCov = [DxCov,ones(length(DxCov),1),AgeCov,SexCov,MotionCov]; % No SiteCov!
    Contrast=zeros(1,size(AllCov,2));
    Contrast(1)=1;
    
    for iMetric = 1:length(FieldList2)-2 % no EglobSet/ElocSet
        GTA_Value = ROIValue.Graph.ValueCombat.(FieldList2{iMetric})(:,GTAIndex);
        GTA_Value = GTA_Value(:,[GroupIndex1;GroupIndex2])';
        for i = 1:size(GTA_Value,2)
            if ROIValue.Graph.Ancova.P_Corrected.(FieldList2{iMetric})(i) <0.05
                [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(GTA_Value(:,i),AllCov,Contrast,'T');
                Stat_GTA_PostHoc.(Labels{iPair}).T.(FieldList2{iMetric})(i) = TF_ForContrast;
                if TF_ForContrast>=0
                    Stat_GTA_PostHoc.(Labels{iPair}).P.(FieldList2{iMetric})(i) = 1-tcdf(TF_ForContrast,length(DxCov)-2);
                else
                    Stat_GTA_PostHoc.(Labels{iPair}).P.(FieldList2{iMetric})(i) = tcdf(TF_ForContrast,length(DxCov)-2);
                end
            else
                Stat_GTA_PostHoc.(Labels{iPair}).T.(FieldList2{iMetric})(i) = 0;
                Stat_GTA_PostHoc.(Labels{iPair}).P.(FieldList2{iMetric})(i) = 1;
            end
        end
    end
    
    %%% NodalEfficiency_gradient %%% 
    GTA_Value = ROIValue.Graph.ValueCombat.NodalEfficiency_AUC(:,GTAIndex);
    GTA_Value = GTA_Value(:,[GroupIndex1;GroupIndex2])';
    for i = 1:size(GTA_Value,2)
        if ROIValue.Graph.Ancova.P_Corrected.NodalEffi_Gradient(i) <0.05
            [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(GTA_Value(:,i),AllCov,Contrast,'T');
            Stat_GTA_PostHoc.(Labels{iPair}).T.NodalEffi_Gradient(i) = TF_ForContrast;
            if TF_ForContrast>=0
                Stat_GTA_PostHoc.(Labels{iPair}).P.NodalEffi_Gradient(i) = 1-tcdf(TF_ForContrast,length(DxCov)-2);
            else
                Stat_GTA_PostHoc.(Labels{iPair}).P.NodalEffi_Gradient(i) = tcdf(TF_ForContrast,length(DxCov)-2);
            end
        else
            Stat_GTA_PostHoc.(Labels{iPair}).T.NodalEffi_Gradient(i) = 0;
            Stat_GTA_PostHoc.(Labels{iPair}).P.NodalEffi_Gradient(i) = 1;
        end
    end
end

ROIValue.Graph.PostHoc = Stat_GTA_PostHoc;
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','ROIValue','SubInfo');






%% %%%%%% Classifier %%%%%%

%% Classifier - PCA - 5-fold randon cross-validation
clc;clear;
AnatData = load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat');
FuncData = load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat');

FuncIndex = find(FuncData.SubInfo.Dx<=3 & FuncData.SubInfo.QCFlag & FuncData.SubInfo.Age <=18 ...
    & FuncData.SubInfo.RestFlag & FuncData.SubInfo.MotionFD<0.2 & FuncData.SubInfo.GTAFlag);

FuncID = FuncData.SubInfo.ID(FuncIndex);
[ID,IA,IB] = intersect(AnatData.SubInfo.ID,FuncID);
FuncIndex = FuncIndex(IB);

% RandIndex = randperm(length(ID))';
% FuncData.SubInfo.RandIndexSVM = RandIndex;
% SubInfo = FuncData.SubInfo;
% ROIValue = FuncData.ROIValue;
% save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat','ROIValue','SubInfo');

%%% Extract raw features %%%
% Phenotype
Sex = FuncData.SubInfo.Sex(FuncIndex);
Dx = FuncData.SubInfo.Dx(FuncIndex);
Age = FuncData.SubInfo.Age(FuncIndex);
MotionFD = FuncData.SubInfo.MotionFD(FuncIndex);
ICV = AnatData.SubInfo.ICV(IA);

% Funcitonal metrics
Data_SeedFC = FuncData.ROIValue.SeedFC.FC(FuncIndex,:);
Data_NetFC = FuncData.ROIValue.NetworkFC.SigFC_FDR05(FuncIndex,:);
Data_ROIFC = FuncData.ROIValue.ROI2ROIFC.FC_Sig(FuncIndex,:);
Data_Connetcome = FuncData.ROIValue.Connetcome.SigFC_FDR005(FuncIndex,:);
SigInd = find(FuncData.ROIValue.Graph.Ancova.P_Corrected.NodalEffi_Gradient<0.05);
Data_Effi = FuncData.ROIValue.Graph.ValueCombat.NodalEfficiency_AUC(SigInd,FuncIndex)';
SigInd = find(FuncData.ROIValue.Graph.Ancova.P_Corrected.Degree_AUC<0.05);
Data_Degree = FuncData.ROIValue.Graph.ValueCombat.Degree_AUC(SigInd,FuncIndex)';
SigInd = find(FuncData.ROIValue.Gradient.Ancova.P_Grad_ANCOVA_Corrected(1,:)<0.05); % Grad1
Data_Grad1 = squeeze(FuncData.ROIValue.Gradient.GradValueCombat(1,SigInd,FuncIndex))';
SigInd = find(FuncData.ROIValue.Gradient.Ancova.P_Grad_ANCOVA_Corrected(2,:)<0.05); % Grad2
Data_Grad2 = squeeze(FuncData.ROIValue.Gradient.GradValueCombat(2,SigInd,FuncIndex))';

% Anatomical metrics
Data_Thickness = AnatData.ROIValue.CorticalMetrics_Sig(IA,:);
Data_Volume = AnatData.ROIValue.SubcortVol_Sig(IA,:);

%%% Do PCA to generate features %%%
% Cortical thickness
[COEFF, SCORE] = pca(Data_Thickness);
Feature_Thickness = SCORE(:,1);

% Subcortical volume
[COEFF, SCORE] = pca(Data_Volume);
Feature_Volume = SCORE(:,1);

% Seed-based FC
[COEFF, SCORE] = pca([Data_SeedFC,Data_ROIFC]);
Feature_FC = SCORE(:,1);

% Connectome
[COEFF, SCORE] = pca([Data_NetFC,Data_Connetcome]);
Feature_Connetcome = SCORE(:,1);

% Gradient 1
[COEFF, SCORE] = pca(Data_Grad1);
Feature_Gradient1 = SCORE(:,1);

% Gradient 2
[COEFF, SCORE] = pca(Data_Grad2);
Feature_Gradient2 = SCORE(:,1);

% Degree
[COEFF, SCORE] = pca(Data_Degree);
Feature_Degree = SCORE(:,1);

% Efficiency
[COEFF, SCORE] = pca(Data_Effi);
Feature_Efficiency = SCORE(:,1);


%%% Five-fold random cross-validation %%%
FeatureAll = [Age,Sex,ICV,MotionFD,Feature_Thickness,Feature_Volume,...
    Feature_FC,Feature_Connetcome,Feature_Gradient1,Feature_Gradient2,...
    Feature_Degree,Feature_Efficiency];
LabelAll = Dx;

% PermIndex = randperm(length(Age));
PermIndex = FuncData.SubInfo.RandIndexSVM;
FeatureAll = FeatureAll(PermIndex,:);
LabelAll = LabelAll(PermIndex);

Pairs = [1,2;3,2;1,3];
PairNames = {'ASD_NC','SCZ_NC','ASD_SCZ'};
for iPair = 1:size(Pairs,1)
    PairIndex = find((LabelAll == Pairs(iPair,1) | LabelAll == Pairs(iPair,2)));
    FeatureSVM = FeatureAll(PairIndex,:);
    LabelRaw = LabelAll(PairIndex);
    LabelSVM = ones(length(PairIndex),1);
    LabelSVM(find(LabelRaw==Pairs(iPair,2)))=0;
    nSub = ceil(length(PairIndex)/5);
    LabelPred = [];
    ScorePred = [];
    for iFold = 1:5
        TestStart = nSub*(iFold-1)+1;
        TestStop = min([nSub*iFold,length(PairIndex)]);
        TrainInd = setdiff([1:length(PairIndex)],[TestStart:TestStop]);
        Classifier.(PairNames{iPair}).(['Fold',num2str(iFold)]) = fitcsvm(FeatureSVM(TrainInd,:),LabelSVM(TrainInd),...
            'KernelFunction','linear','KernelScale','auto','Standardize',true,'ClassNames',[1,0]);%
        [PredTemp,ScoreTemp] = predict(Classifier.(PairNames{iPair}).(['Fold',num2str(iFold)]),FeatureSVM(TestStart:TestStop,:));
        
        [Confusion,~]=confusionmat(LabelSVM(TestStart:TestStop),PredTemp);
        Accuracy(iFold,iPair) = sum(diag(Confusion))/sum(Confusion(:));
        Sensitivity(iFold,iPair) = Confusion(2,2)/sum(Confusion(2,:));
        Specificity(iFold,iPair) = Confusion(1,1)/sum(Confusion(1,:));
        ConfusionMatrix.(PairNames{iPair}).(['Fold',num2str(iFold)]) = Confusion;
        [ROC_X.(PairNames{iPair}).(['Fold',num2str(iFold)]),ROC_Y.(PairNames{iPair}).(['Fold',num2str(iFold)]),~,ROC_AUC(iFold,iPair)]=...
            perfcurve(LabelSVM(TestStart:TestStop),ScoreTemp(:,1),1);
        
        LabelPred = [LabelPred;PredTemp];
        ScorePred = [ScorePred;ScoreTemp];
        disp(['Have down classifier_',PairNames{iPair},' fold',num2str(iFold)])
    end
    [Confusion,~]=confusionmat(LabelSVM,LabelPred);
    Accuracy_All(iPair) = sum(diag(Confusion))/sum(Confusion(:));
    Sensitivity_All(iPair) = Confusion(2,2)/sum(Confusion(2,:));
    Specificity_All(iPair) = Confusion(1,1)/sum(Confusion(1,:));
    ConfusionMatrix.(PairNames{iPair}).All = Confusion;
    [ROC_X.(PairNames{iPair}).All,ROC_Y.(PairNames{iPair}).All,~,ROC_AUC_All(iPair)]=perfcurve(LabelSVM,ScorePred(:,1),1);
end

save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Classifier/Classifier_LinearSVM_PCA_RandomCV.mat',...
    'Classifier','Accuracy','Sensitivity','Specificity','ROC_AUC','ROC_X','ROC_Y','ConfusionMatrix',...
    'Accuracy_All','Sensitivity_All','Specificity_All','ROC_AUC_All');

%%% save the PCA-profile of ASD, NC and SCZ %%%
for iDx = 1:3
    Index = find(Dx==iDx);
    Profile(:,iDx) = mean([Feature_Thickness(Index),Feature_Volume(Index),Feature_FC(Index),Feature_Connetcome(Index),...
        Feature_Gradient1(Index),Feature_Gradient2(Index),Feature_Degree(Index),Feature_Efficiency(Index)],1);
end

zProfile = (Profile-repmat(mean(Profile,2),1,3))./repmat(std(Profile,0,2),1,3);
Signs = repmat([1,-1,1,1,-1,-1,1,1]',1,3); % reverse PC direction according to post-hoc analysis for illustraction
zProfile = zProfile.*Signs;

FeatureName = {'Feature_Thickness','Feature_Volume','Feature_FC','Feature_Connetcome','Feature_Gradient1','Feature_Gradient2','Feature_Degree','Feature_Efficiency'};
save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Classifier/Profile_PCA.mat','Profile','zProfile','FeatureName');



%% Classifier - PCA - leave-site-out cross-validation
clc;clear;
AnatData = load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912.mat');
FuncData = load('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/Data/SubInfo_20230912_Func.mat');

FuncIndex = find(FuncData.SubInfo.Dx<=3 & FuncData.SubInfo.QCFlag & FuncData.SubInfo.Age <=18 ...
    & FuncData.SubInfo.RestFlag & FuncData.SubInfo.MotionFD<0.2 & FuncData.SubInfo.GTAFlag);

FuncID = FuncData.SubInfo.ID(FuncIndex);
[ID,IA,IB] = intersect(AnatData.SubInfo.ID,FuncID);
FuncIndex = FuncIndex(IB);

%%% Extract raw features %%%
% Phenotype
Sex = FuncData.SubInfo.Sex(FuncIndex);
Dx = FuncData.SubInfo.Dx(FuncIndex);
Age = FuncData.SubInfo.Age(FuncIndex);
MotionFD = FuncData.SubInfo.MotionFD(FuncIndex);
Site = FuncData.SubInfo.Site(FuncIndex);
ICV = AnatData.SubInfo.ICV(IA);

% Funcitonal metrics
Data_SeedFC = FuncData.ROIValue.SeedFC.FC(FuncIndex,:);
Data_NetFC = FuncData.ROIValue.NetworkFC.SigFC_FDR05(FuncIndex,:);
Data_ROIFC = FuncData.ROIValue.ROI2ROIFC.FC_Sig(FuncIndex,:);
Data_Connetcome = FuncData.ROIValue.Connetcome.SigFC_FDR005(FuncIndex,:);
SigInd = find(FuncData.ROIValue.Graph.Ancova.P_Corrected.NodalEffi_Gradient<0.05);
Data_Effi = FuncData.ROIValue.Graph.ValueCombat.NodalEfficiency_AUC(SigInd,FuncIndex)';
SigInd = find(FuncData.ROIValue.Graph.Ancova.P_Corrected.Degree_AUC<0.05);
Data_Degree = FuncData.ROIValue.Graph.ValueCombat.Degree_AUC(SigInd,FuncIndex)';
SigInd = find(FuncData.ROIValue.Gradient.Ancova.P_Grad_ANCOVA_Corrected(1,:)<0.05); % Grad1
Data_Grad1 = squeeze(FuncData.ROIValue.Gradient.GradValueCombat(1,SigInd,FuncIndex))';
SigInd = find(FuncData.ROIValue.Gradient.Ancova.P_Grad_ANCOVA_Corrected(2,:)<0.05); % Grad2
Data_Grad2 = squeeze(FuncData.ROIValue.Gradient.GradValueCombat(2,SigInd,FuncIndex))';

% Anatomical metrics
Data_Thickness = AnatData.ROIValue.CorticalMetrics_Sig(IA,:);
Data_Volume = AnatData.ROIValue.SubcortVol_Sig(IA,:);

%%% Do PCA to generate features %%%
% Cortical thickness
[COEFF, SCORE] = pca(Data_Thickness);
Feature_Thickness = SCORE(:,1);

% Subcortical volume
[COEFF, SCORE] = pca(Data_Volume);
Feature_Volume = SCORE(:,1);

% Seed-based FC
[COEFF, SCORE] = pca([Data_SeedFC,Data_ROIFC]);
Feature_FC = SCORE(:,1);

% Connectome
[COEFF, SCORE] = pca([Data_NetFC,Data_Connetcome]);
Feature_Connetcome = SCORE(:,1);

% Gradient 1
[COEFF, SCORE] = pca(Data_Grad1);
Feature_Gradient1 = SCORE(:,1);

% Gradient 2
[COEFF, SCORE] = pca(Data_Grad2);
Feature_Gradient2 = SCORE(:,1);

% Degree
[COEFF, SCORE] = pca(Data_Degree);
Feature_Degree = SCORE(:,1);

% Efficiency
[COEFF, SCORE] = pca(Data_Effi);
Feature_Efficiency = SCORE(:,1);

%%% Leave-site-out cross-validation %%%
FeatureAll = [Age,Sex,ICV,MotionFD,Feature_Thickness,Feature_Volume,...
    Feature_FC,Feature_Connetcome,Feature_Gradient1,Feature_Gradient2,...
    Feature_Degree,Feature_Efficiency];
LabelAll = Dx;

Pairs = [1,2;3,2;1,3];
PairNames = {'ASD_NC','SCZ_NC','ASD_SCZ'};
SiteLists = {[1,3,4],[1,3,4],[1,2,3,4]};% discard small sites, ASD-NC 3 site(folds), NC-SCZ 3 site(folds), ASD-SCZ 4 site(folds)
for iPair = 1:size(Pairs,1)
    PairIndex = find((LabelAll == Pairs(iPair,1) | LabelAll == Pairs(iPair,2))...
        & any(Site == SiteLists{iPair},2));
    SiteSVM = Site(PairIndex);
    FeatureSVM = FeatureAll(PairIndex,:);
    LabelRaw = LabelAll(PairIndex);
    LabelSVM = ones(length(PairIndex),1);
    LabelSVM(find(LabelRaw==Pairs(iPair,2)))=0;
    nSub = ceil(length(PairIndex)/5);
    LabelPred = [];
    ScorePred = [];
    LabelSVM_New = [];
    for iFold = 1:length(SiteLists{iPair})
        TestInd = find(SiteSVM==SiteLists{iPair}(iFold));
        TrainInd = find(SiteSVM~=SiteLists{iPair}(iFold));
        Classifier.(PairNames{iPair}).(['Fold',num2str(iFold)]) = fitcsvm(FeatureSVM(TrainInd,:),LabelSVM(TrainInd),...
            'KernelFunction','linear','KernelScale','auto','Standardize',true,'ClassNames',[1,0]);%
        [PredTemp,ScoreTemp] = predict(Classifier.(PairNames{iPair}).(['Fold',num2str(iFold)]),FeatureSVM(TestInd,:));
        
        [Confusion,~]=confusionmat(LabelSVM(TestInd),PredTemp);
        Accuracy(iFold,iPair) = sum(diag(Confusion))/sum(Confusion(:));
        Sensitivity(iFold,iPair) = Confusion(2,2)/sum(Confusion(2,:));
        Specificity(iFold,iPair) = Confusion(1,1)/sum(Confusion(1,:));
        ConfusionMatrix.(PairNames{iPair}).(['Fold',num2str(iFold)]) = Confusion;
        [ROC_X.(PairNames{iPair}).(['Fold',num2str(iFold)]),ROC_Y.(PairNames{iPair}).(['Fold',num2str(iFold)]),~,ROC_AUC(iFold,iPair)]=...
            perfcurve(LabelSVM(TestInd),ScoreTemp(:,1),1);
        
        LabelPred = [LabelPred;PredTemp];
        ScorePred = [ScorePred;ScoreTemp];
        LabelSVM_New = [LabelSVM_New;LabelSVM(TestInd)];
        disp(['Have down classifier_',PairNames{iPair},' fold',num2str(iFold)])
    end
    [Confusion,~]=confusionmat(LabelSVM_New,LabelPred);
    Accuracy_All(iPair) = sum(diag(Confusion))/sum(Confusion(:));
    Sensitivity_All(iPair) = Confusion(2,2)/sum(Confusion(2,:));
    Specificity_All(iPair) = Confusion(1,1)/sum(Confusion(1,:));
    ConfusionMatrix.(PairNames{iPair}).All = Confusion;
    [ROC_X.(PairNames{iPair}).All,ROC_Y.(PairNames{iPair}).All,~,ROC_AUC_All(iPair)]=perfcurve(LabelSVM_New,ScorePred(:,1),1);
end

save('/mnt/Data6/RfMRILab/Lubin/Project/ASD_SCZ/StatisticalResults/20230912_ICVCovariate/Classifier/Classifier_LinearSVM_PCA_LOSO.mat',...
    'Classifier','Accuracy','Sensitivity','Specificity','ROC_AUC','ROC_X','ROC_Y','ConfusionMatrix',...
    'Accuracy_All','Sensitivity_All','Specificity_All','ROC_AUC_All');



















