

AllDataDir = '/mnt/Data/RfMRILab/RfMRIMaps/Preprocessed/FCPIntermediateFiles';

%load ('/home/data/HeadMotion_YCG/YAN_Work/TransformationProject/SubIDSet_Excluded.mat');
load ('/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/SubInfo/SubIDSet_Excluded_HaveBeijingCambridge.mat');


for iSite=1:length(SiteName)
    
    load /mnt/Data/RfMRILab/RfMRIMaps/Preprocessed/FCPIntermediateFiles/DPARSFA8mmSetting.mat

    AutoDataProcessParameter=Cfg;
    
    AutoDataProcessParameter.DataProcessDir = [AllDataDir,filesep,SiteName{iSite}];
    
    SubID = SubIDSet_ExcludedT1{iSite};
    
    AutoDataProcessParameter.SubjectID=SubID;
    AutoDataProcessParameter.FunctionalSessionNumber=1; 
    
    AutoDataProcessParameter.TimePoints=SiteNTimePoints{iSite}; 
    AutoDataProcessParameter.TR=SiteTR{iSite}; 
 
    DPARSFA_run(AutoDataProcessParameter);
    DPARSFA_RerunWithGSR(AutoDataProcessParameter);
end







%%%
%Get the design matrix

load('/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/SubInfo/SubID_SexAge_All18Sites_895Sub.mat');


Age=double(Age_AllSites);

Sex=zeros(length(SexFM_AllSites),1);
for i=1:length(SexFM_AllSites)
    if strcmpi(SexFM_AllSites{i},'f')
        Sex(i)=-1;
    elseif strcmpi(SexFM_AllSites{i},'m')
        Sex(i)=1;
    end
end

EyeStatus = double(EyeColosed_AllSites);
EyeStatus(find(EyeStatus==0)) = -1; %EC: 1; EO: -1.

TeslaStatus = double(Tesla15T_AllSites);
TeslaStatus(find(TeslaStatus==1)) = -1; %3T: 1; 1.5T: -1.
TeslaStatus(find(TeslaStatus==0)) = 1;  %3T: 1; 1.5T: -1.

Site=SiteID_AllSites;
%Note: NewYork_a, NewYork_b are two different sites now.
%Leiden_2180, Leiden_2200 are two different sites now.

SubID=SubID_AllSites;


%%%1. Deal With Age
AgeRange = [18 32];
%AgeRange = [7 14];

AgeTemp = (Age>=AgeRange(1)) .* (Age<=AgeRange(2));


%%%2. Deal with head motion
MeanFD = MeanFDJenkinson_AllSites;

MeanFDTemp = MeanFD<(mean(MeanFD)+2*std(MeanFD));

WantedSubMatrix = AgeTemp.*MeanFDTemp;



%%%3. Deal with bad coverage

CoverageTemp = VoxelPercentOverlapWithMask_90>(mean(VoxelPercentOverlapWithMask_90)-2*std(VoxelPercentOverlapWithMask_90));

WantedSubMatrix=WantedSubMatrix.*CoverageTemp';

% 
% %%%4. Exclude Sites
% NoBadSite=ones(length(SubID),1);
% NoBadSite(Site==12)=0;%Pitt 2
% NoBadSite(Site==3)=0;%Caltech
% 
% WantedSubMatrix=WantedSubMatrix.*NoBadSite;
% 


%4!!!! Deal with bad motion scrubbed >50%
%WantedSubMatrix=WantedSubMatrix.*(ScrubbedPercentage<0.5);

%5!!!! Error subject
WantedSubMatrix(109)=0;




WantedSubIndex = find(WantedSubMatrix);

%Select subjects
SubID=SubID(WantedSubIndex);

Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
MeanFD=MeanFD(WantedSubIndex);

EyeStatus=EyeStatus(WantedSubIndex);
TeslaStatus=TeslaStatus(WantedSubIndex);

Site=Site(WantedSubIndex);
SiteName_AllSites = SiteName_AllSites(WantedSubIndex);
SiteID_AllSites = SiteID_AllSites(WantedSubIndex);

SiteCov=[];
SiteIndex = unique(Site);
for i=1:length(SiteIndex)-1
    SiteCov(:,i) = Site==SiteIndex(i);
end

%Model 1:
%AllCov = [ones(length(SubID),1), Age,Sex,MeanFD,EyeStatus,TeslaStatus,SiteCov];

%Model 2: don't need EyeStatus, TeslaStatus. Rank deficit
AllCov = [ones(length(SubID),1), Age,Sex,MeanFD,SiteCov];

%Centering: Let the first column (constant) have the mean effect.
AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));

Contrast=zeros(1,size(AllCov,2));
Contrast(3)=1;











PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=1;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'

PALMSettings31=PALMSettings;
PALMSettings31.ClusterFormingThreshold=3.1;
PALMSettings258=PALMSettings;
PALMSettings258.ClusterFormingThreshold=2.58;
PALMSettings329=PALMSettings;
PALMSettings329.ClusterFormingThreshold=3.29;

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

DataDir='/mnt/Data/RfMRILab/RfMRIMaps/Preprocessed/FCPIntermediateFiles';
OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/ReAnalysis_FWMH8mm/FCP716MaleVsFemale';

MaskFile = ['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/SubInfo/FCP_GroupMask90Percent.nii'];


parfor iMeasure=1:length(MeasureSet)
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
    FileList=[];
    FileListGSR=[];
    for iSub=1:length(SubID)
        FileList{iSub,1}=[DataDir,'/',SiteName_AllSites{iSub},'/',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
        FileListGSR{iSub,1}=[DataDir,'/',SiteName_AllSites{iSub},'/',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
    end

    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleT'];
    % [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(DependentVolume,Predictor,OutputName,MaskFile,CovVolume,Contrast,TF_Flag,IsOutputResidual,Header)
    [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(FileList,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MaleVsFemaleT'];
    [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(FileListGSR,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleTPALM23'];
    y_GroupAnalysis_PermutationTest_Image(FileList,AllCov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings);
    %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/Temp'])
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MaleVsFemaleTPALM23'];
    y_GroupAnalysis_PermutationTest_Image(FileListGSR,AllCov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings);
    %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/Temp'])
    
    %PALMSettings.ClusterFormingThreshold=3.1;
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleTPALM31'];
    y_GroupAnalysis_PermutationTest_Image(FileList,AllCov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings31);
    %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/Temp'])
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MaleVsFemaleTPALM31'];
    y_GroupAnalysis_PermutationTest_Image(FileListGSR,AllCov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings31);
    %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/Temp'])
    
    %PALMSettings.ClusterFormingThreshold=3.1;
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleTPALM258'];
    y_GroupAnalysis_PermutationTest_Image(FileList,AllCov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings258);
    %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/Temp'])
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MaleVsFemaleTPALM258'];
    y_GroupAnalysis_PermutationTest_Image(FileListGSR,AllCov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings258);
    %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/Temp'])
    
    %PALMSettings.ClusterFormingThreshold=3.1;
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleTPALM329'];
    y_GroupAnalysis_PermutationTest_Image(FileList,AllCov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings329);
    %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/Temp'])
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MaleVsFemaleTPALM329'];
    y_GroupAnalysis_PermutationTest_Image(FileListGSR,AllCov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings329);
    %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/Temp'])
end


%%%%%%

OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/ReAnalysis_FWMH8mm/FCP716MaleVsFemale';%Copy all measures to ./ECEOPairedT
MaskFile = ['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/SubInfo/FCP_GroupMask90Percent.nii'];

%Then call the code in ECEOStats_BeijingZou.m While chaning ECEOPairedT to MaleVsFemaleT









%Get Reproducibility
CorrectionSet={'OneTailed_AFNI3dClustSim_0.01_0.05','OneTailed_AFNI3dClustSim_0.005_0.05','OneTailed_AFNI3dClustSim_0.001_0.05','OneTailed_AFNI3dClustSim_0.0005_0.05','OneTailed_AFNI3dClustSim_0.01_0.025','OneTailed_AFNI3dClustSim_0.005_0.025','OneTailed_AFNI3dClustSim_0.001_0.025','OneTailed_AFNI3dClustSim_0.0005_0.025','OneTailed_DPABIAlphaSim_0.01_0.05','OneTailed_DPABIAlphaSim_0.005_0.05','OneTailed_DPABIAlphaSim_0.001_0.05','OneTailed_DPABIAlphaSim_0.0005_0.05','OneTailed_DPABIAlphaSim_0.01_0.025','OneTailed_DPABIAlphaSim_0.005_0.025','OneTailed_DPABIAlphaSim_0.001_0.025','OneTailed_DPABIAlphaSim_0.0005_0.025','OneTailedGRF_0.01_0.05','OneTailedGRF_0.005_0.05','OneTailedGRF_0.001_0.05','OneTailedGRF_0.0005_0.05','OneTailedGRF_0.01_0.025','OneTailedGRF_0.005_0.025','OneTailedGRF_0.001_0.025','OneTailedGRF_0.0005_0.025','PALMCluster23','PALMCluster258','PALMCluster31','PALMCluster329','PALMTFCE','PALMVox','FDR'};

DataDirUp='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/ReAnalysis_FWMH8mm/FCP716MaleVsFemale/Binarized';
DataDirUp2='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/ReAnalysis_FWMH8mm/MaleVsFemaleT_420/MaleVsFemaleT/Binarized';
DataDirUp3='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/ReAnalysis_FWMH8mm/MaleVsFemaleT_420/S2_MaleVsFemaleT/Binarized';

ReproducibilitySet=zeros(10,length(CorrectionSet));
JaccardSet=zeros(10,length(CorrectionSet));
DiceSet=zeros(10,length(CorrectionSet));
for iCorrection=1:length(CorrectionSet)
    DataDir=[DataDirUp,'/',CorrectionSet{iCorrection}];
    DataDir2=[DataDirUp2,'/',CorrectionSet{iCorrection}];
    DataDir3=[DataDirUp3,'/',CorrectionSet{iCorrection}];
    
    MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
    ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
    ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
    Reproducibility=[];
    ReproducibilityGSR=[];
    Jaccard=[];
    Dice=[];
    JaccardGSR=[];
    DiceGSR=[];
    for iMeasure=1:length(MeasureSet)
        File=[DataDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'.nii'];
        Data=y_ReadRPI(File);
        File=[DataDir2,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'.nii'];
        Data2=y_ReadRPI(File);
        File=[DataDir3,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'.nii'];
        Data3=y_ReadRPI(File);
        
        Data23 = (Data2+Data3)==2;
        Data=Data+Data23;
        Reproducibility(iMeasure,1)=mean(Data(find(Data>0)));
        Jaccard(iMeasure,1)=length(find(Data==2))/length(find(Data>=1));
        Dice(iMeasure,1)=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));

        File=[DataDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'.nii'];
        Data=y_ReadRPI(File);
        File=[DataDir2,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'.nii'];
        Data2=y_ReadRPI(File);
        File=[DataDir3,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'.nii'];
        Data3=y_ReadRPI(File);
        
        Data23 = (Data2+Data3)==2;
        Data=Data+Data23;
        ReproducibilityGSR(iMeasure,1)=mean(Data(find(Data>0)));
        JaccardGSR(iMeasure,1)=length(find(Data==2))/length(find(Data>=1));
        DiceGSR(iMeasure,1)=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
    end
    
    Reproducibility=[Reproducibility;ReproducibilityGSR];
    Jaccard=[Jaccard;JaccardGSR];
    Dice=[Dice;DiceGSR];
    
    ReproducibilitySet(:,iCorrection)=Reproducibility;
    JaccardSet(:,iCorrection)=Jaccard;
    DiceSet(:,iCorrection)=Dice;
end

DiceSet_Selected=DiceSet(:,[8 16 24 25:end]);

save(['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/Reproducibility/Revision_8mm_MaleFemaleReplicability_CORR_BetweenSession_FCP.mat'],'ReproducibilitySet','JaccardSet','DiceSet','DiceSet_Selected');

[p tbl stats]=friedman(DiceSet)
[c m]=multcompare(stats)





