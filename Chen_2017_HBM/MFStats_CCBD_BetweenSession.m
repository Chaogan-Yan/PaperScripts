
load('/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/SubInfo/SubInfo_420.mat')

SiteUnique=unique(Site);
SiteRegressor=zeros(size(Site,1),length(SiteUnique));
for iSite=1:length(SiteUnique)
    SiteRegressor(find(Site==SiteUnique(iSite)),iSite)=1;
end
SiteRegressor=SiteRegressor(:,1:end-1);

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

SessionSet={'','S2_'};

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

DataDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/ReAnalysis_FWMH8mm/AllMaps_8mm';
OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/ReAnalysis_FWMH8mm/MaleVsFemaleT_420';

MaskFile = '/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/SubInfo/GroupMask_90percent_429AfterExcluding.nii';

for iSession=1:length(SessionSet)
    mkdir([OutDir,'/',SessionSet{iSession},'MaleVsFemaleT']);
    parfor iMeasure=1:length(MeasureSet)
        mkdir([OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
        mkdir([OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
        FileList=[];
        FileListGSR=[];
        
        for iSub=1:length(SubID)
            FileList{iSub,1}=[DataDir,'/',SessionSet{iSession},ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
            FileListGSR{iSub,1}=[DataDir,'/',SessionSet{iSession},ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
        end
        Cov=[Sex,Age,Motion(:,iSession),ones(length(SubID),1),SiteRegressor];
        Contrast=zeros(1,size(Cov,2));
        Contrast(1)=1;
        OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleT'];
        % [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(DependentVolume,Predictor,OutputName,MaskFile,CovVolume,Contrast,TF_Flag,IsOutputResidual,Header)
        [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0);
        OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MaleVsFemaleT'];
        [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(FileListGSR,Cov,OutputName,MaskFile,[],Contrast,'T',0);

        OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleTPALM23'];
        y_GroupAnalysis_PermutationTest_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings);
        %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/Temp'])
        OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MaleVsFemaleTPALM23'];
        y_GroupAnalysis_PermutationTest_Image(FileListGSR,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings);
        %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/Temp'])
        
        %PALMSettings.ClusterFormingThreshold=3.1;
        OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleTPALM31'];
        y_GroupAnalysis_PermutationTest_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings31);
        %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/Temp'])
        OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MaleVsFemaleTPALM31'];
        y_GroupAnalysis_PermutationTest_Image(FileListGSR,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings31);
        %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/Temp'])
        
        %PALMSettings.ClusterFormingThreshold=3.1;
        OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleTPALM258'];
        y_GroupAnalysis_PermutationTest_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings258);
        %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/Temp'])
        OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MaleVsFemaleTPALM258'];
        y_GroupAnalysis_PermutationTest_Image(FileListGSR,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings258);
        %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/Temp'])
        
        %PALMSettings.ClusterFormingThreshold=3.1;
        OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleTPALM329'];
        y_GroupAnalysis_PermutationTest_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings329);
        %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/Temp'])
        OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MaleVsFemaleTPALM329'];
        y_GroupAnalysis_PermutationTest_Image(FileListGSR,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings329);
        %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/Temp'])
    end
end




%%%%%%

OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/ReAnalysis_FWMH8mm/MaleVsFemaleT_420/MaleVsFemaleT';%Copy all measures to ./ECEOPairedT
MaskFile = '/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/SubInfo/GroupMask_90percent_429AfterExcluding.nii';

%Then call the code in ECEOStats_BeijingZou.m While chaning ECEOPairedT to MaleVsFemaleT






%Get Reproducibility
CorrectionSet={'OneTailed_AFNI3dClustSim_0.01_0.05','OneTailed_AFNI3dClustSim_0.005_0.05','OneTailed_AFNI3dClustSim_0.001_0.05','OneTailed_AFNI3dClustSim_0.0005_0.05','OneTailed_AFNI3dClustSim_0.01_0.025','OneTailed_AFNI3dClustSim_0.005_0.025','OneTailed_AFNI3dClustSim_0.001_0.025','OneTailed_AFNI3dClustSim_0.0005_0.025','OneTailed_DPABIAlphaSim_0.01_0.05','OneTailed_DPABIAlphaSim_0.005_0.05','OneTailed_DPABIAlphaSim_0.001_0.05','OneTailed_DPABIAlphaSim_0.0005_0.05','OneTailed_DPABIAlphaSim_0.01_0.025','OneTailed_DPABIAlphaSim_0.005_0.025','OneTailed_DPABIAlphaSim_0.001_0.025','OneTailed_DPABIAlphaSim_0.0005_0.025','OneTailedGRF_0.01_0.05','OneTailedGRF_0.005_0.05','OneTailedGRF_0.001_0.05','OneTailedGRF_0.0005_0.05','OneTailedGRF_0.01_0.025','OneTailedGRF_0.005_0.025','OneTailedGRF_0.001_0.025','OneTailedGRF_0.0005_0.025','PALMCluster23','PALMCluster258','PALMCluster31','PALMCluster329','PALMTFCE','PALMVox','FDR'};

DataDirUp='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/ReAnalysis_FWMH8mm/MaleVsFemaleT_420/MaleVsFemaleT/Binarized';
DataDirUp2='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/ReAnalysis_FWMH8mm/MaleVsFemaleT_420/S2_MaleVsFemaleT/Binarized';

ReproducibilitySet=zeros(10,length(CorrectionSet));
JaccardSet=zeros(10,length(CorrectionSet));
DiceSet=zeros(10,length(CorrectionSet));
for iCorrection=1:length(CorrectionSet)
    DataDir=[DataDirUp,'/',CorrectionSet{iCorrection}];
    DataDir2=[DataDirUp2,'/',CorrectionSet{iCorrection}];
    
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
        Data=Data+Data2;
        Reproducibility(iMeasure,1)=mean(Data(find(Data>0)));
        Jaccard(iMeasure,1)=length(find(Data==2))/length(find(Data>=1));
        Dice(iMeasure,1)=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));

        File=[DataDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'.nii'];
        Data=y_ReadRPI(File);
        File=[DataDir2,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'.nii'];
        Data2=y_ReadRPI(File);
        Data=Data+Data2;
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

save(['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/Reproducibility/Revision_8mm_TRTReliability_CORR_BetweenSession.mat'],'ReproducibilitySet','JaccardSet','DiceSet','DiceSet_Selected');

[p tbl stats]=friedman(DiceSet)
[c m]=multcompare(stats)



