
%%%STATS!!!

load /mnt/Data/RfMRILab/RfMRIMaps/Processing/ECEO/ECEOProcessing/Dong_Beijing_dpabi/DPARSFA_AutoSave_2017_3_19_16_15.mat

SubID=Cfg.SubjectID;

SessionSet={'S2_','S3_'};
MotionDir='/mnt/Data/RfMRILab/RfMRIMaps/Processing/ECEO/ECEOProcessing/Dong_Beijing_dpabi/RealignParameter';
Motion=[];
for iSession=1:length(SessionSet)
    for iSub=1:length(SubID)
        Temp=load([MotionDir,'/',SubID{iSub},'/',SessionSet{iSession},'FD_Jenkinson_',SubID{iSub},'.txt']);
        Motion{iSession}(iSub,1)=mean(Temp);
    end
end




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

DataDir='/mnt/Data/RfMRILab/RfMRIMaps/Processing/ECEO/ECEOProcessing/Dong_Beijing_dpabi';
OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Dong_Beijing/ECEO_PairedT/ECEO_PairedT_8mm';

MaskFile ='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Dong_Beijing/SubInfo/GroupMask_90percent.nii';


parfor iMeasure=1:length(MeasureSet)
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
    FileList=[];
    FileListGSR=[];
    FileList_S2=[];
    FileListGSR_S2=[];
    for iSub=1:length(SubID)
        FileList{iSub,1}=[DataDir,'/S2_',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
        FileListGSR{iSub,1}=[DataDir,'/S2_',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
        FileList_S2{iSub,1}=[DataDir,'/S3_',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
        FileListGSR_S2{iSub,1}=[DataDir,'/S3_',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
    end

    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedT'];
    %[TTestPaired_T,Header] = y_TTestPaired_Image(DependentDirs,OutputName,MaskFile,CovariateDirs,OtherCovariates,PALMSettings)
    [TTestPaired_T,Header] = y_TTestPaired_Image({FileList;FileList_S2},OutputName,MaskFile,[],Motion);
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedT'];
    [TTestPaired_T,Header] = y_TTestPaired_Image({FileListGSR;FileListGSR_S2},OutputName,MaskFile,[],Motion);
    
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedTPALM23'];
    [TTestPaired_T,Header] = y_TTestPaired_Image({FileList;FileList_S2},OutputName,MaskFile,[],Motion,PALMSettings);
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedTPALM23'];
    [TTestPaired_T,Header] = y_TTestPaired_Image({FileListGSR;FileListGSR_S2},OutputName,MaskFile,[],Motion,PALMSettings);

    
    %PALMSettings.ClusterFormingThreshold=3.1;
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedTPALM31'];
    [TTestPaired_T,Header] = y_TTestPaired_Image({FileList;FileList_S2},OutputName,MaskFile,[],Motion,PALMSettings31);
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedTPALM31'];
    [TTestPaired_T,Header] = y_TTestPaired_Image({FileListGSR;FileListGSR_S2},OutputName,MaskFile,[],Motion,PALMSettings31);
    
    %PALMSettings.ClusterFormingThreshold=3.1;
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedTPALM258'];
    [TTestPaired_T,Header] = y_TTestPaired_Image({FileList;FileList_S2},OutputName,MaskFile,[],Motion,PALMSettings258);
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedTPALM258'];
    [TTestPaired_T,Header] = y_TTestPaired_Image({FileListGSR;FileListGSR_S2},OutputName,MaskFile,[],Motion,PALMSettings258);
    
    %PALMSettings.ClusterFormingThreshold=3.1;
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedTPALM329'];
    [TTestPaired_T,Header] = y_TTestPaired_Image({FileList;FileList_S2},OutputName,MaskFile,[],Motion,PALMSettings329);
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedTPALM329'];
    [TTestPaired_T,Header] = y_TTestPaired_Image({FileListGSR;FileListGSR_S2},OutputName,MaskFile,[],Motion,PALMSettings329);
end









%%%%%%

OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Dong_Beijing/ECEO_PairedT/ECEO_PairedT_8mm';%Copy all measures to ./ECEOPairedT
MaskFile ='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Dong_Beijing/SubInfo/GroupMask_90percent.nii';

%Then call the code in ECEOStats_BeijingZou.m
