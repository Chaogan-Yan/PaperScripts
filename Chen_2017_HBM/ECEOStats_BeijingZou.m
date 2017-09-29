



%%%STATS!!!

load /mnt/Data/RfMRILab/RfMRIMaps/Processing/ECEO/ECEOProcessing/Zou_Beijing_dpabi/DPARSFA_AutoSave_2017_5_5_14_58.mat

SubID=Cfg.SubjectID;

SessionSet={'','S2_'};
MotionDir='/mnt/Data/RfMRILab/RfMRIMaps/Processing/ECEO/ECEOProcessing/Zou_Beijing_dpabi/RealignParameter';
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

DataDir='/mnt/Data/RfMRILab/RfMRIMaps/Processing/ECEO/ECEOProcessing/Zou_Beijing_dpabi';
OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Zou_Beijing/ECEO_PairedT/ECEO_PairedT8mm';

MaskFile ='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Zou_Beijing/SubInfo/GroupMask_90percent.nii';


parfor iMeasure=1:length(MeasureSet)
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
    FileList=[];
    FileListGSR=[];
    FileList_S2=[];
    FileListGSR_S2=[];
    for iSub=1:length(SubID)
        FileList{iSub,1}=[DataDir,'/',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
        FileListGSR{iSub,1}=[DataDir,'/',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
        FileList_S2{iSub,1}=[DataDir,'/S2_',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
        FileListGSR_S2{iSub,1}=[DataDir,'/S2_',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
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









%DO
OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Zou_Beijing/ECEO_PairedT/ECEO_PairedT8mm';
MaskFile ='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Zou_Beijing/SubInfo/GroupMask_90percent.nii';


%%%As a function: BEGIN

OutDir_AlphaSim=[OutDir,'/AlphaSim'];
mkdir(OutDir_AlphaSim)
%Get smoothness
%%%

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

MaskData = y_ReadRPI(MaskFile);
MaskIndex = find(MaskData);
nVoxels = length(MaskIndex);

FWHMAllSet=[];
GSRFWHMAllSet=[];
for iMeasure=1:length(MeasureSet)
    FWHMAll=zeros(3,1);
    GSRFWHMAll=zeros(3,1);
    
    FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedT'];
    [Data Header]=y_Read(FileName);
    Header=w_ReadDLH(Header);
    FWHMAll(:,1)=[Header.FWHMx,Header.FWHMy,Header.FWHMz]';
    
    FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedT'];
    [Data Header]=y_Read(FileName);
    Header=w_ReadDLH(Header);
    GSRFWHMAll(:,1)=[Header.FWHMx,Header.FWHMy,Header.FWHMz]';
    
    FWHMAllSet(:,iMeasure)=FWHMAll;
    GSRFWHMAllSet(:,iMeasure)=GSRFWHMAll;
end
save([OutDir_AlphaSim,'/FWHM.mat'],'FWHMAllSet','GSRFWHMAllSet')



%DPABI AlphaSim
load([OutDir_AlphaSim,'/FWHM.mat'])

ClusterSize_OneTailed_NN26Set_AllMeasures=zeros(4,4,5);
ClusterSize_TwoTailed_NN26Set_AllMeasures=zeros(4,4,5);
GSRClusterSize_OneTailed_NN26Set_AllMeasures=zeros(4,4,5);
GSRClusterSize_TwoTailed_NN26Set_AllMeasures=zeros(4,4,5);

parfor iMeasure=1:5
    
    FWHMAll=FWHMAllSet(:,iMeasure);
    GSRFWHMAll=GSRFWHMAllSet(:,iMeasure);
    
    ClusterSize_OneTailed_NN26Set=zeros(4,4);
    ClusterSize_TwoTailed_NN26Set=zeros(4,4);
    GSRClusterSize_OneTailed_NN26Set=zeros(4,4);
    GSRClusterSize_TwoTailed_NN26Set=zeros(4,4);
    
    ClusterSize_OneTailed_NN26Temp=[];
    ClusterSize_TwoTailed_NN26Temp=[];
    [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.01,1000,'fwhm',FWHMAll(:,1)');
    ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
    ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
    
    [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.005,1000,'fwhm',FWHMAll(:,1)');
    ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
    ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
    
    [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.001,1000,'fwhm',FWHMAll(:,1)');
    ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
    ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
    
    [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.0005,1000,'fwhm',FWHMAll(:,1)');
    ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
    ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
    
    ClusterSize_OneTailed_NN26Set(:,:,1)=ClusterSize_OneTailed_NN26Temp;
    ClusterSize_TwoTailed_NN26Set(:,:,1)=ClusterSize_TwoTailed_NN26Temp;
    

    ClusterSize_OneTailed_NN26Temp=[];
    ClusterSize_TwoTailed_NN26Temp=[];
    [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.01,1000,'fwhm',GSRFWHMAll(:,1)');
    ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
    ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
    
    [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.005,1000,'fwhm',GSRFWHMAll(:,1)');
    ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
    ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
    
    [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.001,1000,'fwhm',GSRFWHMAll(:,1)');
    ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
    ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
    
    [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.0005,1000,'fwhm',GSRFWHMAll(:,1)');
    ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
    ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
    
    GSRClusterSize_OneTailed_NN26Set(:,:,1)=ClusterSize_OneTailed_NN26Temp;
    GSRClusterSize_TwoTailed_NN26Set(:,:,1)=ClusterSize_TwoTailed_NN26Temp;
    
    
    ClusterSize_OneTailed_NN26Set_AllMeasures(:,:,iMeasure)=ClusterSize_OneTailed_NN26Set;
    ClusterSize_TwoTailed_NN26Set_AllMeasures(:,:,iMeasure)=ClusterSize_TwoTailed_NN26Set;
    
    GSRClusterSize_OneTailed_NN26Set_AllMeasures(:,:,iMeasure)=GSRClusterSize_OneTailed_NN26Set;
    GSRClusterSize_TwoTailed_NN26Set_AllMeasures(:,:,iMeasure)=GSRClusterSize_TwoTailed_NN26Set;
end
save([OutDir_AlphaSim,'/ClusterSize_DPABIAlphaSim.mat'],'ClusterSize_OneTailed_NN26Set_AllMeasures','ClusterSize_TwoTailed_NN26Set_AllMeasures','GSRClusterSize_OneTailed_NN26Set_AllMeasures','GSRClusterSize_TwoTailed_NN26Set_AllMeasures')




%Apply DPABI AlphaSim Correction
load([OutDir_AlphaSim,'/ClusterSize_DPABIAlphaSim.mat'])

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

IsTwoTailed=1; %!!!
VoxelPThresholdSet_OneTailed=[0.01 0.005 0.001 0.0005 0.01 0.005 0.001 0.0005];
ClusterPThresholdSet_OneTailed=[0.05 0.05 0.05 0.05 0.025 0.025 0.025 0.025];

for iMeasure=1:length(MeasureSet)
    
    for iP=1:length(VoxelPThresholdSet_OneTailed)
        VoxelPThreshold=2*VoxelPThresholdSet_OneTailed(iP); %For y_ApplyAlphaSimThreshold, input two-tailed p values.
        
        iRow=1;
        iColumn=iP;
        if iP>=5
            iRow=2;
            iColumn=iP-4;
        end
        
        CorrectionName=['OneTailed_DPABIAlphaSim_',num2str(VoxelPThresholdSet_OneTailed(iP)),'_',num2str(ClusterPThresholdSet_OneTailed(iP))];
        mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/']);
        mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/']);
        
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedT'];
        OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/ECEOPairedT'];
        ClusterSize_OneTailed_NN26Set = ClusterSize_OneTailed_NN26Set_AllMeasures(:,:,iMeasure);
        ClusterSize=ClusterSize_OneTailed_NN26Set(iRow,iColumn);
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,OutName,MaskFile);

        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedT'];
        OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/ECEOPairedT'];
        ClusterSize_OneTailed_NN26Set = GSRClusterSize_OneTailed_NN26Set_AllMeasures(:,:,iMeasure);
        ClusterSize=ClusterSize_OneTailed_NN26Set(iRow,iColumn);
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,OutName,MaskFile);
    end
end










%AFNI 3dClustSim
%3dClustSim -mask /mnt/diske/Yan/Analysis/MultipleComparison/FCPData/FCP_GroupMask90Percent.nii -pthr 0.01 0.005 0.001 0.0005 -athr 0.05 0.025 -iter 1000 -fwhmxyz ${FWHM} -prefix /mnt/diske/Yan/Analysis/MultipleComparison/FCPData/Beijing/AlphaSim/Revision/AFNI3dClusterSim/FWHM_1_Rand${i}`
load([OutDir_AlphaSim,'/FWHM.mat'])
mkdir([OutDir_AlphaSim,'/AFNI3dClusterSim'])

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

for iMeasure=1:5
    FWHMAll=FWHMAllSet(:,iMeasure);
    GSRFWHMAll=GSRFWHMAllSet(:,iMeasure);

    eval(['!3dClustSim -mask ',MaskFile,' -pthr 0.01 0.005 0.001 0.0005 -athr 0.05 0.025 -iter 1000 -fwhmxyz ',num2str(FWHMAll(1)),' ',num2str(FWHMAll(2)),' ',num2str(FWHMAll(3)),' -prefix ',OutDir_AlphaSim,'/AFNI3dClusterSim/',MeasureSet{iMeasure},ConditionSet{iMeasure}])
    eval(['!3dClustSim -mask ',MaskFile,' -pthr 0.01 0.005 0.001 0.0005 -athr 0.05 0.025 -iter 1000 -fwhmxyz ',num2str(GSRFWHMAll(1)),' ',num2str(GSRFWHMAll(2)),' ',num2str(GSRFWHMAll(3)),' -prefix ',OutDir_AlphaSim,'/AFNI3dClusterSim/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}])
end



%Get the thresholds for AFNI 3dClustSim
ClusterSize_OneTailed_NN26Set_AllMeasures=zeros(2,4,5);
ClusterSize_TwoTailed_NN26Set_AllMeasures=zeros(2,4,5);
GSRClusterSize_OneTailed_NN26Set_AllMeasures=zeros(2,4,5);
GSRClusterSize_TwoTailed_NN26Set_AllMeasures=zeros(2,4,5);
ClusterSimDataDir=[OutDir_AlphaSim,'/AFNI3dClusterSim'];
for iMeasure=1:5
    ClusterSize_OneTailed_NN26Set=zeros(2,4);
    ClusterSize_TwoTailed_NN26Set=zeros(2,4);
    GSRClusterSize_OneTailed_NN26Set=zeros(2,4);
    GSRClusterSize_TwoTailed_NN26Set=zeros(2,4);
    
    File=[ClusterSimDataDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'.NN3_1sided.1D'];
    fid = fopen(File);
    StringFilter = '%f%f%f\n';
    %StringFilter = [StringFilter,'%*[^\n]']; %Skip the else till end of the line
    for i=1:8
        tline = fgetl(fid); %Skip the title line
    end
    ThresholdTemp = textscan(fid,StringFilter);
    fclose(fid);
    Threshold=[ThresholdTemp{2},ThresholdTemp{3}]';
    ClusterSize_OneTailed_NN26Set(:,:,1)=Threshold;
    
    File=[ClusterSimDataDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'.NN3_2sided.1D'];
    fid = fopen(File);
    StringFilter = '%f%f%f\n';
    %StringFilter = [StringFilter,'%*[^\n]']; %Skip the else till end of the line
    for i=1:8
        tline = fgetl(fid); %Skip the title line
    end
    ThresholdTemp = textscan(fid,StringFilter);
    fclose(fid);
    Threshold=[ThresholdTemp{2},ThresholdTemp{3}]';
    ClusterSize_TwoTailed_NN26Set(:,:,1)=Threshold;
    
    File=[ClusterSimDataDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'.NN3_1sided.1D'];
    fid = fopen(File);
    StringFilter = '%f%f%f\n';
    %StringFilter = [StringFilter,'%*[^\n]']; %Skip the else till end of the line
    for i=1:8
        tline = fgetl(fid); %Skip the title line
    end
    ThresholdTemp = textscan(fid,StringFilter);
    fclose(fid);
    Threshold=[ThresholdTemp{2},ThresholdTemp{3}]';
    GSRClusterSize_OneTailed_NN26Set(:,:,1)=Threshold;
    
    File=[ClusterSimDataDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'.NN3_2sided.1D'];
    fid = fopen(File);
    StringFilter = '%f%f%f\n';
    %StringFilter = [StringFilter,'%*[^\n]']; %Skip the else till end of the line
    for i=1:8
        tline = fgetl(fid); %Skip the title line
    end
    ThresholdTemp = textscan(fid,StringFilter);
    fclose(fid);
    Threshold=[ThresholdTemp{2},ThresholdTemp{3}]';
    GSRClusterSize_TwoTailed_NN26Set(:,:,1)=Threshold;
    
    ClusterSize_OneTailed_NN26Set_AllMeasures(:,:,iMeasure)=ClusterSize_OneTailed_NN26Set;
    ClusterSize_TwoTailed_NN26Set_AllMeasures(:,:,iMeasure)=ClusterSize_TwoTailed_NN26Set;
    
    GSRClusterSize_OneTailed_NN26Set_AllMeasures(:,:,iMeasure)=GSRClusterSize_OneTailed_NN26Set;
    GSRClusterSize_TwoTailed_NN26Set_AllMeasures(:,:,iMeasure)=GSRClusterSize_TwoTailed_NN26Set;
end

save([OutDir_AlphaSim,'/ClusterSize_AFNI3dClustSim.mat'],'ClusterSize_OneTailed_NN26Set_AllMeasures','ClusterSize_TwoTailed_NN26Set_AllMeasures','GSRClusterSize_OneTailed_NN26Set_AllMeasures','GSRClusterSize_TwoTailed_NN26Set_AllMeasures')






%Apply AFNI 3dClustSim Correction
load([OutDir_AlphaSim,'/ClusterSize_AFNI3dClustSim.mat'])

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

IsTwoTailed=1; %!!!
VoxelPThresholdSet_OneTailed=[0.01 0.005 0.001 0.0005 0.01 0.005 0.001 0.0005];
ClusterPThresholdSet_OneTailed=[0.05 0.05 0.05 0.05 0.025 0.025 0.025 0.025];

for iMeasure=1:length(MeasureSet)
    
    for iP=1:length(VoxelPThresholdSet_OneTailed)
        VoxelPThreshold=2*VoxelPThresholdSet_OneTailed(iP); %For y_ApplyAlphaSimThreshold, input two-tailed p values.
        
        iRow=1;
        iColumn=iP;
        if iP>=5
            iRow=2;
            iColumn=iP-4;
        end
        
        CorrectionName=['OneTailed_AFNI3dClustSim_',num2str(VoxelPThresholdSet_OneTailed(iP)),'_',num2str(ClusterPThresholdSet_OneTailed(iP))];
        mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/']);
        mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/']);
        
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedT'];
        OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/ECEOPairedT'];
        ClusterSize_OneTailed_NN26Set = ClusterSize_OneTailed_NN26Set_AllMeasures(:,:,iMeasure);
        ClusterSize=ClusterSize_OneTailed_NN26Set(iRow,iColumn);
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,OutName,MaskFile);

        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedT'];
        OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/ECEOPairedT'];
        ClusterSize_OneTailed_NN26Set = GSRClusterSize_OneTailed_NN26Set_AllMeasures(:,:,iMeasure);
        ClusterSize=ClusterSize_OneTailed_NN26Set(iRow,iColumn);
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,OutName,MaskFile);
    end
end






%GRF Correction

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};



IsTwoTailed=1; %!!!
VoxelPThresholdSet_OneTailed=[0.01 0.005 0.001 0.0005 0.01 0.005 0.001 0.0005];
ClusterPThresholdSet_OneTailed=[0.05 0.05 0.05 0.05 0.025 0.025 0.025 0.025];


for iMeasure=1:length(MeasureSet)
    
    for iP=1:length(VoxelPThresholdSet_OneTailed)
        VoxelPThreshold=2*VoxelPThresholdSet_OneTailed(iP); %For y_GRF_Threshold, input two-tailed p values.
        ClusterPThreshold=2*ClusterPThresholdSet_OneTailed(iP); %For y_GRF_Threshold, input two-tailed p values.
        
        CorrectionName=['OneTailedGRF_',num2str(VoxelPThresholdSet_OneTailed(iP)),'_',num2str(ClusterPThresholdSet_OneTailed(iP))];
        
        mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/']);
        mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/']);
        
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedT'];
        OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/ECEOPairedT'];
        y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutName,MaskFile);
        
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedT'];
        OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/ECEOPairedT'];
        y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutName,MaskFile);
    end
end





%FDR Correction
qThreshold=0.05;

for iMeasure=1:length(MeasureSet)
    
    CorrectionName='FDR';
    
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/']);
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/']);
    
    FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedT'];
    OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/ECEOPairedT'];
    y_FDR_Image(FileName,qThreshold,OutName,MaskFile);
    
    FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedT'];
    OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/ECEOPairedT'];
    y_FDR_Image(FileName,qThreshold,OutName,MaskFile);
    
end




%Test Overlap for Cluster-based method
OverlapDirUp=[OutDir,'/Binarized'];
PrefixSet={'OneTailedGRF_','OneTailed_DPABIAlphaSim_','OneTailed_AFNI3dClustSim_'};
for iPrefix=1:length(PrefixSet)
    for iP=1:length(VoxelPThresholdSet_OneTailed)
        CorrectionName=[PrefixSet{iPrefix},num2str(VoxelPThresholdSet_OneTailed(iP)),'_',num2str(ClusterPThresholdSet_OneTailed(iP))];
        OverlapDir=[OverlapDirUp,'/',CorrectionName];
        mkdir(OverlapDir)
        for iMeasure=1:length(MeasureSet)
            DataSum=0;
            DataGSRSum=0;
            
            FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/ClusterThresholded_ECEOPairedT'];
            [Data Head]=y_Read(FileName);
            DataSum=DataSum+(Data~=0);
            
            FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/ClusterThresholded_ECEOPairedT'];
            [Data Head]=y_Read(FileName);
            DataGSRSum=DataGSRSum+(Data~=0);
            
            y_Write(DataSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
            y_Write(DataGSRSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
        end
    end
end


%Test Overlap for FDR
CorrectionName='FDR';
OverlapDir=[OverlapDirUp,'/',CorrectionName];
mkdir(OverlapDir)
for iMeasure=1:length(MeasureSet)
    DataSum=0;
    DataGSRSum=0;
    for iSession=1:1
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/ECEOPairedT'];
        [Data Head]=y_Read(FileName);
        DataSum=DataSum+(Data~=0);
        
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/ECEOPairedT'];
        [Data Head]=y_Read(FileName);
        DataGSRSum=DataGSRSum+(Data~=0);
    end
    y_Write(DataSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
    y_Write(DataGSRSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
end




%Test Overlap for PALM
for iMeasure=1:length(MeasureSet)
    OverlapDir=[OverlapDirUp,'/','PALMVox'];
    mkdir(OverlapDir)
    DataSum=0;
    DataGSRSum=0;
    for iSession=1:1
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedTPALM23_vox_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataSum=DataSum+(Data<=0.05);
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedTPALM23_vox_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataGSRSum=DataGSRSum+(Data<=0.05);
    end
    y_Write(DataSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
    y_Write(DataGSRSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
    
    
    OverlapDir=[OverlapDirUp,'/','PALMTFCE'];
    mkdir(OverlapDir)
    DataSum=0;
    DataGSRSum=0;
    for iSession=1:1
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedTPALM23_tfce_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataSum=DataSum+(Data<=0.05);
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedTPALM23_tfce_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataGSRSum=DataGSRSum+(Data<=0.05);
    end
    y_Write(DataSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
    y_Write(DataGSRSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
    
    OverlapDir=[OverlapDirUp,'/','PALMCluster23'];
    mkdir(OverlapDir)
    DataSum=0;
    DataGSRSum=0;
    for iSession=1:1
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedTPALM23_clustere_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataSum=DataSum+(Data<=0.05);
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedTPALM23_clustere_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataGSRSum=DataGSRSum+(Data<=0.05);
    end
    y_Write(DataSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
    y_Write(DataGSRSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
    
    OverlapDir=[OverlapDirUp,'/','PALMCluster258'];
    mkdir(OverlapDir)
    DataSum=0;
    DataGSRSum=0;
    for iSession=1:1
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedTPALM258_clustere_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataSum=DataSum+(Data<=0.05);
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedTPALM258_clustere_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataGSRSum=DataGSRSum+(Data<=0.05);
    end
    y_Write(DataSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
    y_Write(DataGSRSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
    
    OverlapDir=[OverlapDirUp,'/','PALMCluster31'];
    mkdir(OverlapDir)
    DataSum=0;
    DataGSRSum=0;
    for iSession=1:1
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedTPALM31_clustere_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataSum=DataSum+(Data<=0.05);
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedTPALM31_clustere_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataGSRSum=DataGSRSum+(Data<=0.05);
    end
    y_Write(DataSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
    y_Write(DataGSRSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
    
    OverlapDir=[OverlapDirUp,'/','PALMCluster329'];
    mkdir(OverlapDir)
    DataSum=0;
    DataGSRSum=0;
    for iSession=1:1
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/ECEOPairedTPALM329_clustere_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataSum=DataSum+(Data<=0.05);
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/ECEOPairedTPALM329_clustere_tstat_fwep.nii'];
        [Data Head]=y_Read(FileName);
        Data(find(Data==0))=999;
        DataGSRSum=DataGSRSum+(Data<=0.05);
    end
    y_Write(DataSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
    y_Write(DataGSRSum,Head,[OverlapDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
    
end


%%%As a function: END



%Get Reproducibility
%CorrectionSet={'OneTailedGRF_0.01_0.05','OneTailedGRF_0.005_0.05','OneTailedGRF_0.001_0.05','OneTailedGRF_0.0005_0.05','OneTailedGRF_0.01_0.025','OneTailedGRF_0.005_0.025','OneTailedGRF_0.001_0.025','OneTailedGRF_0.0005_0.025','PALMCluster23','PALMCluster258','PALMCluster31','PALMCluster329','PALMTFCE','PALMVox','FDR'};

CorrectionSet={'OneTailed_AFNI3dClustSim_0.01_0.05','OneTailed_AFNI3dClustSim_0.005_0.05','OneTailed_AFNI3dClustSim_0.001_0.05','OneTailed_AFNI3dClustSim_0.0005_0.05','OneTailed_AFNI3dClustSim_0.01_0.025','OneTailed_AFNI3dClustSim_0.005_0.025','OneTailed_AFNI3dClustSim_0.001_0.025','OneTailed_AFNI3dClustSim_0.0005_0.025','OneTailed_DPABIAlphaSim_0.01_0.05','OneTailed_DPABIAlphaSim_0.005_0.05','OneTailed_DPABIAlphaSim_0.001_0.05','OneTailed_DPABIAlphaSim_0.0005_0.05','OneTailed_DPABIAlphaSim_0.01_0.025','OneTailed_DPABIAlphaSim_0.005_0.025','OneTailed_DPABIAlphaSim_0.001_0.025','OneTailed_DPABIAlphaSim_0.0005_0.025','OneTailedGRF_0.01_0.05','OneTailedGRF_0.005_0.05','OneTailedGRF_0.001_0.05','OneTailedGRF_0.0005_0.05','OneTailedGRF_0.01_0.025','OneTailedGRF_0.005_0.025','OneTailedGRF_0.001_0.025','OneTailedGRF_0.0005_0.025','PALMCluster23','PALMCluster258','PALMCluster31','PALMCluster329','PALMTFCE','PALMVox','FDR'};

DataDirUp='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Zou_Beijing/ECEO_PairedT/ECEO_PairedT8mm/Binarized';
DataDirUp2='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Dong_Beijing/ECEO_PairedT/ECEO_PairedT_8mm/Binarized';

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

save(['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Reproducibility/Revision_8mm_Dice_ZouBeijing_DongBeijing.mat'],'ReproducibilitySet','JaccardSet','DiceSet','DiceSet_Selected');

[p tbl stats]=friedman(DiceSet)
[c m]=multcompare(stats)









%Get the number of significant voxels
CorrectionSet={'OneTailed_AFNI3dClustSim_0.01_0.05','OneTailed_AFNI3dClustSim_0.005_0.05','OneTailed_AFNI3dClustSim_0.001_0.05','OneTailed_AFNI3dClustSim_0.0005_0.05','OneTailed_AFNI3dClustSim_0.01_0.025','OneTailed_AFNI3dClustSim_0.005_0.025','OneTailed_AFNI3dClustSim_0.001_0.025','OneTailed_AFNI3dClustSim_0.0005_0.025','OneTailed_DPABIAlphaSim_0.01_0.05','OneTailed_DPABIAlphaSim_0.005_0.05','OneTailed_DPABIAlphaSim_0.001_0.05','OneTailed_DPABIAlphaSim_0.0005_0.05','OneTailed_DPABIAlphaSim_0.01_0.025','OneTailed_DPABIAlphaSim_0.005_0.025','OneTailed_DPABIAlphaSim_0.001_0.025','OneTailed_DPABIAlphaSim_0.0005_0.025','OneTailedGRF_0.01_0.05','OneTailedGRF_0.005_0.05','OneTailedGRF_0.001_0.05','OneTailedGRF_0.0005_0.05','OneTailedGRF_0.01_0.025','OneTailedGRF_0.005_0.025','OneTailedGRF_0.001_0.025','OneTailedGRF_0.0005_0.025','PALMCluster23','PALMCluster258','PALMCluster31','PALMCluster329','PALMTFCE','PALMVox','FDR'};
DataDirUp='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Zou_Beijing/ECEO_PairedT/ECEO_PairedT8mm/Binarized';
DataDirUp2='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/EOEC/Dong_Beijing/ECEO_PairedT/ECEO_PairedT_8mm/Binarized';

SignificantVoxels1Set=zeros(10,length(CorrectionSet));
SignificantVoxels2Set=zeros(10,length(CorrectionSet));
for iCorrection=1:length(CorrectionSet)
    DataDir=[DataDirUp,'/',CorrectionSet{iCorrection}];
    DataDir2=[DataDirUp2,'/',CorrectionSet{iCorrection}];
    
    MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
    ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
    ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
    SignificantVoxels1=[];
    SignificantVoxels1GSR=[];
    SignificantVoxels2=[];
    SignificantVoxels2GSR=[];
    for iMeasure=1:length(MeasureSet)
        File=[DataDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'.nii'];
        Data=y_ReadRPI(File);
        File=[DataDir2,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'.nii'];
        Data2=y_ReadRPI(File);
        SignificantVoxels1(iMeasure,1)=length(find(Data));
        SignificantVoxels2(iMeasure,1)=length(find(Data2));

        File=[DataDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'.nii'];
        Data=y_ReadRPI(File);
        File=[DataDir2,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'.nii'];
        Data2=y_ReadRPI(File);
        SignificantVoxels1GSR(iMeasure,1)=length(find(Data));
        SignificantVoxels2GSR(iMeasure,1)=length(find(Data2));
    end
    
    SignificantVoxels1=[SignificantVoxels1;SignificantVoxels1GSR];
    SignificantVoxels2=[SignificantVoxels2;SignificantVoxels2GSR];
    
    SignificantVoxels1Set(:,iCorrection)=SignificantVoxels1;
    SignificantVoxels2Set(:,iCorrection)=SignificantVoxels2;
end


