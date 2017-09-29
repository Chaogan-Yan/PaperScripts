% %Generate RandPerms
% for i=1:100
%     RandPermSet_Male116(i,:)=randperm(116);
%     RandPermSet_FeMale105(i,:)=randperm(105);
% end
% save('/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/SampleSizeAnalysis/MaleVsFemale/RandPermSet.mat','RandPermSet_Male116','RandPermSet_FeMale105')


load('/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/SubInfo/SubInfo_420.mat')


Sex_Site8=Sex(find(Site==8));
Age_Site8=Age(find(Site==8));
Motion_Site8=Motion(find(Site==8),:);
SubID_Site8=SubID(find(Site==8));

MaleIndex=find(Sex_Site8==1); %116
FemaleIndex=find(Sex_Site8==-1); %105

load('/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/SampleSizeAnalysis/MaleVsFemale/RandPermSet.mat')



PALMSettings.nPerm = 1000;
PALMSettings.ClusterInference=0;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'


SessionSet={'','S2_'};

MeasureSet={'ALFF'};
%MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

DataDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/ReAnalysis_FWMH8mm/AllMaps_8mm';
OutDirUp='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/SampleSizeAnalysis/MaleVsFemale/MaleVsFemaleT_Site8_8mmReAnalysis';

MaskFile = '/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/SubInfo/GroupMask_90percent_429AfterExcluding.nii';


SampleSizeSet=[15:5:50,60:10:100];

for iSampleSize=1:length(SampleSizeSet)
    parfor iRand=1:100
        Sex=Sex_Site8([MaleIndex(RandPermSet_Male116(iRand,1:SampleSizeSet(iSampleSize))) ; FemaleIndex(RandPermSet_FeMale105(iRand,1:SampleSizeSet(iSampleSize)))]);
        Age=Age_Site8([MaleIndex(RandPermSet_Male116(iRand,1:SampleSizeSet(iSampleSize))) ; FemaleIndex(RandPermSet_FeMale105(iRand,1:SampleSizeSet(iSampleSize)))]);
        Motion=Motion_Site8([MaleIndex(RandPermSet_Male116(iRand,1:SampleSizeSet(iSampleSize))) ; FemaleIndex(RandPermSet_FeMale105(iRand,1:SampleSizeSet(iSampleSize)))],:);
        SubID=SubID_Site8([MaleIndex(RandPermSet_Male116(iRand,1:SampleSizeSet(iSampleSize))) ; FemaleIndex(RandPermSet_FeMale105(iRand,1:SampleSizeSet(iSampleSize)))]);
        
        OutDir=[OutDirUp,'/',num2str(SampleSizeSet(iSampleSize)),'/',num2str(iRand)];
        
        for iSession=1:length(SessionSet)
            mkdir([OutDir,'/',SessionSet{iSession},'MaleVsFemaleT']);
            for iMeasure=1:length(MeasureSet)
                mkdir([OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
                mkdir([OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
                FileList=[];
                FileListGSR=[];
                
                for iSub=1:length(SubID)
                    FileList{iSub,1}=[DataDir,'/',SessionSet{iSession},ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
                    FileListGSR{iSub,1}=[DataDir,'/',SessionSet{iSession},ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
                end
                Cov=[Sex,Age,Motion(:,iSession),ones(length(SubID),1)];
                Contrast=zeros(1,size(Cov,2));
                Contrast(1)=1;
                
                OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleTPALMTFCE'];
                y_GroupAnalysis_PermutationTest_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings);
                %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/Temp'])
                %         OutputName=[OutDir,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MaleVsFemaleTPALMTFCE'];
                %         y_GroupAnalysis_PermutationTest_Image(FileListGSR,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings);
                %         %eval(['!rm -rf ',OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/Temp'])
            end
        end
        
    end
    
    eval(['!rm -rf ',OutDirUp,'/',num2str(SampleSizeSet(iSampleSize)),'/*/*/*/Temp'])
    
end








%Get number of reproducible voxels
SampleSizeSet=[15:5:50,60:10:100];
CorrectionSet={'PALMTFCE'};

DataDirUp='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/ReAnalysis_FWMH8mm/FCP716MaleVsFemale/Binarized';
DataDirUp2='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/ReAnalysis_FWMH8mm/MaleVsFemaleT_420/MaleVsFemaleT/Binarized';
DataDirUp3='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/ReAnalysis_FWMH8mm/MaleVsFemaleT_420/S2_MaleVsFemaleT/Binarized';



DataDirUp_SampleSizeTest='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/SampleSizeAnalysis/MaleVsFemale/MaleVsFemaleT_Site8_8mmReAnalysis';
MaskFile = '/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/CORR/BetweenSession/SubInfo/GroupMask_90percent_429AfterExcluding.nii';

SessionSet={'','S2_'};

iCorrection=1;

MeasureSet={'ALFF'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};

VoxelsOverlapped=[];
VoxelsCCBD=[];
VoxelsCore=[];

iMeasure=1
File=[DataDirUp,'/',CorrectionSet{iCorrection},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'.nii'];
DataFCP=y_ReadRPI(File);
File=[DataDirUp2,'/',CorrectionSet{iCorrection},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'.nii'];
Data_BetweenSession1=y_ReadRPI(File);
File=[DataDirUp3,'/',CorrectionSet{iCorrection},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'.nii'];
Data_BetweenSession2=y_ReadRPI(File);
VoxelsCore=length(find(DataFCP+Data_BetweenSession1+Data_BetweenSession2==3));

for iSampleSize=1:length(SampleSizeSet)
    for iRand=1:100
        DataDir_SampleSizeTest=[DataDirUp_SampleSizeTest,'/',num2str(SampleSizeSet(iSampleSize)),'/',num2str(iRand)];

        DataSum=0;
        for iSession=1:length(SessionSet)
            FileName=[DataDir_SampleSizeTest,'/',SessionSet{iSession},'MaleVsFemaleT/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MaleVsFemaleTPALMTFCE_tfce_tstat_fwep.nii'];
            [Data Head]=y_Read(FileName);
            Data(find(Data==0))=999;
            DataSig=Data<=0.05;
            DataSum=DataSum+(DataSig);
            
            VoxelsSampleSizeTest(iSampleSize,iRand,iSession)=length(find(DataSig));
            VoxelsOverlappedWithCore(iSampleSize,iRand,iSession)=length(find( (DataFCP+Data_BetweenSession1+Data_BetweenSession2==3) + (DataSig==1) ==2));
        end
        VoxelsSampleSizeTest(iSampleSize,iRand,3)=length(find(DataSum>=1));
        VoxelsOverlappedWithCore(iSampleSize,iRand,3)=length(find( (DataFCP+Data_BetweenSession1+Data_BetweenSession2==3) + (DataSum>=1) ==2));
        VoxelsSampleSizeTest(iSampleSize,iRand,4)=length(find(DataSum==2));
        VoxelsOverlappedWithCore(iSampleSize,iRand,4)=length(find( (DataFCP+Data_BetweenSession1+Data_BetweenSession2==3) + (DataSum==2) ==2));
        
        Jaccard(iSampleSize,iRand)=length(find(DataSum==2))/length(find(DataSum>=1));
        Dice(iSampleSize,iRand)=2*length(find(DataSum==2))/(length(find(DataSum>=1))+length(find(DataSum==2)));
    end
end


% nVoxels=50480;
% Specificity=(nVoxels-VoxelsSetCCBD-VoxelsSetCore+VoxelsSetOverlapped)./(nVoxels-VoxelsSetCore);
Sensitivity=VoxelsOverlappedWithCore/VoxelsCore;
PositivePredictiveValue=VoxelsOverlappedWithCore./VoxelsSampleSizeTest;

save('/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/Reproducibility/SampleSizeReanalysis/SampleSizeEffect.mat','VoxelsSampleSizeTest','VoxelsOverlappedWithCore','VoxelsSampleSizeTest','VoxelsOverlappedWithCore','VoxelsSampleSizeTest','VoxelsOverlappedWithCore','Jaccard','Dice','Sensitivity','PositivePredictiveValue')
