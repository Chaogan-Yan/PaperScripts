%To increase readbility, here is an example for ALFF after 8mm FWMH smoothing


%RANDOM TEST
%Group Analysis
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

Site=SiteID_AllSites;
SubID=SubID_AllSites;
MeanFD = MeanFDJenkinson_AllSites;

WantedSubIndex = find((Site==17).*(Sex==-1));
%Select subjects
SubID=SubID(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
MeanFD=MeanFD(WantedSubIndex);
Site=Site(WantedSubIndex);
SiteName_AllSites = SiteName_AllSites(WantedSubIndex);
SiteID_AllSites = SiteID_AllSites(WantedSubIndex);


PALMSettings.nPerm = 1000; %5000;
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

DataDir='/mnt/Data/RfMRILab/RfMRIMaps/Preprocessed/FCPIntermediateFiles/Beijing';
OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/RandPerm8mm';
MaskFile = ['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/SubInfo/FCP_GroupMask90Percent.nii'];

FakeContrast=[ones(20,1);-1*ones(20,1)]; %FakeContrast=[ones(40,1);-1*ones(40,1)];


load('/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/RandPerm/RandPermVectorSet.mat')

parfor iPerm = 1:1000
    RandPermVector=RandPermVectorSet(iPerm,:);
    
    RandName=['Rand',num2str(iPerm)];
    
    mkdir([OutDir,'/',RandName]);
    for iMeasure=1:1
        mkdir([OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
        FileList=[];
        FileListGSR=[];
        for iSub=1:length(FakeContrast)
            FileList{iSub,1}=[DataDir,'/',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{RandPermVector(iSub)},'.nii'];
            FileListGSR{iSub,1}=[DataDir,'/',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{RandPermVector(iSub)},'.nii'];
        end
        AgeFakeSet=Age(RandPermVector((1:length(FakeContrast))));
        MeanFDFakeSet=MeanFD(RandPermVector((1:length(FakeContrast))));
        Cov=[FakeContrast,AgeFakeSet,MeanFDFakeSet,ones(length(FakeContrast),1)];
        Contrast=zeros(1,size(Cov,2));
        Contrast(1)=1;
        
        OutputName=[OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrast'];
        % [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(DependentVolume,Predictor,OutputName,MaskFile,CovVolume,Contrast,TF_Flag,IsOutputResidual,Header)
        [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0);
        
        OutputName=[OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrastPALM23'];
        y_GroupAnalysis_PermutationTest_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings);

        OutputName=[OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrastPALM31'];
        y_GroupAnalysis_PermutationTest_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings31);

        OutputName=[OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrastPALM258'];
        y_GroupAnalysis_PermutationTest_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings258);

        OutputName=[OutDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrastPALM329'];
        y_GroupAnalysis_PermutationTest_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0,[],PALMSettings329);
    end
    
end






%Get smoothness

DataDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/RandPerm8mm';

MeasureSet={'ALFF'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

MaskFile = ['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/SubInfo/FCP_GroupMask90Percent.nii'];
MaskData = y_ReadRPI(MaskFile);
MaskIndex = find(MaskData);
nVoxels = length(MaskIndex);

FWHMAllRandSet=[];
for iMeasure=1:1
    FWHMAllRand=zeros(3,1000);
    for iRand=1:1000
        RandName=['Rand',num2str(iRand)];
        FileName=[DataDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrast.nii'];
        [Data Header]=y_Read(FileName);
        Header=w_ReadDLH(Header);
        FWHMAllRand(:,iRand)=[Header.FWHMx,Header.FWHMy,Header.FWHMz]';

    end
    FWHMAllRandSet(:,:,iMeasure)=FWHMAllRand;
end





%Do DPABI AlphaSim
addpath /mnt/Data/RfMRILab/Yan/YAN_Program/y_MyOwn
MaskFile = ['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/SubInfo/FCP_GroupMask90Percent.nii'];
load /mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/AlphaSim/FWHM.mat

ClusterSize_OneTailed_NN26Set_AllMeasures=zeros(4,4,size(FWHMAllRandSet,2),5);
ClusterSize_TwoTailed_NN26Set_AllMeasures=zeros(4,4,size(FWHMAllRandSet,2),5);

for iMeasure=1:1
    
    FWHMAll=FWHMAllRandSet(:,:,iMeasure);
    
    ClusterSize_OneTailed_NN26Set=zeros(4,4,size(FWHMAllRandSet,2));
    ClusterSize_TwoTailed_NN26Set=zeros(4,4,size(FWHMAllRandSet,2));

    parfor iRand=1:size(FWHMAllRandSet,2)
        disp(iRand)
        ClusterSize_OneTailed_NN26Temp=[];
        ClusterSize_TwoTailed_NN26Temp=[];
        [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.01,1000,'fwhm',FWHMAll(:,iRand)');
        ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
        ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
        
        [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.005,1000,'fwhm',FWHMAll(:,iRand)');
        ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
        ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
        
        [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.001,1000,'fwhm',FWHMAll(:,iRand)');
        ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
        ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
        
        [ClusterSize_OneTailed_NN6 ClusterSize_OneTailed_NN18 ClusterSize_OneTailed_NN26 ClusterSize_TwoTailed_NN6 ClusterSize_TwoTailed_NN18 ClusterSize_TwoTailed_NN26]=y_AlphaSim_Threshold(MaskFile,[],[],0.0005,1000,'fwhm',FWHMAll(:,iRand)');
        ClusterSize_OneTailed_NN26Temp=[ClusterSize_OneTailed_NN26Temp,ClusterSize_OneTailed_NN26];
        ClusterSize_TwoTailed_NN26Temp=[ClusterSize_TwoTailed_NN26Temp,ClusterSize_TwoTailed_NN26];
        
        ClusterSize_OneTailed_NN26Set(:,:,iRand)=ClusterSize_OneTailed_NN26Temp;
        ClusterSize_TwoTailed_NN26Set(:,:,iRand)=ClusterSize_TwoTailed_NN26Temp;

    end
    ClusterSize_OneTailed_NN26Set_AllMeasures(:,:,:,iMeasure)=ClusterSize_OneTailed_NN26Set;
    ClusterSize_TwoTailed_NN26Set_AllMeasures(:,:,:,iMeasure)=ClusterSize_TwoTailed_NN26Set;
end




%Get FWER for DPABI AlphaSim
DataDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/RandPerm8mm';
TempDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/AlphaSim/ParForOK/Temp';

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

MaskFile = ['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/SubInfo/FCP_GroupMask90Percent.nii'];
MaskData = y_ReadRPI(MaskFile);
MaskIndex = find(MaskData);
nVoxels = length(MaskIndex);

load /mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/AlphaSim/ClusterSize_DPABIAlphaSim.mat

FDRateSetAllRand_AllMeasures=zeros(8,1000,length(MeasureSet));
ClusterSizeSetAllRand_AllMeasures=zeros(8,1000,length(MeasureSet));
for iMeasure=1:1 %length(MeasureSet)
    FDRateSetAllRand=zeros(8,1000);
    ClusterSizeSetAllRand=zeros(8,1000);
    parfor iRand=1:1000
        RandName=['Rand',num2str(iRand)];
        FDRateSet=zeros(8,1);
        ClusterSizeSet=zeros(8,1);
        
        FileName=[DataDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrast.nii'];
        ClusterSize_OneTailed_NN26Set = ClusterSize_OneTailed_NN26Set_AllMeasures(:,:,:,iMeasure);
        ClusterSize=ClusterSize_OneTailed_NN26Set(1,1,iRand);
        VoxelPThreshold=0.02;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(1,1)=FDRate;
        ClusterSizeSet(1,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(1,2,iRand);
        VoxelPThreshold=0.01;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(2,1)=FDRate;
        ClusterSizeSet(2,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(1,3,iRand);
        VoxelPThreshold=0.002;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(3,1)=FDRate;
        ClusterSizeSet(3,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(1,4,iRand);
        VoxelPThreshold=0.001;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(4,1)=FDRate;
        ClusterSizeSet(4,1)=ClusterSize;
        
        
        %The each tail p<0.025
        ClusterSize=ClusterSize_OneTailed_NN26Set(2,1,iRand);
        VoxelPThreshold=0.02;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(5,1)=FDRate;
        ClusterSizeSet(5,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(2,2,iRand);
        VoxelPThreshold=0.01;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(6,1)=FDRate;
        ClusterSizeSet(6,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(2,3,iRand);
        VoxelPThreshold=0.002;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(7,1)=FDRate;
        ClusterSizeSet(7,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(2,4,iRand);
        VoxelPThreshold=0.001;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(8,1)=FDRate;
        ClusterSizeSet(8,1)=ClusterSize;
        
        
        FDRateSetAllRand(:,iRand)=FDRateSet;
        ClusterSizeSetAllRand(:,iRand)=ClusterSizeSet;
        
    end
    
    FDRateSetAllRand_AllMeasures(:,:,iMeasure)=FDRateSetAllRand;
    ClusterSizeSetAllRand_AllMeasures(:,:,iMeasure)=ClusterSizeSetAllRand;
end

for iMeasure=1:length(MeasureSet)
    FWERate = sum(FDRateSetAllRand_AllMeasures~=0,2)/size(FDRateSetAllRand_AllMeasures,2);
end

save('AllMeasures_FWERate_ClusterSize_DPABIAlphaSim.mat','FDRateSetAllRand_AllMeasures','ClusterSizeSetAllRand_AllMeasures','FWERate');






%For AFNI

%save FWHM files
load /mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/AlphaSim/FWHM.mat
FWHMAll=FWHMAllRandSet(:,:,1)';
save(['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/AlphaSim/FWHM_1.txt'], 'FWHMAll', '-ASCII', '-DOUBLE','-TABS')


% Then perform on Cluster
#! /bin/bash
FWHM_list=/mnt/diske/Yan/Analysis/MultipleComparison/FCPData/Beijing/AlphaSim/Revision/FWHM_1.txt
i=1
cat ${FWHM_list} | while read FWHM
do
	echo "i ${i}"
	cmd=` echo 3dClustSim -mask /mnt/diske/Yan/Analysis/MultipleComparison/FCPData/FCP_GroupMask90Percent.nii -pthr 0.01 0.005 0.001 0.0005 -athr 0.05 0.025 -iter 1000 -fwhmxyz ${FWHM} -prefix /mnt/diske/Yan/Analysis/MultipleComparison/FCPData/Beijing/AlphaSim/Revision/AFNI3dClusterSim/FWHM_1_Rand${i}`
	echo "${cmd}" |  qsub -N FWHM_1-${i} -m n -l nodes=1:ppn=1 -q wanzhao
	sleep 0.1
	let i=i+1
done





%Get the thresholds

load /mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/AlphaSim/FWHM.mat
ClusterSize_OneTailed_NN26Set_AllMeasures=zeros(2,4,size(FWHMAllRandSet,2),1);
ClusterSize_TwoTailed_NN26Set_AllMeasures=zeros(2,4,size(FWHMAllRandSet,2),1);

DataDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/AlphaSim/AFNI3dClusterSim';
for iMeasure=1:1
    ClusterSize_OneTailed_NN26Set=zeros(2,4,size(FWHMAllRandSet,2));
    ClusterSize_TwoTailed_NN26Set=zeros(2,4,size(FWHMAllRandSet,2));

    for iRand=1:1000
        File=[DataDir,'/FWHM_',num2str(iMeasure),'_Rand',num2str(iRand),'.NN3_1sided.1D'];
        fid = fopen(File);
        StringFilter = '%f%f%f\n';
        %StringFilter = [StringFilter,'%*[^\n]']; %Skip the else till end of the line
        for i=1:8
            tline = fgetl(fid); %Skip the title line
        end
        ThresholdTemp = textscan(fid,StringFilter);
        fclose(fid);
        Threshold=[ThresholdTemp{2},ThresholdTemp{3}]';
        ClusterSize_OneTailed_NN26Set(:,:,iRand)=Threshold;
    end
    
    for iRand=1:1000
        File=[DataDir,'/FWHM_',num2str(iMeasure),'_Rand',num2str(iRand),'.NN3_2sided.1D'];
        fid = fopen(File);
        StringFilter = '%f%f%f\n';
        %StringFilter = [StringFilter,'%*[^\n]']; %Skip the else till end of the line
        for i=1:8
            tline = fgetl(fid); %Skip the title line
        end
        ThresholdTemp = textscan(fid,StringFilter);
        fclose(fid);
        Threshold=[ThresholdTemp{2},ThresholdTemp{3}]';
        ClusterSize_TwoTailed_NN26Set(:,:,iRand)=Threshold;
    end
    
    ClusterSize_OneTailed_NN26Set_AllMeasures(:,:,:,iMeasure)=ClusterSize_OneTailed_NN26Set;
    ClusterSize_TwoTailed_NN26Set_AllMeasures(:,:,:,iMeasure)=ClusterSize_TwoTailed_NN26Set;
end






%Get FWER for AFNI 3dClusterSim

DataDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/RandPerm8mm';

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

MaskFile = ['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/SubInfo/FCP_GroupMask90Percent.nii'];
MaskData = y_ReadRPI(MaskFile);
MaskIndex = find(MaskData);
nVoxels = length(MaskIndex);

load /mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/AlphaSim/ClusterSize_AFNI3dClustSim.mat

FDRateSetAllRand_AllMeasures=zeros(8,1000,length(MeasureSet));
ClusterSizeSetAllRand_AllMeasures=zeros(8,1000,length(MeasureSet));
for iMeasure=1:1
    FDRateSetAllRand=zeros(8,1000);
    ClusterSizeSetAllRand=zeros(8,1000);
    parfor iRand=1:1000
        RandName=['Rand',num2str(iRand)];
        FDRateSet=zeros(8,1);
        ClusterSizeSet=zeros(8,1);
        
        FileName=[DataDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrast.nii'];
        ClusterSize_OneTailed_NN26Set = ClusterSize_OneTailed_NN26Set_AllMeasures(:,:,:,iMeasure);
        ClusterSize=ClusterSize_OneTailed_NN26Set(1,1,iRand);
        VoxelPThreshold=0.02;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(1,1)=FDRate;
        ClusterSizeSet(1,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(1,2,iRand);
        VoxelPThreshold=0.01;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(2,1)=FDRate;
        ClusterSizeSet(2,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(1,3,iRand);
        VoxelPThreshold=0.002;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(3,1)=FDRate;
        ClusterSizeSet(3,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(1,4,iRand);
        VoxelPThreshold=0.001;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(4,1)=FDRate;
        ClusterSizeSet(4,1)=ClusterSize;
        
        
        %The each tail p<0.025
        ClusterSize=ClusterSize_OneTailed_NN26Set(2,1,iRand);
        VoxelPThreshold=0.02;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(5,1)=FDRate;
        ClusterSizeSet(5,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(2,2,iRand);
        VoxelPThreshold=0.01;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(6,1)=FDRate;
        ClusterSizeSet(6,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(2,3,iRand);
        VoxelPThreshold=0.002;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(7,1)=FDRate;
        ClusterSizeSet(7,1)=ClusterSize;
        
        ClusterSize=ClusterSize_OneTailed_NN26Set(2,4,iRand);
        VoxelPThreshold=0.001;
        IsTwoTailed=1;
        %[Data_Corrected, ClusterSize, Header]=y_ApplyAlphaSimThreshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterSize,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header)
        [Data_Corrected ClusterSize] = y_ApplyAlphaSimThreshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterSize,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(8,1)=FDRate;
        ClusterSizeSet(8,1)=ClusterSize;
        
        
        FDRateSetAllRand(:,iRand)=FDRateSet;
        ClusterSizeSetAllRand(:,iRand)=ClusterSizeSet;
    end
    
    FDRateSetAllRand_AllMeasures(:,:,iMeasure)=FDRateSetAllRand;
    ClusterSizeSetAllRand_AllMeasures(:,:,iMeasure)=ClusterSizeSetAllRand;
end


for iMeasure=1:length(MeasureSet)
    FWERate = sum(FDRateSetAllRand_AllMeasures~=0,2)/size(FDRateSetAllRand_AllMeasures,2);
end

save('AllMeasures_FWERate_ClusterSize_AFNI3dClusterSim.mat','FDRateSetAllRand_AllMeasures','ClusterSizeSetAllRand_AllMeasures','FWERate');







%Get FWE for GRF

DataDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/RandPerm8mm';

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

MaskFile = ['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/SubInfo/FCP_GroupMask90Percent.nii'];
MaskData = y_ReadRPI(MaskFile);
MaskIndex = find(MaskData);
nVoxels = length(MaskIndex);


FDRateSetAllRand_AllMeasures=zeros(8,1000,length(MeasureSet));
ClusterSizeSetAllRand_AllMeasures=zeros(8,1000,length(MeasureSet));
GSRFDRateSetAllRand_AllMeasures=zeros(8,1000,length(MeasureSet));
GSRClusterSizeSetAllRand_AllMeasures=zeros(8,1000,length(MeasureSet));
for iMeasure=1:1 %length(MeasureSet)
    FDRateSetAllRand=zeros(8,1000);
    ClusterSizeSetAllRand=zeros(8,1000);
    GSRFDRateSetAllRand=zeros(8,1000);
    GSRClusterSizeSetAllRand=zeros(8,1000);
    parfor iRand=1:1000
        RandName=['Rand',num2str(iRand)];
        FDRateSet=zeros(8,1);
        ClusterSizeSet=zeros(8,1);
        
        FileName=[DataDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrast.nii'];
        
        VoxelPThreshold=0.02;
        IsTwoTailed=1;
        ClusterPThreshold=0.1;
        %[Data_Corrected, ClusterSize, Header]=y_GRF_Threshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header, dLh)
        [Data_Corrected ClusterSize] = y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(1,1)=FDRate;
        ClusterSizeSet(1,1)=ClusterSize;
        
        VoxelPThreshold=0.01;
        IsTwoTailed=1;
        ClusterPThreshold=0.1;
        %[Data_Corrected, ClusterSize, Header]=y_GRF_Threshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header, dLh)
        [Data_Corrected ClusterSize] = y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(2,1)=FDRate;
        ClusterSizeSet(2,1)=ClusterSize;
        
        VoxelPThreshold=0.002;
        IsTwoTailed=1;
        ClusterPThreshold=0.1;
        %[Data_Corrected, ClusterSize, Header]=y_GRF_Threshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header, dLh)
        [Data_Corrected ClusterSize] = y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(3,1)=FDRate;
        ClusterSizeSet(3,1)=ClusterSize;
        
        VoxelPThreshold=0.001;
        IsTwoTailed=1;
        ClusterPThreshold=0.1;
        %[Data_Corrected, ClusterSize, Header]=y_GRF_Threshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header, dLh)
        [Data_Corrected ClusterSize] = y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(4,1)=FDRate;
        ClusterSizeSet(4,1)=ClusterSize;
        
        
        %Then each tail p<0.025
        VoxelPThreshold=0.02;
        IsTwoTailed=1;
        ClusterPThreshold=0.05;
        %[Data_Corrected, ClusterSize, Header]=y_GRF_Threshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header, dLh)
        [Data_Corrected ClusterSize] = y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(5,1)=FDRate;
        ClusterSizeSet(5,1)=ClusterSize;
        
        VoxelPThreshold=0.01;
        IsTwoTailed=1;
        ClusterPThreshold=0.05;
        %[Data_Corrected, ClusterSize, Header]=y_GRF_Threshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header, dLh)
        [Data_Corrected ClusterSize] = y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(6,1)=FDRate;
        ClusterSizeSet(6,1)=ClusterSize;
        
        VoxelPThreshold=0.002;
        IsTwoTailed=1;
        ClusterPThreshold=0.05;
        %[Data_Corrected, ClusterSize, Header]=y_GRF_Threshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header, dLh)
        [Data_Corrected ClusterSize] = y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(7,1)=FDRate;
        ClusterSizeSet(7,1)=ClusterSize;
        
        VoxelPThreshold=0.001;
        IsTwoTailed=1;
        ClusterPThreshold=0.05;
        %[Data_Corrected, ClusterSize, Header]=y_GRF_Threshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header, dLh)
        [Data_Corrected ClusterSize] = y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(8,1)=FDRate;
        ClusterSizeSet(8,1)=ClusterSize;

        FDRateSetAllRand(:,iRand)=FDRateSet;
        ClusterSizeSetAllRand(:,iRand)=ClusterSizeSet;
        
        
    end
    FDRateSetAllRand_AllMeasures(:,:,iMeasure)=FDRateSetAllRand;
    ClusterSizeSetAllRand_AllMeasures(:,:,iMeasure)=ClusterSizeSetAllRand;
end

for iMeasure=1:length(MeasureSet)
    FWERate = sum(FDRateSetAllRand_AllMeasures~=0,2)/size(FDRateSetAllRand_AllMeasures,2);
end

save('AllMeasures_FWERate_ClusterSize_GRF.mat','FDRateSetAllRand_AllMeasures','ClusterSizeSetAllRand_AllMeasures','FWERate');







%Get PALM and FDR

DataDir='/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FWER/RandPerm8mm';

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

MaskFile = ['/mnt/Data/RfMRILab/Yan/YAN_Work/MultipleComparison/MaleFemale/FCP/SubInfo/FCP_GroupMask90Percent.nii'];
MaskData = y_ReadRPI(MaskFile);
MaskIndex = find(MaskData);
nVoxels = length(MaskIndex);

FDRateSetAllRand_AllMeasures=zeros(7,1000,length(MeasureSet));
for iMeasure=1:1 %length(MeasureSet)
    FDRateSetAllRand=zeros(7,1000);
    parfor iRand=1:1000
        RandName=['Rand',num2str(iRand)];
        
        
        FDRateSet=zeros(7,1);
        FileName=[DataDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrastPALM23_clustere_tstat_fwep.nii'];
        Data=y_Read(FileName);
        Data(Data==0)=1;
        FDRate = length(find(Data<0.05))/nVoxels;
        FDRateSet(1,1)=FDRate;
        
        FileName=[DataDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrastPALM258_clustere_tstat_fwep.nii'];
        Data=y_Read(FileName);
        Data(Data==0)=1;
        FDRate = length(find(Data<0.05))/nVoxels;
        FDRateSet(2,1)=FDRate;
        
        FileName=[DataDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrastPALM31_clustere_tstat_fwep.nii'];
        Data=y_Read(FileName);
        Data(Data==0)=1;
        FDRate = length(find(Data<0.05))/nVoxels;
        FDRateSet(3,1)=FDRate;
        
        FileName=[DataDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrastPALM329_clustere_tstat_fwep.nii'];
        Data=y_Read(FileName);
        Data(Data==0)=1;
        FDRate = length(find(Data<0.05))/nVoxels;
        FDRateSet(4,1)=FDRate;
        
        FileName=[DataDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrastPALM23_tfce_tstat_fwep.nii'];
        Data=y_Read(FileName);
        Data(Data==0)=1;
        FDRate = length(find(Data<0.05))/nVoxels;
        FDRateSet(5,1)=FDRate;
        
        
        FileName=[DataDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrastPALM23_vox_tstat_fwep.nii'];
        Data=y_Read(FileName);
        Data(Data==0)=1;
        FDRate = length(find(Data<0.05))/nVoxels;
        FDRateSet(6,1)=FDRate;
        
        
        %FDR
        FileName=[DataDir,'/',RandName,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/FakeContrast.nii'];
        [Data_Corrected, Header]=y_FDR_Image(FileName,0.05,[TempDir,'/',RandName],MaskFile);
        FDRate = length(find(Data_Corrected~=0))/nVoxels;
        FDRateSet(7,1)=FDRate;

        
        FDRateSetAllRand(:,iRand)=FDRateSet;
    end
    
    FDRateSetAllRand_AllMeasures(:,:,iMeasure)=FDRateSetAllRand;
end


for iMeasure=1:length(MeasureSet)
    FWERate = sum(FDRateSetAllRand_AllMeasures~=0,2)/size(FDRateSetAllRand_AllMeasures,2);
end

save('AllMeasures_FWERate_PALM_FDR.mat','FDRateSetAllRand_AllMeasures','FWERate');








