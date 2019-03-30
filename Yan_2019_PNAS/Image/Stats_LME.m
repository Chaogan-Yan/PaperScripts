

%Select subjects
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_943_846/DataQC/CorrSet.mat
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/Info_Final1789_943_846.mat

ReHoGood = (CorrSet_All(:,3) >= 0.6); %Exclude ReHo Correlation < 0.6

%Exclude Site with N<10
SubjectNumberPerSite=[];
SiteIndex = unique(Site);
WantedSubMatrix=ones(length(SubID),1);
for i=1:length(SiteIndex)
    DxTemp=Dx(find((Site==SiteIndex(i)).*ReHoGood)); %DxTemp=Dx(find(Site==SiteIndex(i)));
    SubjectNumberPerSite(i,:)=[SiteIndex(i),length(find(DxTemp==1)),length(find(DxTemp==-1))];
    if (length(find(DxTemp==1))<10)||(length(find(DxTemp==-1))<10)
        WantedSubMatrix(find(Site==SiteIndex(i)))=0;
    end
end

WantedSubMatrix = WantedSubMatrix.*ReHoGood;

%Select subjects
WantedSubIndex = find(WantedSubMatrix);
SubID=SubID(WantedSubIndex);
Dx=Dx(WantedSubIndex);
Age=Age(WantedSubIndex);
Sex=Sex(WantedSubIndex);
Edu=Edu(WantedSubIndex);
Site=Site(WantedSubIndex);
Motion=Motion(WantedSubIndex,:);


X=[ones(size(Dx)),Dx,Age,Sex,Edu,Motion];
Z={ones(size(Dx)),Dx};
G={Site,Site};


MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};

DataDir='/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/AllMaps';
%OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_943_846/Stats_All';
OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_848_794/Stats_All_LME';

MaskFile ='/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/RMM_GroupMask_90percent_1849.nii';

for iMeasure=1:length(MeasureSet)
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure}]);
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure}]);
    FileList=[];
    FileListGSR=[];
    for iSub=1:length(SubID)
        FileList{iSub,1}=[DataDir,'/',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
        FileListGSR{iSub,1}=[DataDir,'/',ResultsSet{iMeasure},'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
    end

    
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MDDVsControlT'];
    DependentVolume=FileList;

    [DependentVolume,VoxelSize,theImgFileList, Header] = y_ReadAll(DependentVolume);
    DependentVolume(find(isnan(DependentVolume)))=0;
    [nDim1,nDim2,nDim3,nDim4]=size(DependentVolume);
    [MaskData,MaskVox,MaskHead]=y_ReadRPI(MaskFile);
    MaskData = any(DependentVolume,4) .* MaskData; % skip the voxels with all zeros

    TMatrix=zeros(nDim1,nDim2,nDim3);
    PMatrix=zeros(nDim1,nDim2,nDim3);
    for i=1:nDim1
        fprintf('.');
        for j=1:nDim2
            parfor k=1:nDim3
                addpath(fullfile(matlabroot,'toolbox','stats','stats'))
                if MaskData(i,j,k)
                    y=double(squeeze(DependentVolume(i,j,k,:)));
                    if any(y)
                        lme = fitlmematrix(X,y,Z,G);
                        TMatrix(i,j,k)=lme.Coefficients{2,4}; %
                        PMatrix(i,j,k)=lme.Coefficients{2,6};
                    end
                end
            end
        end
    end
    
    y_Write(TMatrix,Header,[OutputName]);  %y_Write(TF_ForContrast_brain,HeaderTWithDOF,[OutputName,'_',TF_Flag,'_ForContrast','.nii']);
    y_Write(PMatrix,Header,[OutputName,'_p','.nii']); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size)
    
    
    
    
    
    %For GSR
    OutputName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MDDVsControlT'];
    DependentVolume=FileListGSR;

    [DependentVolume,VoxelSize,theImgFileList, Header] = y_ReadAll(DependentVolume);
    DependentVolume(find(isnan(DependentVolume)))=0;
    [nDim1,nDim2,nDim3,nDim4]=size(DependentVolume);
    [MaskData,MaskVox,MaskHead]=y_ReadRPI(MaskFile);
    MaskData = any(DependentVolume,4) .* MaskData; % skip the voxels with all zeros

    TMatrix=zeros(nDim1,nDim2,nDim3);
    PMatrix=zeros(nDim1,nDim2,nDim3);
    for i=1:nDim1
        fprintf('.');
        for j=1:nDim2
            parfor k=1:nDim3
                addpath(fullfile(matlabroot,'toolbox','stats','stats'))
                if MaskData(i,j,k)
                    y=double(squeeze(DependentVolume(i,j,k,:)));
                    if any(y)
                        lme = fitlmematrix(X,y,Z,G);
                        TMatrix(i,j,k)=lme.Coefficients{2,4}; %
                        PMatrix(i,j,k)=lme.Coefficients{2,6};
                    end
                end
            end
        end
    end
    
    y_Write(TMatrix,Header,[OutputName]);  %y_Write(TF_ForContrast_brain,HeaderTWithDOF,[OutputName,'_',TF_Flag,'_ForContrast','.nii']);
    y_Write(PMatrix,Header,[OutputName,'_p','.nii']); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size)
    
end









%DO
OutDir='/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_848_794/Stats_All_LME';
MaskFile ='/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/SubInfo/RMM_GroupMask_90percent_1849.nii';


%GRF Correction

MeasureSet={'ALFF','fALFF','ReHo','DegreeCentrality','VMHC'};
MeasurePrefixSet={'szALFFMap_','szfALFFMap_','szReHoMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','zVMHCMap_'};
ConditionSet={'_FunImgARCW','_FunImgARCW','_FunImgARCWF','_FunImgARCWF','_FunImgARCWFsymS'};
ConditionGSRSet={'_FunImgARglobalCW','_FunImgARglobalCW','_FunImgARglobalCWF','_FunImgARglobalCWF','_FunImgARglobalCWFsymS'};
ResultsSet={'ResultsS','ResultsS','ResultsS','ResultsS','Results'};



IsTwoTailed=1; %!!!
VoxelPThresholdSet_OneTailed=[0.0005];
ClusterPThresholdSet_OneTailed=[0.025];


for iMeasure=1:length(MeasureSet)
    
    for iP=1:length(VoxelPThresholdSet_OneTailed)
        VoxelPThreshold=2*VoxelPThresholdSet_OneTailed(iP); %For y_GRF_Threshold, input two-tailed p values.
        ClusterPThreshold=2*ClusterPThresholdSet_OneTailed(iP); %For y_GRF_Threshold, input two-tailed p values.
        
        CorrectionName=['OneTailedGRF_',num2str(VoxelPThresholdSet_OneTailed(iP)),'_',num2str(ClusterPThresholdSet_OneTailed(iP))];
        
        mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/']);
        mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/']);
        
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MDDVsControlT'];
        OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/MDDVsControlT'];
        y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutName,MaskFile,'T',1636);
        %[Data_Corrected, ClusterSize, Header]=y_GRF_Threshold(StatsImgFile,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutputName,MaskFile,Flag,Df1,Df2,VoxelSize,Header, dLh)

        
        FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MDDVsControlT'];
        OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/MDDVsControlT'];
        y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutName,MaskFile,'T',1636);
    end
end





%FDR Correction
qThreshold=0.05;

for iMeasure=1:length(MeasureSet)
    
    CorrectionName='FDR';
    
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/']);
    mkdir([OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/']);
    
    FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/MDDVsControlT'];
    OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionSet{iMeasure},'/',CorrectionName,'/MDDVsControlT'];
    y_FDR_Image(FileName,qThreshold,OutName,MaskFile,'T',1636);
    
    FileName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/MDDVsControlT'];
    OutName=[OutDir,'/',MeasureSet{iMeasure},ConditionGSRSet{iMeasure},'/',CorrectionName,'/MDDVsControlT'];
    y_FDR_Image(FileName,qThreshold,OutName,MaskFile,'T',1636);
    
end





