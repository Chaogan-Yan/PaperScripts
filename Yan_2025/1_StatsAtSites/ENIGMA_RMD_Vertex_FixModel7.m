function ENIGMA_RMD_Vertex_FixModel7(DataDir,OutDir,Table,MinimumSubjectNumber,IsDPABISurfDataStyle)
% FORMAT ENIGMA_RMD_Vertex_Models(DataDir,OutDir,Table,MinimumSubjectNumber,IsDPABISurfDataStyle)
% Perform models for ENIGMA & REST-meta-MDD collaborative studies on vertex Thickness and Area within each site: 10 models.
% Input:
%   DataDir - Data Dir for stats. Should have AnatSurfLH AnatSurfRH under this dir.
%   OutDir - Output Dir for stats. 
%   Table - Covariates Table.  SHOULD have a column of eTIV!!!
%   MinimumSubjectNumber - Minimum Subject Number in each group. If the number of availabe subjects are less than this number, then the stats will be skipped. default: 10
%   IsDPABISurfDataStyle - If the data is directly come from DPABISurf, then 1. If the ENIGMA script was used, then 0. default: 0
%
% Output:
%   The stats under OutDir.
%___________________________________________________________________________
% Written by YAN Chao-Gan 210519.
% CAS Key Laboratory of Behavioral Science, Institute of Psychology, Beijing, China;
% International Big-Data Research Center for Depression (IBRCD), Institute of Psychology, Chinese Academy of Sciences, Beijing, China;
% Magnetic Resonance Imaging Research Center, Institute of Psychology, Chinese Academy of Sciences, Beijing, China.
% ycg.yan@gmail.com


if ~exist('MinimumSubjectNumber','var') || isempty(MinimumSubjectNumber)
    MinimumSubjectNumber=10;
end

if ~exist('ParallelWorkersNumber','var') || isempty(ParallelWorkersNumber)
    ParallelWorkersNumber=1;
end

if ~exist('IsDPABISurfDataStyle','var') || isempty(IsDPABISurfDataStyle)
    IsDPABISurfDataStyle=0;
end



[DPABIPath, fileN, extn] = fileparts(which('DPABI.m'));


if IsDPABISurfDataStyle
    MeasureStringSuffix='/fsaverage';
else
    MeasureStringSuffix=[];
end


% Do Stats Analysis.






%7. Model 7: Age of Onset ANOVA analysis!

%Age of onset as three groups: 1) adolescent age of onset (< or = 21) and 2) adult age of onset (>21 - i.e., age 22 and above) 3) HC (group coding: 2,1,0). Note: not to be run in the adolescent group only as they will have only adolescent age of onset
Table.AOGroup = Table.AO <= 21;
Table.AOGroup = Table.AOGroup +1;
%Table.AOGroup = Table.AOGroup .* Table.Dx;

Table.AOGroup(find(isnan(Table.AO)))=NaN; %Deal with NaN

Table.AOGroup(find(Table.Dx==0))=0; %Deal with NaN


%First count available subjects
AllCov = [Table.AOGroup, Table.Sex, Table.Age];
%Deal with NaN
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN),:)=[];
NumberOfAdolescentAO = length(find(AllCov(:,1)==2));
NumberOfAdultAO = length(find(AllCov(:,1)==1));
NumberOfHC = length(find(AllCov(:,1)==0));


%7.1. ANOVA for Age of Onset status

% delete([OutDir,'/AnatSurfLH/Thickness',MeasureStringSuffix,'/M7_AOGroup_ANOVA_F.gii']);
% delete([OutDir,'/AnatSurfLH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii']);
% delete([OutDir,'/AnatSurfLH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsHC.gii']);
% delete([OutDir,'/AnatSurfLH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdultAOVsHC.gii']);
% delete([OutDir,'/AnatSurfRH/Thickness',MeasureStringSuffix,'/M7_AOGroup_ANOVA_F.gii']);
% delete([OutDir,'/AnatSurfRH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii']);
% delete([OutDir,'/AnatSurfRH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsHC.gii']);
% delete([OutDir,'/AnatSurfRH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdultAOVsHC.gii']);
% 
% delete([OutDir,'/AnatSurfLH/Area',MeasureStringSuffix,'/M7_AOGroup_ANOVA_F.gii']);
% delete([OutDir,'/AnatSurfLH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii']);
% delete([OutDir,'/AnatSurfLH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsHC.gii']);
% delete([OutDir,'/AnatSurfLH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdultAOVsHC.gii']);
% delete([OutDir,'/AnatSurfRH/Area',MeasureStringSuffix,'/M7_AOGroup_ANOVA_F.gii']);
% delete([OutDir,'/AnatSurfRH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii']);
% delete([OutDir,'/AnatSurfRH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsHC.gii']);
% delete([OutDir,'/AnatSurfRH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdultAOVsHC.gii']);
% 
% 
% delete([OutDir,'/AnatSurfLH/Thickness',MeasureStringSuffix,'/M7_AOGroup_ANOVA_F_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfLH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfLH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsHC_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfLH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdultAOVsHC_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfRH/Thickness',MeasureStringSuffix,'/M7_AOGroup_ANOVA_F_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfRH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfRH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsHC_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfRH/Thickness',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdultAOVsHC_Cohen_f2.gii']);
% 
% delete([OutDir,'/AnatSurfLH/Area',MeasureStringSuffix,'/M7_AOGroup_ANOVA_F_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfLH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfLH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsHC_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfLH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdultAOVsHC_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfRH/Area',MeasureStringSuffix,'/M7_AOGroup_ANOVA_F_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfRH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfRH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdolescentAOVsHC_Cohen_f2.gii']);
% delete([OutDir,'/AnatSurfRH/Area',MeasureStringSuffix,'/M7_AOGroup_PostHoc_AdultAOVsHC_Cohen_f2.gii']);



if (NumberOfAdolescentAO>=MinimumSubjectNumber) && (NumberOfAdultAO>=MinimumSubjectNumber) && (NumberOfHC>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    GroupLabel=Table.AOGroup;
    GroupLabelUnique=[2 1 0];
    Df_Group=length(GroupLabelUnique)-1;
    GroupDummyVariable=zeros(length(Table.SubjID),Df_Group);
    for i=1:Df_Group
        GroupDummyVariable(:,i)=GroupLabel==GroupLabelUnique(i);
    end

    AllCov = [GroupDummyVariable, ones(length(Table.SubjID),1), Table.Sex, Table.Age];
    
    %Deal with NaN
    AllCov(find(HasNaN),:)=[];

    Contrast=zeros(1,size(AllCov,2));
    Contrast(1:Df_Group) = 1;
    
    FileListLH=[];
    FileListRH=[];
    for iSub=1:length(Table.SubjID)
        FileListLH{iSub,1}=sprintf('%s/AnatSurfLH/%s/s%s_space-fsaverage_hemi-L.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
        FileListRH{iSub,1}=sprintf('%s/AnatSurfRH/%s/s%s_space-fsaverage_hemi-R.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
    end
    FileListLH(find(HasNaN),:)=[];
    FileListRH(find(HasNaN),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M7_AOGroup_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M7_AOGroup_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);

    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    GroupLabel=Table.AOGroup;
    GroupLabelUnique=[2 1 0];
    Df_Group=length(GroupLabelUnique)-1;
    GroupDummyVariable=zeros(length(Table.SubjID),Df_Group);
    for i=1:Df_Group
        GroupDummyVariable(:,i)=GroupLabel==GroupLabelUnique(i);
    end

    AllCov = [GroupDummyVariable, ones(length(Table.SubjID),1), Table.Sex, Table.Age, Table.eTIV];
    
    %Deal with NaN
    AllCov(find(HasNaN),:)=[];

    Contrast=zeros(1,size(AllCov,2));
    Contrast(1:Df_Group) = 1;
    
    FileListLH=[];
    FileListRH=[];
    for iSub=1:length(Table.SubjID)
        FileListLH{iSub,1}=sprintf('%s/AnatSurfLH/%s/s%s_space-fsaverage_hemi-L.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
        FileListRH{iSub,1}=sprintf('%s/AnatSurfRH/%s/s%s_space-fsaverage_hemi-R.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
    end
    FileListLH(find(HasNaN),:)=[];
    FileListRH(find(HasNaN),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M7_AOGroup_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M7_AOGroup_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
end


%7.2. Posthoc for Age of Onset status: AdolescentAO vs. AdultAO (2 vs. 1)
if (NumberOfAdolescentAO>=MinimumSubjectNumber) && (NumberOfAdultAO>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AOGroup, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AOGroup==0),:)=[];
    
    %Centering: Let the first column (constant) have the mean effect.
    AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
    
    Contrast=zeros(1,size(AllCov,2));
    Contrast(2)=1;
    
    FileListLH=[];
    FileListRH=[];
    for iSub=1:length(Table.SubjID)
        FileListLH{iSub,1}=sprintf('%s/AnatSurfLH/%s/s%s_space-fsaverage_hemi-L.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
        FileListRH{iSub,1}=sprintf('%s/AnatSurfRH/%s/s%s_space-fsaverage_hemi-R.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
    end
    FileListLH(find(HasNaN | Table.AOGroup==0),:)=[];
    FileListRH(find(HasNaN | Table.AOGroup==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);

    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AOGroup, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AOGroup==0),:)=[];
    
    %Centering: Let the first column (constant) have the mean effect.
    AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
    
    Contrast=zeros(1,size(AllCov,2));
    Contrast(2)=1;
    
    FileListLH=[];
    FileListRH=[];
    for iSub=1:length(Table.SubjID)
        FileListLH{iSub,1}=sprintf('%s/AnatSurfLH/%s/s%s_space-fsaverage_hemi-L.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
        FileListRH{iSub,1}=sprintf('%s/AnatSurfRH/%s/s%s_space-fsaverage_hemi-R.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
    end
    FileListLH(find(HasNaN | Table.AOGroup==0),:)=[];
    FileListRH(find(HasNaN | Table.AOGroup==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end

%7.3. Posthoc for Age of Onset status: AdolescentAO vs. HC (2 vs. 0)
if (NumberOfAdolescentAO>=MinimumSubjectNumber) && (NumberOfHC>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AOGroup, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AOGroup==1),:)=[];
    
    %Centering: Let the first column (constant) have the mean effect.
    AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
    
    Contrast=zeros(1,size(AllCov,2));
    Contrast(2)=1;
    
    FileListLH=[];
    FileListRH=[];
    for iSub=1:length(Table.SubjID)
        FileListLH{iSub,1}=sprintf('%s/AnatSurfLH/%s/s%s_space-fsaverage_hemi-L.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
        FileListRH{iSub,1}=sprintf('%s/AnatSurfRH/%s/s%s_space-fsaverage_hemi-R.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
    end
    FileListLH(find(HasNaN | Table.AOGroup==1),:)=[];
    FileListRH(find(HasNaN | Table.AOGroup==1),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M7_AOGroup_PostHoc_AdolescentAOVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M7_AOGroup_PostHoc_AdolescentAOVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AOGroup, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AOGroup==1),:)=[];
    
    %Centering: Let the first column (constant) have the mean effect.
    AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
    
    Contrast=zeros(1,size(AllCov,2));
    Contrast(2)=1;
    
    FileListLH=[];
    FileListRH=[];
    for iSub=1:length(Table.SubjID)
        FileListLH{iSub,1}=sprintf('%s/AnatSurfLH/%s/s%s_space-fsaverage_hemi-L.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
        FileListRH{iSub,1}=sprintf('%s/AnatSurfRH/%s/s%s_space-fsaverage_hemi-R.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
    end
    FileListLH(find(HasNaN | Table.AOGroup==1),:)=[];
    FileListRH(find(HasNaN | Table.AOGroup==1),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M7_AOGroup_PostHoc_AdolescentAOVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M7_AOGroup_PostHoc_AdolescentAOVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);

end

%7.4. Posthoc for Age of Onset status: AdultAO vs. HC (1 vs. 0)
if (NumberOfAdultAO>=MinimumSubjectNumber) && (NumberOfHC>=MinimumSubjectNumber)
    
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AOGroup, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AOGroup==2),:)=[];
    
    %Centering: Let the first column (constant) have the mean effect.
    AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
    
    Contrast=zeros(1,size(AllCov,2));
    Contrast(2)=1;
    
    FileListLH=[];
    FileListRH=[];
    for iSub=1:length(Table.SubjID)
        FileListLH{iSub,1}=sprintf('%s/AnatSurfLH/%s/s%s_space-fsaverage_hemi-L.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
        FileListRH{iSub,1}=sprintf('%s/AnatSurfRH/%s/s%s_space-fsaverage_hemi-R.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
    end
    FileListLH(find(HasNaN | Table.AOGroup==2),:)=[];
    FileListRH(find(HasNaN | Table.AOGroup==2),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M7_AOGroup_PostHoc_AdultAOVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M7_AOGroup_PostHoc_AdultAOVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AOGroup, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AOGroup==2),:)=[];
    
    %Centering: Let the first column (constant) have the mean effect.
    AllCov(:,2:end) = (AllCov(:,2:end)-repmat(mean(AllCov(:,2:end)),size(AllCov(:,2:end),1),1));
    
    Contrast=zeros(1,size(AllCov,2));
    Contrast(2)=1;
    
    FileListLH=[];
    FileListRH=[];
    for iSub=1:length(Table.SubjID)
        FileListLH{iSub,1}=sprintf('%s/AnatSurfLH/%s/s%s_space-fsaverage_hemi-L.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
        FileListRH{iSub,1}=sprintf('%s/AnatSurfRH/%s/s%s_space-fsaverage_hemi-R.%s.gii',DataDir,MeasureString, Table.SubjID{iSub},MeasureStringLower);
    end
    FileListLH(find(HasNaN | Table.AOGroup==2),:)=[];
    FileListRH(find(HasNaN | Table.AOGroup==2),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M7_AOGroup_PostHoc_AdultAOVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M7_AOGroup_PostHoc_AdultAOVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end




%Save Summary Information for Model 7

AllCov = [Table.AOGroup, Table.Sex, Table.Age];
%Deal with NaN
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN),:)=[];

Race_Ethnicity=Table.Race_Ethnicity;
Race_Ethnicity(find(HasNaN),:)=[];

InfoFile=[OutDir,filesep,'ModelSummaryInfo.txt'];
fid = fopen(InfoFile,'at+');
fprintf(fid,'Summary Information for Model 7: Age of Onset Status\n');
fprintf(fid,'Number of Subjects\t%g\n',size(AllCov,1));
fprintf(fid,'Number of Adolescent Age Onset MDDs\t%g\n',length(find(AllCov(:,1)==2)));
fprintf(fid,'Number of Adult Age Onset MDDs\t%g\n',length(find(AllCov(:,1)==1)));
fprintf(fid,'Number of HCs\t%g\n',length(find(AllCov(:,1)==0)));

fprintf(fid,'Number of Males\t%g\n',length(find(AllCov(:,2)==1)));
fprintf(fid,'Number of Females\t%g\n',length(find(AllCov(:,2)==2)));
fprintf(fid,'Number of Males in Adolescent Age Onset MDD\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Females in Adolescent Age Onset MDD\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Males in Adult Age Onset MDD\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Females in Adult Age Onset MDD\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Males in HC\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Females in HC\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==0))));

fprintf(fid,'Mean Age\t%g\n',mean(AllCov(:,3)));
fprintf(fid,'STD Age\t%g\n',std(AllCov(:,3)));
fprintf(fid,'Mean Age in Adolescent Age Onset MDD\t%g\n',mean(AllCov(find(AllCov(:,1)==2),3)));
fprintf(fid,'STD Age in Adolescent Age Onset MDD\t%g\n',std(AllCov(find(AllCov(:,1)==2),3)));
fprintf(fid,'Mean Age in Adult Age Onset MDD\t%g\n',mean(AllCov(find(AllCov(:,1)==1),3)));
fprintf(fid,'STD Age in Adult Age Onset MDD\t%g\n',std(AllCov(find(AllCov(:,1)==1),3)));
fprintf(fid,'Mean Age in HC\t%g\n',mean(AllCov(find(AllCov(:,1)==0),3)));
fprintf(fid,'STD Age in HC\t%g\n',std(AllCov(find(AllCov(:,1)==0),3)));

fprintf(fid,'Number of White/caucasian\t%g\n',length(find(Race_Ethnicity==1)));
fprintf(fid,'Number of Black/African\t%g\n',length(find(Race_Ethnicity==2)));
fprintf(fid,'Number of Asian\t%g\n',length(find(Race_Ethnicity==3)));
fprintf(fid,'Number of Other Race Ethnicity\t%g\n',length(find(Race_Ethnicity==4)));
fprintf(fid,'Number of NaN Race Ethnicity\t%g\n',length(find(isnan(Race_Ethnicity))));
fprintf(fid,'Number of White/caucasian in Adolescent Age Onset MDD\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Black/African in Adolescent Age Onset MDD\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Asian in Adolescent Age Onset MDD\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Other Race Ethnicity in Adolescent Age Onset MDD\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==2))));
fprintf(fid,'Number of NaN Race Ethnicity in Adolescent Age Onset MDD\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==2))));
fprintf(fid,'Number of White/caucasian in Adult Age Onset MDD\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Black/African in Adult Age Onset MDD\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Asian in Adult Age Onset MDD\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Other Race Ethnicity in Adult Age Onset MDD\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==1))));
fprintf(fid,'Number of NaN Race Ethnicity in Adult Age Onset MDD\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==1))));
fprintf(fid,'Number of White/caucasian in HC\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Black/African in HC\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Asian in HC\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Other Race Ethnicity in HC\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==0))));
fprintf(fid,'Number of NaN Race Ethnicity in HC\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==0))));

fprintf(fid,'\n\n\n');
fclose(fid);





fprintf('\n\tPerform Stats for ENIGMA & REST-meta-MDD collaborative studies on vertex Thickness and Area within each site: fix Model 7: Finished!!!\n');


