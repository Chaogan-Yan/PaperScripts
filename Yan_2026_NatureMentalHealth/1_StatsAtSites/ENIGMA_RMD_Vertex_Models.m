function ENIGMA_RMD_Vertex_Models(DataDir,OutDir,Table,MinimumSubjectNumber,IsDPABISurfDataStyle)
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

%1. Model 1: Dx


NumberOfMDD = length(find(Table.Dx==1));
NumberOfHC = length(find(Table.Dx==0));

if (NumberOfMDD>=MinimumSubjectNumber) && (NumberOfHC>=MinimumSubjectNumber)
    
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Dx, Table.Sex, Table.Age];
    
    %Deal with NaN
    HasNaN = isnan(sum(AllCov,2));
    AllCov(find(HasNaN),:)=[];
    
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
    FileListLH(find(HasNaN),:)=[];
    FileListRH(find(HasNaN),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M1_Dx.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M1_Dx.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Dx, Table.Sex, Table.Age, Table.eTIV];
    
    %Deal with NaN
    HasNaN = isnan(sum(AllCov,2));
    AllCov(find(HasNaN),:)=[];
    
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
    FileListLH(find(HasNaN),:)=[];
    FileListRH(find(HasNaN),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M1_Dx.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M1_Dx.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    
    %2. Model 2: Dx by Age
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), (Table.Dx-nanmean(Table.Dx)).*(Table.Age-nanmean(Table.Age)), Table.Dx, Table.Sex, Table.Age];
    
    %Deal with NaN
    HasNaN = isnan(sum(AllCov,2));
    AllCov(find(HasNaN),:)=[];
    
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
    FileListLH(find(HasNaN),:)=[];
    FileListRH(find(HasNaN),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M2_DxByAge.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M2_DxByAge.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), (Table.Dx-nanmean(Table.Dx)).*(Table.Age-nanmean(Table.Age)), Table.Dx, Table.Sex, Table.Age, Table.eTIV];
    
    %Deal with NaN
    HasNaN = isnan(sum(AllCov,2));
    AllCov(find(HasNaN),:)=[];
    
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
    FileListLH(find(HasNaN),:)=[];
    FileListRH(find(HasNaN),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M2_DxByAge.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M2_DxByAge.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    
    
    %3. Model 3: Dx by Sex
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), (Table.Dx-nanmean(Table.Dx)).*(Table.Sex-nanmean(Table.Sex)), Table.Dx, Table.Sex, Table.Age];
    
    %Deal with NaN
    HasNaN = isnan(sum(AllCov,2));
    AllCov(find(HasNaN),:)=[];
    
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
    FileListLH(find(HasNaN),:)=[];
    FileListRH(find(HasNaN),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M3_DxBySex.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M3_DxBySex.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), (Table.Dx-nanmean(Table.Dx)).*(Table.Sex-nanmean(Table.Sex)), Table.Dx, Table.Sex, Table.Age, Table.eTIV];
    
    %Deal with NaN
    HasNaN = isnan(sum(AllCov,2));
    AllCov(find(HasNaN),:)=[];
    
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
    FileListLH(find(HasNaN),:)=[];
    FileListRH(find(HasNaN),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M3_DxBySex.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M3_DxBySex.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
end


%Save Summary Information for Models 1~3

AllCov = [Table.Dx, Table.Sex, Table.Age];
%Deal with NaN
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN),:)=[];

Race_Ethnicity=Table.Race_Ethnicity;
Race_Ethnicity(find(HasNaN),:)=[];

InfoFile=[OutDir,filesep,'ModelSummaryInfo.txt'];
fid = fopen(InfoFile,'at+');
fprintf(fid,'Summary Information for Models 1~3: Dx; Dx by Age; Dx by Sex\n');
fprintf(fid,'Number of Subjects\t%g\n',size(AllCov,1));
fprintf(fid,'Number of MDDs\t%g\n',length(find(AllCov(:,1)==1)));
fprintf(fid,'Number of HCs\t%g\n',length(find(AllCov(:,1)==0)));

fprintf(fid,'Number of Males\t%g\n',length(find(AllCov(:,2)==1)));
fprintf(fid,'Number of Females\t%g\n',length(find(AllCov(:,2)==2)));
fprintf(fid,'Number of Males in MDD\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Females in MDD\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Males in HC\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Females in HC\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==0))));

fprintf(fid,'Mean Age\t%g\n',mean(AllCov(:,3)));
fprintf(fid,'STD Age\t%g\n',std(AllCov(:,3)));
fprintf(fid,'Mean Age in MDD\t%g\n',mean(AllCov(find(AllCov(:,1)==1),3)));
fprintf(fid,'STD Age in MDD\t%g\n',std(AllCov(find(AllCov(:,1)==1),3)));
fprintf(fid,'Mean Age in HC\t%g\n',mean(AllCov(find(AllCov(:,1)==0),3)));
fprintf(fid,'STD Age in HC\t%g\n',std(AllCov(find(AllCov(:,1)==0),3)));

fprintf(fid,'Number of White/caucasian\t%g\n',length(find(Race_Ethnicity==1)));
fprintf(fid,'Number of Black/African\t%g\n',length(find(Race_Ethnicity==2)));
fprintf(fid,'Number of Asian\t%g\n',length(find(Race_Ethnicity==3)));
fprintf(fid,'Number of Other Race Ethnicity\t%g\n',length(find(Race_Ethnicity==4)));
fprintf(fid,'Number of NaN Race Ethnicity\t%g\n',length(find(isnan(Race_Ethnicity))));
fprintf(fid,'Number of White/caucasian in MDD\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Black/African in MDD\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Asian in MDD\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Other Race Ethnicity in MDD\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==1))));
fprintf(fid,'Number of NaN Race Ethnicity in MDD\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==1))));
fprintf(fid,'Number of White/caucasian in HC\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Black/African in HC\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Asian in HC\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Other Race Ethnicity in HC\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==0))));
fprintf(fid,'Number of NaN Race Ethnicity in HC\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==0))));

fprintf(fid,'\n\n\n');
fclose(fid);


%4. Model 4: Recurrence Status

%First count available subjects
AllCov = [Table.Recur, Table.Sex, Table.Age];
%Deal with NaN
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN),:)=[];
NumberOfRecurrentMDD = length(find(AllCov(:,1)==2));
NumberOfFirstEpisodeMDD = length(find(AllCov(:,1)==1));
NumberOfHC = length(find(AllCov(:,1)==0));

%4.1. ANOVA for Recurrence status

if (NumberOfRecurrentMDD>=MinimumSubjectNumber) && (NumberOfFirstEpisodeMDD>=MinimumSubjectNumber) && (NumberOfHC>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    GroupLabel=Table.Recur;
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
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M4_Recur_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M4_Recur_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);

    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    GroupLabel=Table.Recur;
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
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M4_Recur_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M4_Recur_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
end


%4.2. Posthoc for Recurrence status: RecurrentMDD vs. FirstEpisodeMDD (2 vs. 1)
if (NumberOfRecurrentMDD>=MinimumSubjectNumber) && (NumberOfFirstEpisodeMDD>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Recur, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Recur==0),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Recur==0),:)=[];
    FileListRH(find(HasNaN | Table.Recur==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M4_Recur_PostHoc_RecurrentVsFirstEpisode.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M4_Recur_PostHoc_RecurrentVsFirstEpisode.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Recur, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Recur==0),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Recur==0),:)=[];
    FileListRH(find(HasNaN | Table.Recur==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M4_Recur_PostHoc_RecurrentVsFirstEpisode.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M4_Recur_PostHoc_RecurrentVsFirstEpisode.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end

%4.3. Posthoc for Recurrence status: RecurrentMDD vs. HC (2 vs. 0)
if (NumberOfRecurrentMDD>=MinimumSubjectNumber) && (NumberOfHC>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Recur, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Recur==1),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Recur==1),:)=[];
    FileListRH(find(HasNaN | Table.Recur==1),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M4_Recur_PostHoc_RecurrentVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M4_Recur_PostHoc_RecurrentVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Recur, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Recur==1),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Recur==1),:)=[];
    FileListRH(find(HasNaN | Table.Recur==1),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M4_Recur_PostHoc_RecurrentVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M4_Recur_PostHoc_RecurrentVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);

end

%4.4. Posthoc for Recurrence status: FirstEpisodeMDD vs. HC (1 vs. 0)
if (NumberOfRecurrentMDD>=MinimumSubjectNumber) && (NumberOfFirstEpisodeMDD>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Recur, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Recur==2),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Recur==2),:)=[];
    FileListRH(find(HasNaN | Table.Recur==2),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M4_Recur_PostHoc_FirstEpisodeVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M4_Recur_PostHoc_FirstEpisodeVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Recur, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Recur==2),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Recur==2),:)=[];
    FileListRH(find(HasNaN | Table.Recur==2),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M4_Recur_PostHoc_FirstEpisodeVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M4_Recur_PostHoc_FirstEpisodeVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end


%Save Summary Information for Model 4

AllCov = [Table.Recur, Table.Sex, Table.Age];
%Deal with NaN
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN),:)=[];

Race_Ethnicity=Table.Race_Ethnicity;
Race_Ethnicity(find(HasNaN),:)=[];

InfoFile=[OutDir,filesep,'ModelSummaryInfo.txt'];
fid = fopen(InfoFile,'at+');
fprintf(fid,'Summary Information for Model 4: Recurrence Status\n');
fprintf(fid,'Number of Subjects\t%g\n',size(AllCov,1));
fprintf(fid,'Number of Recurrent MDDs\t%g\n',length(find(AllCov(:,1)==2)));
fprintf(fid,'Number of First Episode MDDs\t%g\n',length(find(AllCov(:,1)==1)));
fprintf(fid,'Number of HCs\t%g\n',length(find(AllCov(:,1)==0)));

fprintf(fid,'Number of Males\t%g\n',length(find(AllCov(:,2)==1)));
fprintf(fid,'Number of Females\t%g\n',length(find(AllCov(:,2)==2)));
fprintf(fid,'Number of Males in Recurrent MDD\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Females in Recurrent MDD\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Males in First Episode MDD\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Females in First Episode MDD\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Males in HC\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Females in HC\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==0))));

fprintf(fid,'Mean Age\t%g\n',mean(AllCov(:,3)));
fprintf(fid,'STD Age\t%g\n',std(AllCov(:,3)));
fprintf(fid,'Mean Age in Recurrent MDD\t%g\n',mean(AllCov(find(AllCov(:,1)==2),3)));
fprintf(fid,'STD Age in Recurrent MDD\t%g\n',std(AllCov(find(AllCov(:,1)==2),3)));
fprintf(fid,'Mean Age in First Episode MDD\t%g\n',mean(AllCov(find(AllCov(:,1)==1),3)));
fprintf(fid,'STD Age in First Episode MDD\t%g\n',std(AllCov(find(AllCov(:,1)==1),3)));
fprintf(fid,'Mean Age in HC\t%g\n',mean(AllCov(find(AllCov(:,1)==0),3)));
fprintf(fid,'STD Age in HC\t%g\n',std(AllCov(find(AllCov(:,1)==0),3)));

fprintf(fid,'Number of White/caucasian\t%g\n',length(find(Race_Ethnicity==1)));
fprintf(fid,'Number of Black/African\t%g\n',length(find(Race_Ethnicity==2)));
fprintf(fid,'Number of Asian\t%g\n',length(find(Race_Ethnicity==3)));
fprintf(fid,'Number of Other Race Ethnicity\t%g\n',length(find(Race_Ethnicity==4)));
fprintf(fid,'Number of NaN Race Ethnicity\t%g\n',length(find(isnan(Race_Ethnicity))));
fprintf(fid,'Number of White/caucasian in Recurrent MDD\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Black/African in Recurrent MDD\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Asian in Recurrent MDD\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Other Race Ethnicity in Recurrent MDD\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==2))));
fprintf(fid,'Number of NaN Race Ethnicity in Recurrent MDD\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==2))));
fprintf(fid,'Number of White/caucasian in First Episode MDD\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Black/African in First Episode MDD\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Asian in First Episode MDD\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Other Race Ethnicity in First Episode MDD\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==1))));
fprintf(fid,'Number of NaN Race Ethnicity in First Episode MDD\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==1))));
fprintf(fid,'Number of White/caucasian in HC\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Black/African in HC\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Asian in HC\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Other Race Ethnicity in HC\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==0))));
fprintf(fid,'Number of NaN Race Ethnicity in HC\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==0))));

fprintf(fid,'\n\n\n');
fclose(fid);




%5. Model 5: Remission status

%First count available subjects
AllCov = [Table.Rem, Table.Sex, Table.Age];
%Deal with NaN
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN),:)=[];
NumberOfAcutelyDepressed = length(find(AllCov(:,1)==2));
NumberOfRemittedeMDD = length(find(AllCov(:,1)==1));
NumberOfHC = length(find(AllCov(:,1)==0));


%5.1. ANOVA for Remission status

if (NumberOfAcutelyDepressed>=MinimumSubjectNumber) && (NumberOfRemittedeMDD>=MinimumSubjectNumber) && (NumberOfHC>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    GroupLabel=Table.Rem;
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
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M5_Rem_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M5_Rem_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);

    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    GroupLabel=Table.Rem;
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
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M5_Rem_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M5_Rem_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
end


%5.2. Posthoc for Remission status: AcutelyDepressed vs. RemittedeMDD (2 vs. 1)
if (NumberOfAcutelyDepressed>=MinimumSubjectNumber) && (NumberOfRemittedeMDD>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Rem, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Rem==0),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Rem==0),:)=[];
    FileListRH(find(HasNaN | Table.Rem==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M5_Rem_PostHoc_AcutelyDepressedVsRemittedeMDD.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M5_Rem_PostHoc_AcutelyDepressedVsRemittedeMDD.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);

    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Rem, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Rem==0),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Rem==0),:)=[];
    FileListRH(find(HasNaN | Table.Rem==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M5_Rem_PostHoc_AcutelyDepressedVsRemittedeMDD.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M5_Rem_PostHoc_AcutelyDepressedVsRemittedeMDD.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end

%5.3. Posthoc for Remission status: AcutelyDepressed vs. HC (2 vs. 0)
if (NumberOfAcutelyDepressed>=MinimumSubjectNumber) && (NumberOfHC>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Rem, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Rem==1),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Rem==1),:)=[];
    FileListRH(find(HasNaN | Table.Rem==1),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M5_Rem_PostHoc_AcutelyDepressedVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M5_Rem_PostHoc_AcutelyDepressedVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Rem, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Rem==1),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Rem==1),:)=[];
    FileListRH(find(HasNaN | Table.Rem==1),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M5_Rem_PostHoc_AcutelyDepressedVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M5_Rem_PostHoc_AcutelyDepressedVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);

end

%5.4. Posthoc for Remission status: RemittedeMDD vs. HC (1 vs. 0)
if (NumberOfAcutelyDepressed>=MinimumSubjectNumber) && (NumberOfRemittedeMDD>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Rem, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Rem==2),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Rem==2),:)=[];
    FileListRH(find(HasNaN | Table.Rem==2),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M5_Rem_PostHoc_RemittedeMDDVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M5_Rem_PostHoc_RemittedeMDDVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.Rem, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.Rem==2),:)=[];
    
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
    FileListLH(find(HasNaN | Table.Rem==2),:)=[];
    FileListRH(find(HasNaN | Table.Rem==2),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M5_Rem_PostHoc_RemittedeMDDVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M5_Rem_PostHoc_RemittedeMDDVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end




%Save Summary Information for Model 5

AllCov = [Table.Rem, Table.Sex, Table.Age];
%Deal with NaN
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN),:)=[];

Race_Ethnicity=Table.Race_Ethnicity;
Race_Ethnicity(find(HasNaN),:)=[];

InfoFile=[OutDir,filesep,'ModelSummaryInfo.txt'];
fid = fopen(InfoFile,'at+');
fprintf(fid,'Summary Information for Model 5: Remission Status\n');
fprintf(fid,'Number of Subjects\t%g\n',size(AllCov,1));
fprintf(fid,'Number of Acutely Depressed MDDs\t%g\n',length(find(AllCov(:,1)==2)));
fprintf(fid,'Number of Remitted MDDs\t%g\n',length(find(AllCov(:,1)==1)));
fprintf(fid,'Number of HCs\t%g\n',length(find(AllCov(:,1)==0)));

fprintf(fid,'Number of Males\t%g\n',length(find(AllCov(:,2)==1)));
fprintf(fid,'Number of Females\t%g\n',length(find(AllCov(:,2)==2)));
fprintf(fid,'Number of Males in Acutely Depressed MDD\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Females in Acutely Depressed MDD\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Males in Remitted MDD\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Females in Remitted MDD\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Males in HC\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Females in HC\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==0))));

fprintf(fid,'Mean Age\t%g\n',mean(AllCov(:,3)));
fprintf(fid,'STD Age\t%g\n',std(AllCov(:,3)));
fprintf(fid,'Mean Age in Acutely Depressed MDD\t%g\n',mean(AllCov(find(AllCov(:,1)==2),3)));
fprintf(fid,'STD Age in Acutely Depressed MDD\t%g\n',std(AllCov(find(AllCov(:,1)==2),3)));
fprintf(fid,'Mean Age in Remitted MDD\t%g\n',mean(AllCov(find(AllCov(:,1)==1),3)));
fprintf(fid,'STD Age in Remitted MDD\t%g\n',std(AllCov(find(AllCov(:,1)==1),3)));
fprintf(fid,'Mean Age in HC\t%g\n',mean(AllCov(find(AllCov(:,1)==0),3)));
fprintf(fid,'STD Age in HC\t%g\n',std(AllCov(find(AllCov(:,1)==0),3)));

fprintf(fid,'Number of White/caucasian\t%g\n',length(find(Race_Ethnicity==1)));
fprintf(fid,'Number of Black/African\t%g\n',length(find(Race_Ethnicity==2)));
fprintf(fid,'Number of Asian\t%g\n',length(find(Race_Ethnicity==3)));
fprintf(fid,'Number of Other Race Ethnicity\t%g\n',length(find(Race_Ethnicity==4)));
fprintf(fid,'Number of NaN Race Ethnicity\t%g\n',length(find(isnan(Race_Ethnicity))));
fprintf(fid,'Number of White/caucasian in Acutely Depressed MDD\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Black/African in Acutely Depressed MDD\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Asian in Acutely Depressed MDD\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Other Race Ethnicity in Acutely Depressed MDD\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==2))));
fprintf(fid,'Number of NaN Race Ethnicity in Acutely Depressed MDD\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==2))));
fprintf(fid,'Number of White/caucasian in Remitted MDD\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Black/African in Remitted MDD\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Asian in Remitted MDD\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Other Race Ethnicity in Remitted MDD\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==1))));
fprintf(fid,'Number of NaN Race Ethnicity in Remitted MDD\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==1))));
fprintf(fid,'Number of White/caucasian in HC\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Black/African in HC\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Asian in HC\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Other Race Ethnicity in HC\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==0))));
fprintf(fid,'Number of NaN Race Ethnicity in HC\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==0))));

fprintf(fid,'\n\n\n');
fclose(fid);





%6. Model 6: Age of Onset

Measure='Thickness';
MeasureString = [Measure,MeasureStringSuffix];
MeasureStringLower = lower(Measure);

AllCov = [ones(length(Table.SubjID),1), Table.AO, Table.Sex, Table.Age];

%Deal with NaN and remove HCs
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN | Table.Dx==0),:)=[];

if size(AllCov,1) >= MinimumSubjectNumber
    
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
    FileListLH(find(HasNaN | Table.Dx==0),:)=[];
    FileListRH(find(HasNaN | Table.Dx==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M6_AgeOfOnset.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M6_AgeOfOnset.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end

Measure='Area';
MeasureString = [Measure,MeasureStringSuffix];
MeasureStringLower = lower(Measure);

AllCov = [ones(length(Table.SubjID),1), Table.AO, Table.Sex, Table.Age, Table.eTIV];

%Deal with NaN and remove HCs
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN | Table.Dx==0),:)=[];

if size(AllCov,1) >= MinimumSubjectNumber
    
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
    FileListLH(find(HasNaN | Table.Dx==0),:)=[];
    FileListRH(find(HasNaN | Table.Dx==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M6_AgeOfOnsetgii.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M6_AgeOfOnset.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end



%Save Summary Information for Model 6
AllCov = [Table.AO, Table.Sex, Table.Age];

%Deal with NaN and remove HCs
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN | Table.Dx==0),:)=[];

Race_Ethnicity=Table.Race_Ethnicity;
Race_Ethnicity(find(HasNaN | Table.Dx==0),:)=[];

InfoFile=[OutDir,filesep,'ModelSummaryInfo.txt'];
fid = fopen(InfoFile,'at+');
fprintf(fid,'Summary Information for Model 6: Age of Onset\n');
fprintf(fid,'Number of Subjects\t%g\n',size(AllCov,1));

fprintf(fid,'Mean Age of Onset\t%g\n',mean(AllCov(:,1)));
fprintf(fid,'STD Age of Onset\t%g\n',std(AllCov(:,1)));

fprintf(fid,'Number of Males\t%g\n',length(find(AllCov(:,2)==1)));
fprintf(fid,'Number of Females\t%g\n',length(find(AllCov(:,2)==2)));

fprintf(fid,'Mean Age\t%g\n',mean(AllCov(:,3)));
fprintf(fid,'STD Age\t%g\n',std(AllCov(:,3)));

fprintf(fid,'Number of White/caucasian\t%g\n',length(find(Race_Ethnicity==1)));
fprintf(fid,'Number of Black/African\t%g\n',length(find(Race_Ethnicity==2)));
fprintf(fid,'Number of Asian\t%g\n',length(find(Race_Ethnicity==3)));
fprintf(fid,'Number of Other Race Ethnicity\t%g\n',length(find(Race_Ethnicity==4)));
fprintf(fid,'Number of NaN Race Ethnicity\t%g\n',length(find(isnan(Race_Ethnicity))));

fprintf(fid,'\n\n\n');
fclose(fid);









%7. Model 7: Age of Onset ANOVA analysis!

%Age of onset as three groups: 1) adolescent age of onset (< or = 21) and 2) adult age of onset (>21 - i.e., age 22 and above) 3) HC (group coding: 2,1,0). Note: not to be run in the adolescent group only as they will have only adolescent age of onset
Table.AOGroup = Table.AO <= 21;
Table.AOGroup = Table.AOGroup +1;
Table.AOGroup = Table.AOGroup .* Table.Dx;


Table.AOGroup(find(isnan(Table.AO)))=NaN; %Deal with NaN %YAN Chao-Gan, Added on 20230803. Fixed for Model 7!!!


%First count available subjects
AllCov = [Table.AOGroup, Table.Sex, Table.Age];
%Deal with NaN
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN),:)=[];
NumberOfAdolescentAO = length(find(AllCov(:,1)==2));
NumberOfAdultAO = length(find(AllCov(:,1)==1));
NumberOfHC = length(find(AllCov(:,1)==0));


%7.1. ANOVA for Age of Onset status

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
if (NumberOfAdolescentAO>=MinimumSubjectNumber) && (NumberOfAdultAO>=MinimumSubjectNumber)
    
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


















%8. Model 8: BDI severity of depressive symptoms

Measure='Thickness';
MeasureString = [Measure,MeasureStringSuffix];
MeasureStringLower = lower(Measure);

AllCov = [ones(length(Table.SubjID),1), Table.BDI, Table.Sex, Table.Age];

%Deal with NaN and remove HCs
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN | Table.Dx==0),:)=[];

if size(AllCov,1) >= MinimumSubjectNumber
    
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
    FileListLH(find(HasNaN | Table.Dx==0),:)=[];
    FileListRH(find(HasNaN | Table.Dx==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M8_BDI.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M8_BDI.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end

Measure='Area';
MeasureString = [Measure,MeasureStringSuffix];
MeasureStringLower = lower(Measure);

AllCov = [ones(length(Table.SubjID),1), Table.BDI, Table.Sex, Table.Age, Table.eTIV];

%Deal with NaN and remove HCs
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN | Table.Dx==0),:)=[];

if size(AllCov,1) >= MinimumSubjectNumber
    
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
    FileListLH(find(HasNaN | Table.Dx==0),:)=[];
    FileListRH(find(HasNaN | Table.Dx==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M8_BDI.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M8_BDI.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end



%Save Summary Information for Model 8
AllCov = [Table.BDI, Table.Sex, Table.Age];

%Deal with NaN and remove HCs
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN | Table.Dx==0),:)=[];

Race_Ethnicity=Table.Race_Ethnicity;
Race_Ethnicity(find(HasNaN | Table.Dx==0),:)=[];

InfoFile=[OutDir,filesep,'ModelSummaryInfo.txt'];
fid = fopen(InfoFile,'at+');
fprintf(fid,'Summary Information for Model 8: BDI severity of depressive symptoms\n');
fprintf(fid,'Number of Subjects\t%g\n',size(AllCov,1));

fprintf(fid,'Mean BDI of Onset\t%g\n',mean(AllCov(:,1)));
fprintf(fid,'STD BDI of Onset\t%g\n',std(AllCov(:,1)));

fprintf(fid,'Number of Males\t%g\n',length(find(AllCov(:,2)==1)));
fprintf(fid,'Number of Females\t%g\n',length(find(AllCov(:,2)==2)));

fprintf(fid,'Mean Age\t%g\n',mean(AllCov(:,3)));
fprintf(fid,'STD Age\t%g\n',std(AllCov(:,3)));

fprintf(fid,'Number of White/caucasian\t%g\n',length(find(Race_Ethnicity==1)));
fprintf(fid,'Number of Black/African\t%g\n',length(find(Race_Ethnicity==2)));
fprintf(fid,'Number of Asian\t%g\n',length(find(Race_Ethnicity==3)));
fprintf(fid,'Number of Other Race Ethnicity\t%g\n',length(find(Race_Ethnicity==4)));
fprintf(fid,'Number of NaN Race Ethnicity\t%g\n',length(find(isnan(Race_Ethnicity))));

fprintf(fid,'\n\n\n');
fclose(fid);





%9. Model 9: HDRS severity of depressive symptoms

Measure='Thickness';
MeasureString = [Measure,MeasureStringSuffix];
MeasureStringLower = lower(Measure);

AllCov = [ones(length(Table.SubjID),1), Table.HDRS, Table.Sex, Table.Age];

%Deal with NaN and remove HCs
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN | Table.Dx==0),:)=[];

if size(AllCov,1) >= MinimumSubjectNumber
    
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
    FileListLH(find(HasNaN | Table.Dx==0),:)=[];
    FileListRH(find(HasNaN | Table.Dx==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M9_HDRS.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M9_HDRS.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end

Measure='Area';
MeasureString = [Measure,MeasureStringSuffix];
MeasureStringLower = lower(Measure);

AllCov = [ones(length(Table.SubjID),1), Table.HDRS, Table.Sex, Table.Age, Table.eTIV];

%Deal with NaN and remove HCs
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN | Table.Dx==0),:)=[];

if size(AllCov,1) >= MinimumSubjectNumber
    
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
    FileListLH(find(HasNaN | Table.Dx==0),:)=[];
    FileListRH(find(HasNaN | Table.Dx==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M9_HDRS.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M9_HDRS.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end


%Save Summary Information for Model 9
AllCov = [Table.HDRS, Table.Sex, Table.Age];

%Deal with NaN and remove HCs
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN | Table.Dx==0),:)=[];

Race_Ethnicity=Table.Race_Ethnicity;
Race_Ethnicity(find(HasNaN | Table.Dx==0),:)=[];

InfoFile=[OutDir,filesep,'ModelSummaryInfo.txt'];
fid = fopen(InfoFile,'at+');
fprintf(fid,'Summary Information for Model 9: HDRS severity of depressive symptoms\n');
fprintf(fid,'Number of Subjects\t%g\n',size(AllCov,1));

fprintf(fid,'Mean HDRS of Onset\t%g\n',mean(AllCov(:,1)));
fprintf(fid,'STD HDRS of Onset\t%g\n',std(AllCov(:,1)));

fprintf(fid,'Number of Males\t%g\n',length(find(AllCov(:,2)==1)));
fprintf(fid,'Number of Females\t%g\n',length(find(AllCov(:,2)==2)));

fprintf(fid,'Mean Age\t%g\n',mean(AllCov(:,3)));
fprintf(fid,'STD Age\t%g\n',std(AllCov(:,3)));

fprintf(fid,'Number of White/caucasian\t%g\n',length(find(Race_Ethnicity==1)));
fprintf(fid,'Number of Black/African\t%g\n',length(find(Race_Ethnicity==2)));
fprintf(fid,'Number of Asian\t%g\n',length(find(Race_Ethnicity==3)));
fprintf(fid,'Number of Other Race Ethnicity\t%g\n',length(find(Race_Ethnicity==4)));
fprintf(fid,'Number of NaN Race Ethnicity\t%g\n',length(find(isnan(Race_Ethnicity))));

fprintf(fid,'\n\n\n');
fclose(fid);






%10. Model 10: Antidepressant medication use at time of scan

%First count available subjects
AllCov = [Table.AD, Table.Sex, Table.Age];
%Deal with NaN
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN),:)=[];
NumberOfAntidepressantUse = length(find(AllCov(:,1)==2));
NumberOfAntidepressantFree = length(find(AllCov(:,1)==1));
NumberOfHC = length(find(AllCov(:,1)==0));


%10.1. ANOVA for Antidepressant status

if (NumberOfAntidepressantUse>=MinimumSubjectNumber) && (NumberOfAntidepressantFree>=MinimumSubjectNumber) && (NumberOfHC>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    GroupLabel=Table.AD;
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
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M10_Antidepressant_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M10_Antidepressant_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);

    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    GroupLabel=Table.AD;
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
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M10_Antidepressant_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M10_Antidepressant_ANOVA_F.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'F',0);
end


%10.2. Posthoc for Antidepressant status: AntidepressantUse vs. AntidepressantFree (2 vs. 1)
if (NumberOfAntidepressantUse>=MinimumSubjectNumber) && (NumberOfAntidepressantFree>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AD, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AD==0),:)=[];
    
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
    FileListLH(find(HasNaN | Table.AD==0),:)=[];
    FileListRH(find(HasNaN | Table.AD==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantUseVsAntidepressantFree.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantUseVsAntidepressantFree.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);

    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AD, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AD==0),:)=[];
    
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
    FileListLH(find(HasNaN | Table.AD==0),:)=[];
    FileListRH(find(HasNaN | Table.AD==0),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantUseVsAntidepressantFree.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantUseVsAntidepressantFree.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end

%10.3. Posthoc for Antidepressant status: AntidepressantUse vs. HC (2 vs. 0)
if (NumberOfAntidepressantUse>=MinimumSubjectNumber) && (NumberOfHC>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AD, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AD==1),:)=[];
    
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
    FileListLH(find(HasNaN | Table.AD==1),:)=[];
    FileListRH(find(HasNaN | Table.AD==1),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantUseVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantUseVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AD, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AD==1),:)=[];
    
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
    FileListLH(find(HasNaN | Table.AD==1),:)=[];
    FileListRH(find(HasNaN | Table.AD==1),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantUseVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantUseVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);

end

%10.4. Posthoc for Antidepressant status: AntidepressantFree vs. HC (1 vs. 0)
if (NumberOfAntidepressantUse>=MinimumSubjectNumber) && (NumberOfAntidepressantFree>=MinimumSubjectNumber)
    
    Measure='Thickness';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AD, Table.Sex, Table.Age];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AD==2),:)=[];
    
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
    FileListLH(find(HasNaN | Table.AD==2),:)=[];
    FileListRH(find(HasNaN | Table.AD==2),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantFreeVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantFreeVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    
    Measure='Area';
    MeasureString = [Measure,MeasureStringSuffix];
    MeasureStringLower = lower(Measure);
    
    AllCov = [ones(length(Table.SubjID),1), Table.AD, Table.Sex, Table.Age, Table.eTIV];

    %Deal with NaN and Remove HC
    AllCov(find(HasNaN | Table.AD==2),:)=[];
    
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
    FileListLH(find(HasNaN | Table.AD==2),:)=[];
    FileListRH(find(HasNaN | Table.AD==2),:)=[];
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfLH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfLH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantFreeVsHC.gii'];
    y_GroupAnalysis_Image(FileListLH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
    
    MaskFile=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    mkdir([OutDir,'/AnatSurfRH/',MeasureString])
    OutputName=[OutDir,'/AnatSurfRH/',MeasureString,'/M10_Antidepressant_PostHoc_AntidepressantFreeVsHC.gii'];
    y_GroupAnalysis_Image(FileListRH,AllCov,OutputName,MaskFile,[],Contrast,'T',0);
end



%Save Summary Information for Model 10

AllCov = [Table.AD, Table.Sex, Table.Age];
%Deal with NaN
HasNaN = isnan(sum(AllCov,2));
AllCov(find(HasNaN),:)=[];

Race_Ethnicity=Table.Race_Ethnicity;
Race_Ethnicity(find(HasNaN),:)=[];

InfoFile=[OutDir,filesep,'ModelSummaryInfo.txt'];
fid = fopen(InfoFile,'at+');
fprintf(fid,'Summary Information for Model 10: Antidepressant medication use at time of scan\n');
fprintf(fid,'Number of Subjects\t%g\n',size(AllCov,1));
fprintf(fid,'Number of Antidepressant Use MDDs\t%g\n',length(find(AllCov(:,1)==2)));
fprintf(fid,'Number of Antidepressant Free MDDs\t%g\n',length(find(AllCov(:,1)==1)));
fprintf(fid,'Number of HCs\t%g\n',length(find(AllCov(:,1)==0)));

fprintf(fid,'Number of Males\t%g\n',length(find(AllCov(:,2)==1)));
fprintf(fid,'Number of Females\t%g\n',length(find(AllCov(:,2)==2)));
fprintf(fid,'Number of Males in Antidepressant Use MDD\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Females in Antidepressant Use MDD\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Males in Antidepressant Free\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Females in Antidepressant Free\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Males in HC\t%g\n',length(find((AllCov(:,2)==1) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Females in HC\t%g\n',length(find((AllCov(:,2)==2) & (AllCov(:,1)==0))));

fprintf(fid,'Mean Age\t%g\n',mean(AllCov(:,3)));
fprintf(fid,'STD Age\t%g\n',std(AllCov(:,3)));
fprintf(fid,'Mean Age in Antidepressant Use MDD\t%g\n',mean(AllCov(find(AllCov(:,1)==2),3)));
fprintf(fid,'STD Age in Antidepressant Use MDD\t%g\n',std(AllCov(find(AllCov(:,1)==2),3)));
fprintf(fid,'Mean Age in Antidepressant Free\t%g\n',mean(AllCov(find(AllCov(:,1)==1),3)));
fprintf(fid,'STD Age in Antidepressant Free\t%g\n',std(AllCov(find(AllCov(:,1)==1),3)));
fprintf(fid,'Mean Age in HC\t%g\n',mean(AllCov(find(AllCov(:,1)==0),3)));
fprintf(fid,'STD Age in HC\t%g\n',std(AllCov(find(AllCov(:,1)==0),3)));

fprintf(fid,'Number of White/caucasian\t%g\n',length(find(Race_Ethnicity==1)));
fprintf(fid,'Number of Black/African\t%g\n',length(find(Race_Ethnicity==2)));
fprintf(fid,'Number of Asian\t%g\n',length(find(Race_Ethnicity==3)));
fprintf(fid,'Number of Other Race Ethnicity\t%g\n',length(find(Race_Ethnicity==4)));
fprintf(fid,'Number of NaN Race Ethnicity\t%g\n',length(find(isnan(Race_Ethnicity))));
fprintf(fid,'Number of White/caucasian in Antidepressant Use MDD\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Black/African in Antidepressant Use MDD\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Asian in Antidepressant Use MDD\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==2))));
fprintf(fid,'Number of Other Race Ethnicity in Antidepressant Use MDD\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==2))));
fprintf(fid,'Number of NaN Race Ethnicity in Antidepressant Use MDD\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==2))));
fprintf(fid,'Number of White/caucasian in Antidepressant Free\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Black/African in Antidepressant Free\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Asian in Antidepressant Free\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==1))));
fprintf(fid,'Number of Other Race Ethnicity in Antidepressant Free\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==1))));
fprintf(fid,'Number of NaN Race Ethnicity in Antidepressant Free\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==1))));
fprintf(fid,'Number of White/caucasian in HC\t%g\n',length(find((Race_Ethnicity==1) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Black/African in HC\t%g\n',length(find((Race_Ethnicity==2) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Asian in HC\t%g\n',length(find((Race_Ethnicity==3) & (AllCov(:,1)==0))));
fprintf(fid,'Number of Other Race Ethnicity in HC\t%g\n',length(find((Race_Ethnicity==4) & (AllCov(:,1)==0))));
fprintf(fid,'Number of NaN Race Ethnicity in HC\t%g\n',length(find(isnan(Race_Ethnicity) & (AllCov(:,1)==0))));

fprintf(fid,'\n\n\n');
fclose(fid);


fprintf('\n\tPerform Stats for ENIGMA & REST-meta-MDD collaborative studies on vertex Thickness and Area within each site: 10 models: Finished!!!\n');


