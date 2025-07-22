function DIRECT_Meta_Regression_Model1(DataDir,OutputDir,AllSiteNames,SubjectsPool, Measure, SiteInfo)

FDRQ=0.05;
FDRSuffix='_FDRCorrected';


[DPABIPath, fileN, extn] = fileparts(which('DPABI.m'));


%1. Meta Regression 1: Mean Age

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:length(AllSiteNames)
    SiteName=AllSiteNames{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         %Rate = (Info.Var2(23)+Info.Var2(28))/sum(Info.Var2(21:30));
         Regressor=[Regressor;Info.Var2(10)];
         
         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression1_MeanAge.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,fiqlesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression1_MeanAge.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end







%2. Meta Regression 2: ProportionFemale

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:length(AllSiteNames)
    SiteName=AllSiteNames{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Rate = (Info.Var2(5))/sum(Info.Var2(4:5));
         Regressor=[Regressor;Rate];
         
         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression2_ProportionFemale.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression2_ProportionFemale.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end







%3. Meta Regression 3: ProportionAsian

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:length(AllSiteNames)
    SiteName=AllSiteNames{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Rate = (Info.Var2(23)+Info.Var2(28))/sum(Info.Var2(21:30));
         Regressor=[Regressor;Rate];
         
         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression3_ProportionAsian.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression3_ProportionAsian.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end





%4. Meta Regression 4: ProportionCaucasian

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:length(AllSiteNames)
    SiteName=AllSiteNames{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Rate = (Info.Var2(21)+Info.Var2(26))/sum(Info.Var2(21:30));
         Regressor=[Regressor;Rate];
         
         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression4_ProportionCaucasian.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression4_ProportionCaucasian.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end





%5. Meta Regression 5: ProportionAfrican

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:length(AllSiteNames)
    SiteName=AllSiteNames{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Rate = (Info.Var2(22)+Info.Var2(27))/sum(Info.Var2(21:30));
         Regressor=[Regressor;Rate];
         
         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression5_ProportionAfrican.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression5_ProportionAfrican.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end










%6. Meta Regression 6: MeanHAMD

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:size(SiteInfo,1)
    SiteName=SiteInfo.Site{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Metric=SiteInfo.MeanHAMD(iSite);
         Regressor=[Regressor;Metric];
         
         if (N1(end)==0) || (N2(end)==0) || isnan(Metric)
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression6_MeanHAMD.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression6_MeanHAMD.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end






%7. Meta Regression 7: ProportionRecurrent

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:size(SiteInfo,1)
    SiteName=SiteInfo.Site{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Metric=SiteInfo.ProportionRecurrent(iSite);
         Regressor=[Regressor;Metric];
         
         if (N1(end)==0) || (N2(end)==0) || isnan(Metric)
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression7_ProportionRecurrent.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression7_ProportionRecurrent.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end






%8. Meta Regression 8: ProportionAntidepressant

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:size(SiteInfo,1)
    SiteName=SiteInfo.Site{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Metric=SiteInfo.ProportionAntidepressant(iSite);
         Regressor=[Regressor;Metric];
         
         if (N1(end)==0) || (N2(end)==0) || isnan(Metric)
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression8_ProportionAntidepressant.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression8_ProportionAntidepressant.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end






%9. Meta Regression 9: FreeSurferVersion (6 or not)

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:size(SiteInfo,1)
    SiteName=SiteInfo.Site{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Metric=SiteInfo.FreeSurferVersion(iSite)==6;
         Regressor=[Regressor;Metric];
         
         if (N1(end)==0) || (N2(end)==0) || isnan(SiteInfo.FreeSurferVersion(iSite))
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression9_FreeSurferVersion.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression9_FreeSurferVersion.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end






%10. Meta Regression 10: ScannerSIEMENS (SIEMENS or not)

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:size(SiteInfo,1)
    SiteName=SiteInfo.Site{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Metric=SiteInfo.Scanner(iSite)==1;
         Regressor=[Regressor;Metric];
         
         if (N1(end)==0) || (N2(end)==0) || isnan(SiteInfo.Scanner(iSite))
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression10_ScannerSIEMENS.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression10_ScannerSIEMENS.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end





%10. Meta Regression 11: ScannerSIEMENSvsGE (SIEMENS or GE)

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:size(SiteInfo,1)
    SiteName=SiteInfo.Site{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Metric=SiteInfo.Scanner(iSite)==1;
         Regressor=[Regressor;Metric];
         
         if (N1(end)==0) || (N2(end)==0) || SiteInfo.Scanner(iSite)==2
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression11_ScannerSIEMENSvsGE.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression11_ScannerSIEMENSvsGE.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end






%12. Meta Regression 12: FieldStrength (3T or not)

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:size(SiteInfo,1)
    SiteName=SiteInfo.Site{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Metric=SiteInfo.FieldStrength(iSite)==3;
         Regressor=[Regressor;Metric];
         
         if (N1(end)==0) || (N2(end)==0) || isnan(Metric)
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression12_FieldStrength.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression12_FieldStrength.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end





%13. Meta Regression 13: SampleType (clinical or not)

TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:size(SiteInfo,1)
    SiteName=SiteInfo.Site{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M1_Dx.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         Metric=SiteInfo.SampleType(iSite)==1;
         Regressor=[Regressor;Metric];
         
         if (N1(end)==0) || (N2(end)==0) || isnan(SiteInfo.SampleType(iSite))
             N1(end)=[];
             N2(end)=[];
             Regressor(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if (length(N1)>=3) && (std(Regressor)~=0)
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_MetaRegression13_SampleType.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_MetaRegression13_SampleType.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end








fprintf('\n\tPerform Meta regression for ENIGMA & REST-meta-MDD collaborative studies on vertex: Finished!!!\n');


