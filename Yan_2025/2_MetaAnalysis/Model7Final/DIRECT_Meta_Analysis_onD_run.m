function DIRECT_Meta_Analysis_onD_run(DataDir,OutputDir,AllSiteNames,SubjectsPool, Measure)


FDRQ=0.05;
FDRSuffix='_FDRCorrected';


[DPABIPath, fileN, extn] = fileparts(which('DPABI.m'));



%7. Model 7: Age of Onset ANOVA analysis!


%7.2. Posthoc for Age of Onset status: AdolescentAO vs. AdultAO (2 vs. 1)
TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:length(AllSiteNames)
    SiteName=AllSiteNames{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
        Style_Stats='Stats';
    else
        Style_fsaverage = [];
        Style_Stats='StatsFixModel7';
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/',Style_Stats,'/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/',Style_Stats,'/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/',Style_Stats,'/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         
         try
             N1=[N1;Info.Var2(2)];
             N2=[N2;Info.Var2(3)];
         catch
         end

         if (N1(end)<=3) || (N2(end)<=3)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end




%7.3. Posthoc for Age of Onset status: AdolescentAO vs. HC (2 vs. 0)
TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:length(AllSiteNames)
    SiteName=AllSiteNames{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
        Style_Stats='Stats';
    else
        Style_fsaverage = [];
        Style_Stats='StatsFixModel7';
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/',Style_Stats,'/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsHC.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/',Style_Stats,'/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsHC.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/',Style_Stats,'/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         
         
         try
             N1=[N1;Info.Var2(2)];
             N2=[N2;Info.Var2(4)];
         catch
         end

         if (N1(end)<=3) || (N2(end)<=3)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsHC.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsHC.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end







%7.4. Posthoc for Age of Onset status: AdultAO vs. HC (1 vs. 0)
TFilesLH={};
TFilesRH={};

N1=[];
N2=[];
Regressor=[];
for iSite=1:length(AllSiteNames)
    SiteName=AllSiteNames{iSite};
    
    if (length(SiteName) >= 3) && strcmp(SiteName(1:3),'IS0')
        Style_fsaverage = [filesep,'fsaverage'];
        Style_Stats='Stats';
    else
        Style_fsaverage = [];
        Style_Stats='StatsFixModel7';
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/',Style_Stats,'/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdultAOVsHC.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/',Style_Stats,'/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdultAOVsHC.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/',Style_Stats,'/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         
         try
             N1=[N1;Info.Var2(3)];
             N2=[N2;Info.Var2(4)];
         catch
         end

         if (N1(end)<=3) || (N2(end)<=3)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M7_AOGroup_PostHoc_AdultAOVsHC.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M7_AOGroup_PostHoc_AdultAOVsHC.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end








fprintf('\n\tPerform Meta Analysis for ENIGMA & REST-meta-MDD collaborative studies on vertex: Model 7: Finished!!!\n');


