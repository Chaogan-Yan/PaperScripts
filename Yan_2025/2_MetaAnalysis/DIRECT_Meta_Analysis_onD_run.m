function DIRECT_Meta_Analysis_onD_run(DataDir,OutputDir,AllSiteNames,SubjectsPool, Measure)


FDRQ=0.05;
FDRSuffix='_FDRCorrected';


[DPABIPath, fileN, extn] = fileparts(which('DPABI.m'));


%1. Model 1: Dx

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

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end

end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M1_Dx_Meta.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M1_Dx_Meta.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end







%2. Model 2: Dx by Age

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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M2_DxByAge.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M2_DxByAge.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M2_DxByAge.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M2_DxByAge.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end









%3. Model 3: Dx by Sex

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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M3_DxBySex.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M3_DxBySex.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(2)];
         N2=[N2;Info.Var2(3)];

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M3_DxBySex.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M3_DxBySex.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end








%4. Model 4: Recurrence Status

%4.2. Posthoc for Recurrence status: RecurrentMDD vs. FirstEpisodeMDD (2 vs. 1)
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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M4_Recur_PostHoc_RecurrentVsFirstEpisode.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M4_Recur_PostHoc_RecurrentVsFirstEpisode.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(33)];
         N2=[N2;Info.Var2(34)];

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M4_Recur_PostHoc_RecurrentVsFirstEpisode.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M4_Recur_PostHoc_RecurrentVsFirstEpisode.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end






%4.3. Posthoc for Recurrence status: RecurrentMDD vs. HC (2 vs. 0)
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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M4_Recur_PostHoc_RecurrentVsHC.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M4_Recur_PostHoc_RecurrentVsHC.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(33)];
         N2=[N2;Info.Var2(35)];

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M4_Recur_PostHoc_RecurrentVsHC.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M4_Recur_PostHoc_RecurrentVsHC.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end





%4.4. Posthoc for Recurrence status: FirstEpisodeMDD vs. HC (1 vs. 0)
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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M4_Recur_PostHoc_FirstEpisodeVsHC.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M4_Recur_PostHoc_FirstEpisodeVsHC.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(34)];
         N2=[N2;Info.Var2(35)];

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M4_Recur_PostHoc_FirstEpisodeVsHC.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M4_Recur_PostHoc_FirstEpisodeVsHC.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end




%5. Model 5: Remission status


%5.2. Posthoc for Remission status: AcutelyDepressed vs. RemittedeMDD (2 vs. 1)
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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M5_Rem_PostHoc_AcutelyDepressedVsRemittedeMDD.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M5_Rem_PostHoc_AcutelyDepressedVsRemittedeMDD.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(74)];
         N2=[N2;Info.Var2(75)];

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M5_Rem_PostHoc_AcutelyDepressedVsRemittedeMDD.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M5_Rem_PostHoc_AcutelyDepressedVsRemittedeMDD.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end




%5.3. Posthoc for Remission status: AcutelyDepressed vs. HC (2 vs. 0)
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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M5_Rem_PostHoc_AcutelyDepressedVsHC.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M5_Rem_PostHoc_AcutelyDepressedVsHC.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(74)];
         N2=[N2;Info.Var2(76)];

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M5_Rem_PostHoc_AcutelyDepressedVsHC.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M5_Rem_PostHoc_AcutelyDepressedVsHC.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end






%5.4. Posthoc for Remission status: RemittedeMDD vs. HC (1 vs. 0)
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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M5_Rem_PostHoc_RemittedeMDDVsHC.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M5_Rem_PostHoc_RemittedeMDDVsHC.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(75)];
         N2=[N2;Info.Var2(76)];

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M5_Rem_PostHoc_RemittedeMDDVsHC.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M5_Rem_PostHoc_RemittedeMDDVsHC.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end








%6. Model 6: Age of Onset

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
    

    if strcmpi(Measure,'Area')
        TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M6_AgeOfOnsetgii.gii'];
    else
        TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M6_AgeOfOnset.gii'];
    end
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M6_AgeOfOnset.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(114)];

         if N1(end)==0
             N1(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M6_AgeOfOnset.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M6_AgeOfOnset.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end






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
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsAdultAO.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         
         try
             N1=[N1;Info.Var2(236)];
             N2=[N2;Info.Var2(237)];
         catch
             N1=[N1;Info.Var2(128)];
             N2=[N2;Info.Var2(129)];
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
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsHC.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdolescentAOVsHC.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         
         
         try
             N1=[N1;Info.Var2(236)];
             N2=[N2;Info.Var2(238)];
         catch
             N1=[N1;Info.Var2(128)];
             N2=[N2;Info.Var2(130)];
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
    else
        Style_fsaverage = [];
    end
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdultAOVsHC.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M7_AOGroup_PostHoc_AdultAOVsHC.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         
         try
             N1=[N1;Info.Var2(237)];
             N2=[N2;Info.Var2(238)];
         catch
             N1=[N1;Info.Var2(129)];
             N2=[N2;Info.Var2(130)];
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








%8. Model 8: BDI severity of depressive symptoms

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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M8_BDI.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M8_BDI.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(168)];

         if N1(end)==0
             N1(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M8_BDI.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M8_BDI.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end







%9. Model 9: HDRS severity of depressive symptoms

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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M9_HDRS.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M9_HDRS.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(181)];

         if N1(end)==0
             N1(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[]; 
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M9_HDRS.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M9_HDRS.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end






%10. Model 10: Antidepressant medication use at time of scan

%10.2. Posthoc for Antidepressant status: AntidepressantUse vs. AntidepressantFree (2 vs. 1)
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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M10_Antidepressant_PostHoc_AntidepressantUseVsAntidepressantFree.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M10_Antidepressant_PostHoc_AntidepressantUseVsAntidepressantFree.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(195)];
         N2=[N2;Info.Var2(196)];

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M10_Antidepressant_PostHoc_AntidepressantUseVsAntidepressantFree.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M10_Antidepressant_PostHoc_AntidepressantUseVsAntidepressantFree.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end





%10.3. Posthoc for Antidepressant status: AntidepressantUse vs. HC (2 vs. 0)
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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M10_Antidepressant_PostHoc_AntidepressantUseVsHC.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M10_Antidepressant_PostHoc_AntidepressantUseVsHC.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(195)];
         N2=[N2;Info.Var2(197)];

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[];
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M10_Antidepressant_PostHoc_AntidepressantUseVsHC.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M10_Antidepressant_PostHoc_AntidepressantUseVsHC.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end




%10.4. Posthoc for Antidepressant status: AntidepressantFree vs. HC (1 vs. 0)
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
    

     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,Style_fsaverage,filesep,'M10_Antidepressant_PostHoc_AntidepressantFreeVsHC.gii'];
     if exist(TempFile)
         TFilesLH=[TFilesLH;{TempFile}];
     end
     
     TempFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,Style_fsaverage,filesep,'M10_Antidepressant_PostHoc_AntidepressantFreeVsHC.gii'];
     if exist(TempFile)
         TFilesRH=[TFilesRH;{TempFile}];
     end
     
     InfoFile=[DataDir,filesep,SiteName,filesep,'/Stats/',SubjectsPool,filesep,'ModelSummaryInfo.txt'];
     
     Info=readtable(InfoFile);
     
     if exist(TempFile)
         N1=[N1;Info.Var2(196)];
         N2=[N2;Info.Var2(197)];

         if (N1(end)==0) || (N2(end)==0)
             N1(end)=[];
             N2(end)=[];
             TFilesLH(end)=[];
             TFilesRH(end)=[]; 
         end
     end
end

if length(N1)>=3
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure])
    OutputNameLH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfLH',filesep,Measure,filesep,'M10_Antidepressant_PostHoc_AntidepressantFreeVsHC.gii'];
    MaskFileLH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_lh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesLH,OutputNameLH,MaskFileLH,N1, N2, Regressor);
    
    
    mkdir([OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure])
    OutputNameRH=[OutputDir,filesep,SubjectsPool,filesep,'AnatSurfRH',filesep,Measure,filesep,'M10_Antidepressant_PostHoc_AntidepressantFreeVsHC.gii'];
    MaskFileRH=fullfile(DPABIPath, 'DPABISurf', 'SurfTemplates','fsaverage_rh_cortex.label.gii');
    y_Meta_Image_CallR(TFilesRH,OutputNameRH,MaskFileRH,N1, N2, Regressor);
    
    y_FDR_SurfaceLHRH_WithP_onD(OutputNameLH,OutputNameRH,MaskFileLH,MaskFileRH,FDRQ,FDRSuffix);
end



fprintf('\n\tPerform Meta Analysis for ENIGMA & REST-meta-MDD collaborative studies on vertex: 10 models: Finished!!!\n');


