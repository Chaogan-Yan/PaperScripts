%Calculate Stability on the surface level
%%%%%%%%%%%%%%%%%%%% main settings%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
DataDir = '/mnt/Data3/RfMRILab/ChenX/Rumination_surf';
OutputDir = '/mnt/Data3/RfMRILab/ChenX/Rumination_surf/Analysis/Stability_Surf/smoothed';
TempDir = '/mnt/Data3/RfMRILab/ChenX/CX_software/DPABI_V4.2_190919/DPABISurf/SurfTemplates';
SiteSet  = {'IPCAS','PKUGE','PKUSIEMENS'};

load([DataDir,'/Analysis/DemographicInfo_V2.mat']);
WindowSize = 32; WindowStep = 2; WindowType = 'rectwin';
InFile_Volume = [];ROIDef = 'VertexToVertex';
AMaskFilename_LH = [TempDir,'/fsaverage5_lh_cortex.label.gii'];
AMaskFilename_RH = [TempDir,'/fsaverage5_rh_cortex.label.gii'];
IsMultipleLabel = 0; IsNeedDetrend = 0; CUTNUMBER = 1;

parpool(7);

%rest
for iSite = 1:3
     parfor iSubject = 1:length(Sublist)
         OutputName_LH = [OutputDir,'/',SiteSet{iSite},'/rest/',Sublist{iSubject},'-LH'];
         OutputName_RH = [OutputDir,'/',SiteSet{iSite},'/rest/',Sublist{iSubject},'-RH'];
         mkdir([OutputDir,'/',SiteSet{iSite},'/rest/']);
         InFile_LH = [DataDir,'/',SiteSet{iSite},'_rest/FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_task-rest_space-fsaverage5_hemi-L.func.gii'];
         InFile_RH = [DataDir,'/',SiteSet{iSite},'_rest/FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_task-rest_space-fsaverage5_hemi-R.func.gii'];
         [~, ~, ~, ~] = y_Stability_Surf_Window(WindowSize, WindowStep, WindowType, InFile_LH, InFile_RH, InFile_Volume, ROIDef, OutputName_LH, OutputName_RH, AMaskFilename_LH, AMaskFilename_RH, IsMultipleLabel, IsNeedDetrend, CUTNUMBER);
     end  
end

%rum
for iSite = 1:3
     parfor iSubject = 1:length(Sublist)
         OutputName_LH = [OutputDir,'/',SiteSet{iSite},'/rum/',Sublist{iSubject},'-LH'];
         OutputName_RH = [OutputDir,'/',SiteSet{iSite},'/rum/',Sublist{iSubject},'-RH'];
         mkdir([OutputDir,'/',SiteSet{iSite},'/rum/']);
         InFile_LH = [DataDir,'/',SiteSet{iSite},'_task/S3_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-3_task-rest_space-fsaverage5_hemi-L.func.gii'];
         InFile_RH = [DataDir,'/',SiteSet{iSite},'_task/S3_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-3_task-rest_space-fsaverage5_hemi-R.func.gii'];
         [~, ~, ~, ~] = y_Stability_Surf_Window(WindowSize, WindowStep, WindowType, InFile_LH, InFile_RH, InFile_Volume, ROIDef, OutputName_LH, OutputName_RH, AMaskFilename_LH, AMaskFilename_RH, IsMultipleLabel, IsNeedDetrend, CUTNUMBER);
     end  
end

%dis
for iSite = 1:3
     parfor iSubject = 1:length(Sublist)
         OutputName_LH = [OutputDir,'/',SiteSet{iSite},'/dis/',Sublist{iSubject},'-LH'];
         OutputName_RH = [OutputDir,'/',SiteSet{iSite},'/dis/',Sublist{iSubject},'-RH'];
         mkdir([OutputDir,'/',SiteSet{iSite},'/dis/']);
         InFile_LH = [DataDir,'/',SiteSet{iSite},'_task/S4_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-4_task-rest_space-fsaverage5_hemi-L.func.gii'];
         InFile_RH = [DataDir,'/',SiteSet{iSite},'_task/S4_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-4_task-rest_space-fsaverage5_hemi-R.func.gii'];
         [~, ~, ~, ~] = y_Stability_Surf_Window(WindowSize, WindowStep, WindowType, InFile_LH, InFile_RH, InFile_Volume, ROIDef, OutputName_LH, OutputName_RH, AMaskFilename_LH, AMaskFilename_RH, IsMultipleLabel, IsNeedDetrend, CUTNUMBER);
     end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WL = 16 TRs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
DataDir = '/mnt/Data3/RfMRILab/ChenX/Rumination_surf';
OutputDir = '/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Stability_Surf/smoothed_WL_16';
TempDir = '/mnt/Data3/RfMRILab/ChenX/CX_software/DPABI_V4.2_190919/DPABISurf/SurfTemplates';
SiteSet  = {'IPCAS','PKUGE','PKUSIEMENS'};

load([DataDir,'/Analysis/DemographicInfo_V2.mat']);
WindowSize = 16; WindowStep = 2; WindowType = 'rectwin';
InFile_Volume = [];ROIDef = 'VertexToVertex';
AMaskFilename_LH = [TempDir,'/fsaverage5_lh_cortex.label.gii'];
AMaskFilename_RH = [TempDir,'/fsaverage5_rh_cortex.label.gii'];
IsMultipleLabel = 0; IsNeedDetrend = 0; CUTNUMBER = 1;

parpool(11);

%rest
for iSite = 1:3
     parfor iSubject = 1:length(Sublist)
         OutputName_LH = [OutputDir,'/',SiteSet{iSite},'/rest/',Sublist{iSubject},'-LH'];
         OutputName_RH = [OutputDir,'/',SiteSet{iSite},'/rest/',Sublist{iSubject},'-RH'];
         mkdir([OutputDir,'/',SiteSet{iSite},'/rest/']);
         InFile_LH = [DataDir,'/',SiteSet{iSite},'_rest/FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_task-rest_space-fsaverage5_hemi-L.func.gii'];
         InFile_RH = [DataDir,'/',SiteSet{iSite},'_rest/FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_task-rest_space-fsaverage5_hemi-R.func.gii'];
         [~, ~, ~, ~] = y_Stability_Surf_Window(WindowSize, WindowStep, WindowType, InFile_LH, InFile_RH, InFile_Volume, ROIDef, OutputName_LH, OutputName_RH, AMaskFilename_LH, AMaskFilename_RH, IsMultipleLabel, IsNeedDetrend, CUTNUMBER);
     end  
end

%rum
for iSite = 1:3
     parfor iSubject = 1:length(Sublist)
         OutputName_LH = [OutputDir,'/',SiteSet{iSite},'/rum/',Sublist{iSubject},'-LH'];
         OutputName_RH = [OutputDir,'/',SiteSet{iSite},'/rum/',Sublist{iSubject},'-RH'];
         mkdir([OutputDir,'/',SiteSet{iSite},'/rum/']);
         InFile_LH = [DataDir,'/',SiteSet{iSite},'_task/S3_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-3_task-rest_space-fsaverage5_hemi-L.func.gii'];
         InFile_RH = [DataDir,'/',SiteSet{iSite},'_task/S3_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-3_task-rest_space-fsaverage5_hemi-R.func.gii'];
         [~, ~, ~, ~] = y_Stability_Surf_Window(WindowSize, WindowStep, WindowType, InFile_LH, InFile_RH, InFile_Volume, ROIDef, OutputName_LH, OutputName_RH, AMaskFilename_LH, AMaskFilename_RH, IsMultipleLabel, IsNeedDetrend, CUTNUMBER);
     end  
end

%dis
for iSite = 1:3
     parfor iSubject = 1:length(Sublist)
         OutputName_LH = [OutputDir,'/',SiteSet{iSite},'/dis/',Sublist{iSubject},'-LH'];
         OutputName_RH = [OutputDir,'/',SiteSet{iSite},'/dis/',Sublist{iSubject},'-RH'];
         mkdir([OutputDir,'/',SiteSet{iSite},'/dis/']);
         InFile_LH = [DataDir,'/',SiteSet{iSite},'_task/S4_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-4_task-rest_space-fsaverage5_hemi-L.func.gii'];
         InFile_RH = [DataDir,'/',SiteSet{iSite},'_task/S4_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-4_task-rest_space-fsaverage5_hemi-R.func.gii'];
         [~, ~, ~, ~] = y_Stability_Surf_Window(WindowSize, WindowStep, WindowType, InFile_LH, InFile_RH, InFile_Volume, ROIDef, OutputName_LH, OutputName_RH, AMaskFilename_LH, AMaskFilename_RH, IsMultipleLabel, IsNeedDetrend, CUTNUMBER);
     end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WL = 64 TRs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
DataDir = '/mnt/Data3/RfMRILab/ChenX/Rumination_surf';
OutputDir = '/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Stability_Surf/smoothed_WL_64';
TempDir = '/mnt/Data3/RfMRILab/ChenX/CX_software/DPABI_V4.2_190919/DPABISurf/SurfTemplates';
SiteSet  = {'IPCAS','PKUGE','PKUSIEMENS'};

load([DataDir,'/Analysis/DemographicInfo_V2.mat']);
WindowSize = 64; WindowStep = 2; WindowType = 'rectwin';
InFile_Volume = [];ROIDef = 'VertexToVertex';
AMaskFilename_LH = [TempDir,'/fsaverage5_lh_cortex.label.gii'];
AMaskFilename_RH = [TempDir,'/fsaverage5_rh_cortex.label.gii'];
IsMultipleLabel = 0; IsNeedDetrend = 0; CUTNUMBER = 1;

parpool(11);

%rest
for iSite = 1:3
     parfor iSubject = 1:length(Sublist)
         OutputName_LH = [OutputDir,'/',SiteSet{iSite},'/rest/',Sublist{iSubject},'-LH'];
         OutputName_RH = [OutputDir,'/',SiteSet{iSite},'/rest/',Sublist{iSubject},'-RH'];
         mkdir([OutputDir,'/',SiteSet{iSite},'/rest/']);
         InFile_LH = [DataDir,'/',SiteSet{iSite},'_rest/FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_task-rest_space-fsaverage5_hemi-L.func.gii'];
         InFile_RH = [DataDir,'/',SiteSet{iSite},'_rest/FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_task-rest_space-fsaverage5_hemi-R.func.gii'];
         [~, ~, ~, ~] = y_Stability_Surf_Window(WindowSize, WindowStep, WindowType, InFile_LH, InFile_RH, InFile_Volume, ROIDef, OutputName_LH, OutputName_RH, AMaskFilename_LH, AMaskFilename_RH, IsMultipleLabel, IsNeedDetrend, CUTNUMBER);
     end  
end

%rum
for iSite = 1:3
     parfor iSubject = 1:length(Sublist)
         OutputName_LH = [OutputDir,'/',SiteSet{iSite},'/rum/',Sublist{iSubject},'-LH'];
         OutputName_RH = [OutputDir,'/',SiteSet{iSite},'/rum/',Sublist{iSubject},'-RH'];
         mkdir([OutputDir,'/',SiteSet{iSite},'/rum/']);
         InFile_LH = [DataDir,'/',SiteSet{iSite},'_task/S3_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-3_task-rest_space-fsaverage5_hemi-L.func.gii'];
         InFile_RH = [DataDir,'/',SiteSet{iSite},'_task/S3_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-3_task-rest_space-fsaverage5_hemi-R.func.gii'];
         [~, ~, ~, ~] = y_Stability_Surf_Window(WindowSize, WindowStep, WindowType, InFile_LH, InFile_RH, InFile_Volume, ROIDef, OutputName_LH, OutputName_RH, AMaskFilename_LH, AMaskFilename_RH, IsMultipleLabel, IsNeedDetrend, CUTNUMBER);
     end  
end

%dis
for iSite = 1:3
     parfor iSubject = 1:length(Sublist)
         OutputName_LH = [OutputDir,'/',SiteSet{iSite},'/dis/',Sublist{iSubject},'-LH'];
         OutputName_RH = [OutputDir,'/',SiteSet{iSite},'/dis/',Sublist{iSubject},'-RH'];
         mkdir([OutputDir,'/',SiteSet{iSite},'/dis/']);
         InFile_LH = [DataDir,'/',SiteSet{iSite},'_task/S4_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-4_task-rest_space-fsaverage5_hemi-L.func.gii'];
         InFile_RH = [DataDir,'/',SiteSet{iSite},'_task/S4_FunSurfWCFS/sub-',Sublist{iSubject},'/ssub-',Sublist{iSubject},'_ses-4_task-rest_space-fsaverage5_hemi-R.func.gii'];
         [~, ~, ~, ~] = y_Stability_Surf_Window(WindowSize, WindowStep, WindowType, InFile_LH, InFile_RH, InFile_Volume, ROIDef, OutputName_LH, OutputName_RH, AMaskFilename_LH, AMaskFilename_RH, IsMultipleLabel, IsNeedDetrend, CUTNUMBER);
     end  
end