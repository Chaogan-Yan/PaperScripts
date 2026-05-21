function DIRECT_Meta_Sites(DataDir,OutputDir,ETable2,MinimumSubjectNumber,ParallelWorkersNumber)
% FORMAT ENIGMA_RMD_Vertex_Stats(DataDir,CovariatesFile,MinimumSubjectNumber,ParallelWorkersNumber,IsDPABISurfDataStyle,FreesurferSubjectsDir)
% Perform Stats for ENIGMA & REST-meta-MDD collaborative studies on vertex Thickness and Area within each site: 10 models.
% !Environment preparation: 1) Add folder for SPM12 (NOT add with subfolders). 2) Add with subfolders for DPABI_V6_ForENIGMA_RMD and Scripts_ENIGMA_RMD; 3) pull DPABI docker: docker pull cgyan/dpabi, OR install freesurfer and GNU parallel on your machine.
% Input:
%   DataDir - Data Dir for stats. Should have AnatSurfLH AnatSurfRH under this dir.
%   CovariatesFile - Covariates File (.xlsx) in ENIGMA-MDD style.  SHOULD have a column of eTIV!!!
%   FreesurferSubjectsDir - Freesurfer dir for smoothing. Should have 'fsaverage' subject under this dir.
%   MinimumSubjectNumber - Minimum Subject Number in each group. If the number of availabe subjects are less than this number, then the stats will be skipped. default: 10
%   ParallelWorkersNumber - Parallel workers for smoothing surfaces. default: 1
%   IsDPABISurfDataStyle - If the data is directly come from DPABISurf, then 1. If the ENIGMA script was used, then 0. default: 0
%
% Output:
%   The stats under DataDir/Stats.
%___________________________________________________________________________
% Written by YAN Chao-Gan 210519.
% CAS Key Laboratory of Behavioral Science, Institute of Psychology, Beijing, China;
% International Big-Data Research Center for Depression (IBRCD), Institute of Psychology, Chinese Academy of Sciences, Beijing, China;
% Magnetic Resonance Imaging Research Center, Institute of Psychology, Chinese Academy of Sciences, Beijing, China.
% ycg.yan@gmail.com



if ~exist('MinimumSubjectNumber','var') || isempty(MinimumSubjectNumber)
    MinimumSubjectNumber=10;
end
if ischar(MinimumSubjectNumber)
    MinimumSubjectNumber=str2num(MinimumSubjectNumber);
end

if ~exist('ParallelWorkersNumber','var') || isempty(ParallelWorkersNumber)
    ParallelWorkersNumber=1;
end
if ischar(ParallelWorkersNumber)
    ParallelWorkersNumber=str2num(ParallelWorkersNumber);
end

IsDPABISurfDataStyle=1;

Site=ETable2.Site;
USite=unique(Site);

for iSite=1:length(USite)
    SiteName=sprintf('IS%.3d', USite(iSite));
    Table=ETable2(find(Site==USite(iSite)),:);



    %3. Do Stats Analysis for ALL subjects.
    OutDir=[OutputDir,'/',SiteName,'/Stats/AllSubjects'];
    mkdir(OutDir);
    ENIGMA_RMD_Vertex_Models(DataDir,OutDir,Table,MinimumSubjectNumber,IsDPABISurfDataStyle);

    %4. Do Stats Analysis for Adolescent subjects.
    OutDir=[OutputDir,'/',SiteName,'/Stats/AdolescentSubjects'];
    mkdir(OutDir);

    TableAdolescent = Table(find(Table.Age<=21),:);

    ENIGMA_RMD_Vertex_Models(DataDir,OutDir,TableAdolescent,MinimumSubjectNumber,IsDPABISurfDataStyle);


    %5. Do Stats Analysis for Adult subjects.
    OutDir=[OutputDir,'/',SiteName,'/Stats/AdultSubjects'];
    mkdir(OutDir);

    TableAdult = Table(find(Table.Age>21),:);

    ENIGMA_RMD_Vertex_Models(DataDir,OutDir,TableAdult,MinimumSubjectNumber,IsDPABISurfDataStyle);


end

fprintf('\n\tCongratulations! Perform ALL Stats for ENIGMA & REST-meta-MDD collaborative studies on vertex Thickness and Area within each site: 10 models: Finished!!!\n');


