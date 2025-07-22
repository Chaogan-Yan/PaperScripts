function ENIGMA_RMD_Vertex_Stats(DataDir,CovariatesFile,MinimumSubjectNumber,ParallelWorkersNumber,IsDPABISurfDataStyle,FreesurferSubjectsDir)
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


setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH'),':/usr/lib/fsl/5.0']);


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

if ~exist('IsDPABISurfDataStyle','var') || isempty(IsDPABISurfDataStyle)
    IsDPABISurfDataStyle=0;
end
if ischar(IsDPABISurfDataStyle)
    IsDPABISurfDataStyle=str2num(IsDPABISurfDataStyle);
end

if ~exist('FreesurferSubjectsDir','var') || isempty(FreesurferSubjectsDir)
    FreesurferSubjectsDir='/opt/freesurfer/subjects';
end


%1. Read the Covariates.xlsx -- % SHOULD have a column of eTIV!!!
if ischar(CovariatesFile)
    Table=readtable(CovariatesFile);

    %Deal with 'NA'
    for iCol=1:length(Table.Properties.VariableNames)
        HasEverNaN=0;
        eval(['Temp=Table.',Table.Properties.VariableNames{iCol},';']);
        for kk=1:size(Temp,1)
            if isequal(Temp(kk,1),{'NA'}) || isequal(Temp(kk,1),{''})
                Temp{kk,1}=NaN;
                HasEverNaN=1;
            end
        end
        if HasEverNaN
            Temp=cell2mat(Temp);
            eval(['Table.',Table.Properties.VariableNames{iCol},'=Temp;']);
        end
    end

else
    Table=CovariatesFile;
end


%2. Smooth the surfaces

[DPABIPath, fileN, extn] = fileparts(which('DPABI.m'));

SubjectIDString=[];
for i=1:length(Table.SubjID)
    SubjectIDString = sprintf('%s %s',SubjectIDString,Table.SubjID{i});
end


HasFreesurfer = system('which mri_surf2surf'); %Test if Freesurfer installed
if HasFreesurfer == 0
    %If yes, call Freesurfer's mri_surf2surf directly
    CommandInit=sprintf('export SUBJECTS_DIR=%s && ', FreesurferSubjectsDir);
    DataDirStr=DataDir;
else
    CommandInit=sprintf('docker run -ti --rm -v %s:/opt/freesurfer/license.txt -v %s:/data -e SUBJECTS_DIR=/opt/freesurfer/subjects cgyan/dpabi', fullfile(DPABIPath, 'DPABISurf', 'FreeSurferLicense', 'license.txt'), DataDir); 
    DataDirStr='/data';
end

if IsDPABISurfDataStyle
    MeasureStringSuffix='/fsaverage';
else
    MeasureStringSuffix=[];
end

Measure='Thickness';
MeasureString = [Measure,MeasureStringSuffix];
MeasureStringLower = lower(Measure);

Command = sprintf('%s parallel -j %g mri_surf2surf --s fsaverage --hemi lh --sval %s/AnatSurfLH/%s/{1}_space-fsaverage_hemi-L.%s.gii  --fwhm 10 --cortex --tval %s/AnatSurfLH/%s/s{1}_space-fsaverage_hemi-L.%s.gii ::: %s', CommandInit, ParallelWorkersNumber, DataDirStr, MeasureString, MeasureStringLower, DataDirStr, MeasureString, MeasureStringLower, SubjectIDString);
system(Command);
Command = sprintf('%s parallel -j %g mri_surf2surf --s fsaverage --hemi rh --sval %s/AnatSurfRH/%s/{1}_space-fsaverage_hemi-R.%s.gii  --fwhm 10 --cortex --tval %s/AnatSurfRH/%s/s{1}_space-fsaverage_hemi-R.%s.gii ::: %s', CommandInit, ParallelWorkersNumber, DataDirStr, MeasureString, MeasureStringLower, DataDirStr, MeasureString, MeasureStringLower, SubjectIDString);
system(Command);


Measure='Area';
MeasureString = [Measure,MeasureStringSuffix];
MeasureStringLower = lower(Measure);

Command = sprintf('%s parallel -j %g mri_surf2surf --s fsaverage --hemi lh --sval %s/AnatSurfLH/%s/{1}_space-fsaverage_hemi-L.%s.gii  --fwhm 10 --cortex --tval %s/AnatSurfLH/%s/s{1}_space-fsaverage_hemi-L.%s.gii ::: %s', CommandInit, ParallelWorkersNumber, DataDirStr, MeasureString, MeasureStringLower, DataDirStr, MeasureString, MeasureStringLower, SubjectIDString);
system(Command);
Command = sprintf('%s parallel -j %g mri_surf2surf --s fsaverage --hemi rh --sval %s/AnatSurfRH/%s/{1}_space-fsaverage_hemi-R.%s.gii  --fwhm 10 --cortex --tval %s/AnatSurfRH/%s/s{1}_space-fsaverage_hemi-R.%s.gii ::: %s', CommandInit, ParallelWorkersNumber, DataDirStr, MeasureString, MeasureStringLower, DataDirStr, MeasureString, MeasureStringLower, SubjectIDString);
system(Command);


%3. Do Stats Analysis for ALL subjects.
OutDir=[DataDir,'/Stats/AllSubjects'];
mkdir(OutDir);
ENIGMA_RMD_Vertex_Models(DataDir,OutDir,Table,MinimumSubjectNumber,IsDPABISurfDataStyle);

%4. Do Stats Analysis for Adolescent subjects.
OutDir=[DataDir,'/Stats/AdolescentSubjects'];
mkdir(OutDir);

TableAdolescent = Table(find(Table.Age<=21),:);

ENIGMA_RMD_Vertex_Models(DataDir,OutDir,TableAdolescent,MinimumSubjectNumber,IsDPABISurfDataStyle);


%5. Do Stats Analysis for Adult subjects.
OutDir=[DataDir,'/Stats/AdultSubjects'];
mkdir(OutDir);

TableAdult = Table(find(Table.Age>21),:);

ENIGMA_RMD_Vertex_Models(DataDir,OutDir,TableAdult,MinimumSubjectNumber,IsDPABISurfDataStyle);


fprintf('\n\tCongratulations! Perform ALL Stats for ENIGMA & REST-meta-MDD collaborative studies on vertex Thickness and Area within each site: 10 models: Finished!!!\n');


