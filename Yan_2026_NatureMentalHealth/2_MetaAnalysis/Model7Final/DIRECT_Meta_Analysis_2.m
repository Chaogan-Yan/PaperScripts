
%Please specify:

%DataDir - Where the meta maps. E.g., Site1/Stats/AllSubjects/AnatSurfRH/Area/fsaverage/M1_Dx.gii under DataDir

%OutputDir - Where the final results.

%AllSiteNames - All the site names. E.g., {'Site1';'Site2';'Site3';'Site4';'Site5'}






%For All Sites Corrected on Cohen's D maps
DataDir='/mnt/Data7/RfMRILab/Yan/YAN_Work/REST_meta-MDD_Surf/Analysis/MetaAnalysis/AfterQCExclusion/FixModel7/DIRECT_ENIGMAData';

OutputDir='/mnt/Data7/RfMRILab/Yan/YAN_Work/REST_meta-MDD_Surf/Analysis/MetaAnalysis/AfterQCExclusion/FixModel7/MetaAnalysis_DIRECT_ENIGMA_AllSites_CohenDMaps_RemoveEmptySites';


% % In y_Meta_Image_CallR.m Add:
% NeedRemoveI=[];
% for iFile=1:size(AllVolume,1)
%     if ~any(AllVolume(iFile,:))
%         NeedRemoveI=[NeedRemoveI;iFile];
%     end
% end
% AllVolume(NeedRemoveI,:)=[];
% N1(NeedRemoveI,:)=[];
% if length(N2)>2
%     N2(NeedRemoveI,:)=[];
% end
% if length(Regressor)>2
%     Regressor(NeedRemoveI,:)=[];
% end
% TFiles(NeedRemoveI,:)=[];
% nDimTimePoints=size(AllVolume,1);
% 



Dir=dir(DataDir);
AllSiteNames={};
for i=3:length(Dir)
    if Dir(i).isdir
        AllSiteNames=[AllSiteNames;{Dir(i).name}];
    end
end


SubjectsPoolSet={'AllSubjects';'AdolescentSubjects';'AdultSubjects'};
%SubjectsPoolSet={'AdolescentSubjects';'AdultSubjects'};
MeasureSet={'Thickness';'Area'};

for iSubjectsPool = 1:length(SubjectsPoolSet)
    for iMeasure = 1:length(MeasureSet)
        
        SubjectsPool=SubjectsPoolSet{iSubjectsPool};
        Measure=MeasureSet{iMeasure};
        
        DIRECT_Meta_Analysis_onD_run(DataDir,OutputDir,AllSiteNames,SubjectsPool, Measure);
    end
end


fprintf('\n\tPerform All Meta Analysis for ENIGMA & REST-meta-MDD collaborative studies: Finished!!!\n');







