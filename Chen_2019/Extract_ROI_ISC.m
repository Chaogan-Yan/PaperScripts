%Extract the network level ROIs' ISC

%initialize
clc;clear all;
DataDir = '/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis4Publish/ISC/ROI_ISC_CorrMaps';
SubList = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/Analysis/IPCAS_Sublist.txt');
ResultDir = '/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis4Publish/ISC';
SiteSet = {'IPCAS','PKUGE','PKUSIMENS'};
SessionNameSet = {'Rest','happy', 'sad', 'rum', 'dis'};

%extract each network's ISC (not ISFC)
for iSite = 1:length(SiteSet)
    ISCMap = {};
    for iSession = 1:length(SessionNameSet) 
        for iSubject = 1:length(SubList)
            MeanISCorrelationMap = [];
            load([DataDir,'/',SessionNameSet{iSession},'/',SiteSet{iSite},'/',SubList{iSubject},'_ISCMap.mat']);
            %extract DMN core
            iCount_temp=0;
            Sum_temp=0;
            for isub = 1:9
                  Sum_temp=Sum_temp+MeanISCorrelationMap(isub,isub);
                  iCount_temp=iCount_temp+1;
            end
            ISCMap{iSession}(iSubject,1) = Sum_temp/iCount_temp;
            %extract DMN dmPFC
            iCount_temp=0;
            Sum_temp=0;
            for isub = 10:18
                  Sum_temp=Sum_temp+MeanISCorrelationMap(isub,isub);
                  iCount_temp=iCount_temp+1;
            end
            ISCMap{iSession}(iSubject,2) = Sum_temp/iCount_temp;
            %extract DMN MTL
            iCount_temp=0;
            Sum_temp=0;
            for isub = 19:24
                  Sum_temp=Sum_temp+MeanISCorrelationMap(isub,isub);
                  iCount_temp=iCount_temp+1;
            end
            ISCMap{iSession}(iSubject,3) = Sum_temp/iCount_temp;
            %extract DMN FPCN
            iCount_temp=0;
            Sum_temp=0;
            for isub = 25:35
                  Sum_temp=Sum_temp+MeanISCorrelationMap(isub,isub);
                  iCount_temp=iCount_temp+1;
            end
            ISCMap{iSession}(iSubject,4) = Sum_temp/iCount_temp;
            %extract DMN SN
            iCount_temp=0;
            Sum_temp=0;
            for isub = 36:59
                  Sum_temp=Sum_temp+MeanISCorrelationMap(isub,isub);
                  iCount_temp=iCount_temp+1;
            end
            ISCMap{iSession}(iSubject,5) = Sum_temp/iCount_temp;
        end 
    end
    save([ResultDir,'/',SiteSet{iSite},'_ISCMap.mat'],'ISCMap');
end

for iSite = 1:length(SiteSet)
    ISFCMap = {};
    for iSession = 1:length(SessionNameSet) 
        for iSubject = 1:length(SubList)
            MeanISCorrelationMap = [];
            load([DataDir,'/',SessionNameSet{iSession},'/',SiteSet{iSite},'/',SubList{iSubject},'_ISCMap.mat']);
            %DMN core and DMN dmPFC
            iCount_temp=0;
            Sum_temp=0;
            for isub = 1:9
                for jsub = 10:18
                  Sum_temp=Sum_temp+MeanISCorrelationMap(isub,jsub);
                  Sum_temp=Sum_temp+MeanISCorrelationMap(jsub,isub);
                  iCount_temp=iCount_temp+2;
                end
            end
            ISFCMap{iSession}(iSubject,1) = Sum_temp/iCount_temp;
            
            %DMN core and DMN MTL
            iCount_temp=0;
            Sum_temp=0;
            for isub = 1:9
                for jsub = 19:24
                  Sum_temp=Sum_temp+MeanISCorrelationMap(isub,jsub);
                  Sum_temp=Sum_temp+MeanISCorrelationMap(jsub,isub);
                  iCount_temp=iCount_temp+2;
                end
            end
            ISFCMap{iSession}(iSubject,2) = Sum_temp/iCount_temp;
           
            %DMN dmPFC and DMN MTL
            iCount_temp=0;
            Sum_temp=0;
            for isub = 10:18
                for jsub = 19:24
                  Sum_temp=Sum_temp+MeanISCorrelationMap(isub,jsub);
                  Sum_temp=Sum_temp+MeanISCorrelationMap(jsub,isub);
                  iCount_temp=iCount_temp+2;
                end
            end
            ISFCMap{iSession}(iSubject,3) = Sum_temp/iCount_temp;
        end
    end
    save([ResultDir,'/',SiteSet{iSite},'_ISFCMap.mat'],'ISFCMap');
end