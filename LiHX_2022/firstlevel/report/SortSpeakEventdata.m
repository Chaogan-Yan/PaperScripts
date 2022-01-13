clc;clear
path = '/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RASanalyisis/firstLevelSpeakEvent/firstlevel';
outpath = '/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RASanalyisis/firstLevelSpeakEvent/SpeakEventBeta/';
SubID=dir([path,'/sub*']);
SpeakEventSub45=load('/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RASanalyisis/firstLevelSpeakEvent/SpeakEvents/SpeakEvent.mat');

for isub =1:size(SubID,1)
    subpath=fullfile(path,SubID(isub).name);
    outpathSub=fullfile(outpath,SubID(isub).name);
    mkdir(outpathSub)
    
    subid=str2double(SubID(isub).name(isstrprop(SubID(isub).name,'digit')));
    events=SpeakEventSub45.SpeakEvent(find(SpeakEventSub45.SpeakEvent(:,1)==subid),2);
   
    for ie = 1:events
        if ie<10
            evenname=['beta_000',num2str(ie),'.nii'];
        else
            evenname=['beta_00',num2str(ie),'.nii'];
        end
        findcon = spm_select('FPList',subpath,evenname);
        copyfile( findcon,outpathSub);
        newname = [SubID(isub).name,'Even',num2str(ie),'_',evenname];
        movefile(fullfile(outpathSub,evenname), fullfile(outpathSub,newname));
    end
    
end
