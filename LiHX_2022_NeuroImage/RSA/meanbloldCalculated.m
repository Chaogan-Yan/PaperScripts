clc;clear;
datadir='/mnt/Data3/RfMRILab/Lihuixian/Processing/dpabiTaskprocess/2020ExperimentReportS2ControlReportDelete20s/FunImgARWS/';  
outpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/MeanBlodAnalysis/MeanBlodSpeakEvent';
sublistpath = '/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RASanalyisis/firstLevelSpeakEvent/firstlevel/';
SubID=dir([sublistpath, '/sub*']);

Maskpath='/mnt/Data3/RfMRILab/Lihuixian/DPABI_V6.0_ForCamp/Templates/BrainMask_05_61x73x61.img';
[maskdata,maskheader]=y_Read(Maskpath);

for isub = 1:size(SubID,1)
    outpathSub=fullfile(outpath,SubID(isub).name);
    mkdir(outpathSub)
    subpath=fullfile(datadir,SubID(isub).name);
    imagefilepath=spm_select('FPList',subpath,'.*\.nii$'); 
    [Datanii,Header]=y_Read(imagefilepath);
    Datanii(find(isnan(Datanii))) = 0;
    DataSize=size(Datanii);
    %speack event
    mulcondition_file_name=[SubID(isub).name,'.mat'];
    mulcondpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/Report/Behaviordata/TimeInfodelete20s';
    mulconditions=load(fullfile(mulcondpath,mulcondition_file_name));
    
    onsetT=mulconditions.SpeakBegin_onsetTime/2; %change to time point TR=2
    endT=mulconditions.SpeakEnd_onsetTime/2;
    
    OnsetTime=floor(onsetT); 
    if OnsetTime(1,1)<1
        OnsetTime(1,1)=1;
    end
    
    EndTime=ceil(endT);
    if EndTime(end,1)>DataSize(1,4)
        EndTime(end,1)=DataSize(1,4);
    end
    
    for i =1:size(OnsetTime,1)
        evendata=mean(Datanii(:,:,:,(OnsetTime(i,1):EndTime(i,1))),4).*maskdata;
        if i<10
            evenname=fullfile(outpathSub,['MeanBlod_Event00',num2str(i),'.nii']);
        else
            evenname=fullfile(outpathSub,['MeanBlod_Event0',num2str(i),'.nii']);
        end
        
        y_Write(evendata,Header,evenname)
    end
end
    
    
    
    
    