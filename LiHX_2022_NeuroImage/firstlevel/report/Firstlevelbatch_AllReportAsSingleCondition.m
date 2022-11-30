clc;clear;
datadir='/mnt/Data3/RfMRILab/Lihuixian/Processing/dpabiTaskprocess/2020ExperimentReportS2ControlReportDelete20s/FunImgARWS/';  

sublistpath = '/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/Report/firstlevel/';
SubID=dir([sublistpath, '/sub*']);

%%%spm process 
spm_get_defaults          
cwd='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/Report/firstlevel/';
spm_jobman('initcfg');

for isub = 1:size(SubID,1)
    jobs{1}.stats{1}.fmri_spec.dir=cellstr([cwd SubID(isub).name]);   
    jobs{1}.stats{1}.fmri_spec.timing.units='secs';   
    jobs{1}.stats{1}.fmri_spec.timing.RT=2;   
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t = 37;    
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t0 = 19;
    
    muldirFiles=[ datadir,SubID(isub).name]; 
    imagefilepath=spm_select('FPList',muldirFiles,'.*\.nii$');
    
    [Datanii,Header]=y_Read(imagefilepath);
    mulimagefile={};
    for idata = 1:size(Datanii,4)
        mulimagefile{idata,1} = [imagefilepath, ',', num2str(idata)];
    end
    jobs{1}.stats{1}.fmri_spec.sess.scans=mulimagefile;
 
    %set condition
    mulcondition_file_name=[SubID(isub).name,'.mat'];
    mulcondpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/Report/Behaviordata/TimeInfodelete20s';
    mulconditions=load(fullfile(mulcondpath,mulcondition_file_name));
    jobs{1}.stats{1}.fmri_spec.sess.cond(1).name='SpeakThoughtOnset'; 
    jobs{1}.stats{1}.fmri_spec.sess.cond(1).onset=mulconditions.SpeakBegin_onsetTime;
    jobs{1}.stats{1}.fmri_spec.sess.cond(1).duration=mulconditions.Duration;
    
    headmotionpath='/mnt/Data3/RfMRILab/Lihuixian/Processing/dpabiTaskprocess/2020ExperimentReportS2ControlReportDelete20s/RealignParameter/';
    headmotionfile=[headmotionpath,SubID(isub).name];
    mulmotionfiles=spm_select('FPList',headmotionfile,'^rp.*\.txt$');
    jobs{1}.stats{1}.fmri_spec.sess.multi_reg=cellstr(mulmotionfiles);
    %set brainmask
    jobs{1}.stats{1}.fmri_spec.mthresh = 0;
    jobs{1}.stats{1}.fmri_spec.mask={'/mnt/Data3/RfMRILab/Lihuixian/DPABI_V5.1_201230/Templates/BrainMask_05_61x73x61.img,1'};
     
    resultpath=[cwd SubID(isub).name];
    jobs{1}.stats{2}.fmri_est.spmmat=cellstr(fullfile(resultpath,'SPM.mat')); 
    save(fullfile(resultpath,'modelspecification.mat'),'jobs'); 
    
    spm_jobman('run',jobs);
    
    clear jobs;
    
    %set contrast
    resultpath=[cwd SubID(isub).name];
    jobs{1}.stats{1}.con.spmmat=cellstr(fullfile(resultpath,'SPM.mat'));
    jobs{1}.stats{1}.con.consess{1}.tcon=struct('name','SpeakThought','convec',[ 1 0 0 0 0 0 0],'sessrep','none');   

    
     spm_jobman('run',jobs);
     clear jobs;
end

    

    
    
    
    
    
    
    
    