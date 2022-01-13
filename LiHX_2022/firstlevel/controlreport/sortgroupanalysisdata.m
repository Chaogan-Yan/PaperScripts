clc;clear
path = '/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/ControlReport/firstlevel/';
outpath1 = '/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/ControlReport/groupanalysis/data/';
mkdir(outpath1)
SubID=dir([path, '/sub*']);

for isub =1:size(SubID,1)
     subpath = [path,SubID(isub).name];
     
     con1 = spm_select('FPList',subpath,'con_0001.nii');
     copyfile( con1,outpath1);
     newname1 = [SubID(isub).name,'_con_0001.nii'];
     movefile(fullfile(outpath1,'con_0001.nii'), fullfile(outpath1,newname1));
     
end
     