clear; clc;

poolsize = 30;
p = gcp('nocreate');
if isempty(p)
    parpool('local', poolsize);
elseif p.NumWorkers ~= poolsize
    delete(p);
    parpool('local', poolsize);
end

% ------------------ YOUR DATA ------------------
dataDir = '/mnt/Data8/RfMRILab/gaoqinglin/2_TMS/DIRECT/AMN_Surf/AMN_DMN/ROISignals/AMN_DMN/ROISignal_TP230/NC/'; 
working_dir ='/mnt/Data8/RfMRILab/gaoqinglin/2_TMS/DIRECT/AMN_Surf/AMN_DMN/CAP/NC/CAP_Results';                
if ~exist(working_dir, 'dir'); mkdir(working_dir); end

% ------------------ K-means params ------------------
K_Means_Para.Range = 2:1:10;                 
K_Means_Para.Num = length(K_Means_Para.Range);
K_Means_Para.DistanceFunc = 'correlation';   
K_Means_Para.ReplicateNum = 50;             

% ------------------ RUN ------------------
CAP_Framework_ROI_TS(working_dir, dataDir, K_Means_Para);
