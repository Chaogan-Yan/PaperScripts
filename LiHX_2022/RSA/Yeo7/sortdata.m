clear;clc
outpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/convergentTwoexperiment/RSAROIanalysis/1paper86Yeo7/';
path={'/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/WYSWYT/RSAanalysisWYSWYT/RSAnewprepocess/ROIRSAanalysis/1paperYeo7'; ...
    '/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RASanalyisis/ROIRSAanalysis/1paperYeo7network'};

for ie=1:size(path,1)
    network=dir([path{ie,:},'/E*']);
    
    for in=1:size(network,1)
        outpath_net=[outpath,'Yeo',num2str(in)];
        if ~exist(outpath_net,'dir'); mkdir(outpath_net);end
        
        copyfile([path{ie,:},'/',network(in).name],outpath_net)
    end
end