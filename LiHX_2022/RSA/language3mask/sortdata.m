clear;clc
outpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/convergentTwoexperiment/RSAROIanalysis/languageMask/1paper86/statisticanalysis/';
path={'/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/WYSWYT/RSAanalysisWYSWYT/RSAnewprepocess/ROIRSAanalysis/languagemasks/1paper86'; ...
    '/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RASanalyisis/ROIRSAanalysis/languageMask/1paper86'};

for ie=1:size(path,1)
    network=dir([path{ie,:},'/E*']);
    
    for in=1:size(network,1)
        outpath_net=[outpath,network(in).name(3:11)];
        if ~exist(outpath_net,'dir'); mkdir(outpath_net);end
        
        copyfile([path{ie,:},'/',network(in).name],outpath_net)
    end
end