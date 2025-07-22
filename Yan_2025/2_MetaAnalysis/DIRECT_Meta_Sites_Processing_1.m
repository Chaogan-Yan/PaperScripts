

addpath /mnt/Data6/RfMRILab/Yan/YAN_Scripts/MDD/REST-meta-MDD_Surf/Analysis/ENIGMACo/MetaStats/AfterQCExclusion/


DataDir='/mnt/Data7/RfMRILab/Yan/YAN_Work/REST_meta-MDD_Surf/Organized/Processing/ResultsS10/';

OutputDir='/mnt/Data7/RfMRILab/Yan/YAN_Work/REST_meta-MDD_Surf/Analysis/MetaAnalysis/AfterQCExclusion/MetaMaps/';

load /mnt/Data7/RfMRILab/Yan/YAN_Work/REST_meta-MDD_Surf/Analysis/SubInfo/ETable_AfterQCExclusion.mat


addpath(genpath('/mnt/Data6/RfMRILab/Yan/YAN_Scripts/MDD/REST-meta-MDD_Surf/Analysis/ENIGMACo/DPABI_V6_ForENIGMA_RMD/'))
addpath /mnt/Data6/RfMRILab/Yan/YAN_Program/spm12/
addpath /mnt/Data6/RfMRILab/Yan/YAN_Scripts/MDD/REST-meta-MDD_Surf/Analysis/ENIGMACo/Scripts_ENIGMA_RMD/


DIRECT_Meta_Sites(DataDir,OutputDir,ETable2);


