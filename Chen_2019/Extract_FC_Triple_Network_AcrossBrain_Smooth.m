%I defined a total of 66 ROIs including triple networks (FPCN, DMN and SN)
%as well as sgACC and 6 subcortical regions
%This script is to calculate the FCs of each networks
%modified from Extract_FC_Dixon
%20181120 ChenXiao

%This script is trying to calculate FCs across brain

clear;clc;

%initialization
Workdir='/mnt/Data/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing';
Resultdir='/MD3860F/RfMRILab/ChenX/Rumination_project/Analysis/Analysis4Publish';
SubList = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/Analysis/IPCAS_Sublist.txt');
SiteSet = {'IPCAS','PKUGE','PKUSIMENS'};

for iSite = 1:length(SiteSet)
    RestFCDir = ['/',SiteSet{iSite},'_rest/Results/ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
    HappyFCDir = ['/',SiteSet{iSite},'_task/Results/ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
    SadFCDir = ['/',SiteSet{iSite},'_task/S2_Results/S2_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
    RumFCDir = ['/',SiteSet{iSite},'_task/S3_Results/S3_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
    DisFCDir = ['/',SiteSet{iSite},'_task/S4_Results/S4_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];

    FCDirset = {RestFCDir, HappyFCDir,SadFCDir,RumFCDir,DisFCDir,};
    %For the serial number for these 66 ROIs, please refer to the table "ROI analysis" in OneNote.

    %1,Rest; 2, Happy; 3, Sad; 4, Rum; 5, Dis;
    AllWithin_Network_FC = {};
    %1,Rest; 2, Happy; 3, Sad; 4, Rum; 5, Dis;
    AllBetween_Network_FC = {};

for iFC = 1:length(FCDirset)
        %set within-network FCs: a n x 6 matrix
        %1, DN_core; 2, DN_dmPFC; 3, DN_MTL; 4, FPCN; 5, SN; 6, DMN_total
        %n, subjects' number
        Within_Network_FC = zeros(length(SubList),6);
        
        %set Between-network FCs: a n x 6 x 6  matrix
        %1, DN_core; 2, DN_dmPFC; 3, DN_MTL; 4, FPCN; 5, SN; 6, DMN_total
        %n, subjects' number
        Between_Network_FC = zeros(6,6,length(SubList));
    for i = 1:length(SubList)
       load([Workdir,FCDirset{iFC},SubList{i},'.mat']);
       
       %%%%%%%%%%%%%%%%%%%%%%% WITHIN_NETWORKS %%%%%%%%%%%%%%%%%%%%%%%%%%%
       %calculate the FC within DN_core
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:8
          for jsub =  isub+1:9
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC(i,1)= Sum_temp/iCount_temp;
       

       %calculate the FC within DN_dmPFC
       iCount_temp=0;
       Sum_temp=0;
       for isub = 10:17
          for jsub =  isub+1:18
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC(i,2)= Sum_temp/iCount_temp;
       
       %calculate the FC within DN_MTL
       iCount_temp=0;
       Sum_temp=0;
       for isub = 19:23
          for jsub =  isub+1:24
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC(i,3)= Sum_temp/iCount_temp;
       
       
        %calculate the FC within FPCN
       iCount_temp=0;
       Sum_temp=0;
       for isub = 25:34
          for jsub =  isub+1:35
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC(i,4)= Sum_temp/iCount_temp;
      
       
        %calculate the FC within SN
       iCount_temp=0;
       Sum_temp=0;
       for isub = 36:58
          for jsub =  isub+1:59
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC(i,5)= Sum_temp/iCount_temp;
     
       
       %the FC within DMN-total: the mean of three within_subsystems' FC
       Within_Network_FC(i,6) = (Within_Network_FC(i,1) + Within_Network_FC(i,2) + Within_Network_FC(i,3))/3;
       
       %%%%%%%%%%%%%%%%%%%%%%% BETWEEN_NETWORKS %%%%%%%%%%%%%%%%%%%%%%%%%%%
      
       %%%%%%%%%%%%%%%%%%%%%%% DN_core and others %%%%%%%%%%%%%%%%%%%%%%%%
       
       %calculate FCs between DN_core and DN_dmPFC
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:9
          for jsub = 10:18
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(1,2,i) = Sum_temp/iCount_temp;
      
       
       %calculate FCs between DN_core and DN_MTL
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:9
          for jsub = 19:24
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(1,3,i) = Sum_temp/iCount_temp;
       
       
       %calculate FCs between DN_core and FPCN
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:9
          for jsub = 25:35
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(1,4,i) = Sum_temp/iCount_temp;
       
       
       %calculate FCs between DN_core and SN
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:9
          for jsub = 36:59
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(1,5,i) = Sum_temp/iCount_temp;
      
       
       %For those between DMN_subsystems and DMN_total, I calculated them
       %and the rest ROIs
       %calculate FCs between DN_core and DMN_total
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:9
          for jsub = 10:24
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(1,6,i) = Sum_temp/iCount_temp;
       
       %%%%%%%%%%%%%%%%%%%% DN_dmPFC and others %%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate FCs between DN_dmPFC and DN_MTL
       iCount_temp=0;
       Sum_temp=0;
       for isub = 10:18
          for jsub = 19:24
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(2,3,i) = Sum_temp/iCount_temp;
      
       
        %calculate FCs between DN_dmPFC and FPCN
       iCount_temp=0;
       Sum_temp=0;
       for isub = 10:18
          for jsub = 25:35
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(2,4,i) = Sum_temp/iCount_temp;
       
        %calculate FCs between DN_dmPFC and SN
       iCount_temp=0;
       Sum_temp=0;
       for isub = 10:18
          for jsub = 36:59
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(2,5,i) = Sum_temp/iCount_temp;
       
       %calculate FCs between DN_dmPFC and DMN_total
        iCount_temp=0;
       Sum_temp=0;
       for isub = 10:18
          for jsub = [1:9, 19:24]
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(2,6,i) = Sum_temp/iCount_temp;
       
       %%%%%%%%%%%%%%%%%%%%%%DN_MTL and others%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate FCs between DN_MTL and FPCN
       iCount_temp=0;
       Sum_temp=0;
       for isub = 19:24
          for jsub = 25:35
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(3,4,i) = Sum_temp/iCount_temp;
       
        %calculate FCs between DN_MTL and SN
       iCount_temp=0;
       Sum_temp=0;
       for isub = 19:24
          for jsub = 36:59
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(3,5,i) = Sum_temp/iCount_temp;
       
       %calculate FCs between DN_MTL and DMN_total
       iCount_temp=0;
       Sum_temp=0;
       for isub = 19:24
          for jsub = 1:18
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(3,6,i) = Sum_temp/iCount_temp;
       
       %%%%%%%%%%%%%%%%%%%%%%% FPCN and others %%%%%%%%%%%%%%%%%%%%%%%%%%%
       %calculate FCs between FPCN and SN
       iCount_temp=0;
       Sum_temp=0;
       for isub = 25:35
          for jsub = 36:59
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(4,5,i) = Sum_temp/iCount_temp;
       
       %calculate FCs between FPCN and DMN_total
       iCount_temp=0;
       Sum_temp=0;
       for isub = 25:35
          for jsub = 1:24
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(4,6,i) = Sum_temp/iCount_temp;
      
       
       %calculate FCs between SN and DMN_total
       iCount_temp=0;
       Sum_temp=0;
       for isub = 36:59
          for jsub = 1:24
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC(5,6,i) = Sum_temp/iCount_temp;
      
    end
    %Calculate the whole brain FCs
    AllWithin_Network_FC{iFC} = Within_Network_FC;
    AllBetween_Network_FC{iFC} = Between_Network_FC;
end
    save([Resultdir,'/',SiteSet{iSite},'_Network_FC_across_brain_Smooth.mat'],'All*');
end