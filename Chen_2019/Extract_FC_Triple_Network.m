%I defined a total of 66 ROIs including triple networks (FPCN, DMN and SN)
%as well as sgACC and 6 subcortical regions
%This script is to calculate the FCs of each networks
%modified from Extract_FC_Dixon
%20181120 ChenXiao

clear all;clc;

%initialization
Workdir='/mnt/Data/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing';
Resultdir='/MD3860F/RfMRILab/ChenX/Rumination_project/Analysis/Analysis4Publish/HemisphereApart';
SubList = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/Analysis/IPCAS_Sublist.txt');
SiteSet = {'IPCAS','PKUGE','PKUSIMENS'};

for iSite = 1:length(SiteSet)
    
RestFCDir = ['/',SiteSet{iSite},'_rest/Results/ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
HappyFCDir = ['/',SiteSet{iSite},'_task/Results/ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
SadFCDir = ['/',SiteSet{iSite},'_task/S2_Results/S2_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
RumFCDir = ['/',SiteSet{iSite},'_task/S3_Results/S3_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
DisFCDir = ['/',SiteSet{iSite},'_task/S4_Results/S4_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];

FCDirset = {RestFCDir, HappyFCDir,SadFCDir,RumFCDir,DisFCDir,};

%Dixon said FC across hemisphere should be excluded because they are
%affected by some comfounding issues. So I calculated each networks' FCs within
%each hemisphere separately, then average them.

%For the serial number for these 43 ROIs, please refer to the table in OneNote.


%1,Rest; 2, Happy; 3, Sad; 4, Rum; 5, Dis;
AllWithin_Network_FC = {};
AllWithin_Network_FC_L = {};
AllWithin_Network_FC_L = {};

%1,Rest; 2, Happy; 3, Sad; 4, Rum; 5, Dis;
AllBetween_Network_FC = {};
AllBetween_Network_FC_L = {};
AllBetween_Network_FC_R = {};

for iFC = 1:length(FCDirset)
    %set within-network FCs: a n x 6 matrix
        %1, DN_core; 2, DN_dmPFC; 3, DN_MTL; 4, FPCN; 5, SN; 6, DMN_total
        %n, subjects' number
        Within_Network_FC_L = zeros(length(SubList),6);
        Within_Network_FC_R = zeros(length(SubList),6);
        Within_Network_FC = zeros(length(SubList),6);
        %set Between-network FCs: a n x 6 x 6  matrix
        %1, DN_core; 2, DN_dmPFC; 3, DN_MTL; 4, FPCN; 5, SN; 6, DMN_total
        %n, subjects' number
        Between_Network_FC_L = zeros(6,6,length(SubList));
        Between_Network_FC_R = zeros(6,6,length(SubList));
        Between_Network_FC = zeros(6,6,length(SubList));
    for i = 1:length(SubList)
        
       load([Workdir,FCDirset{iFC},SubList{i},'.mat']);
       
       %%%%%%%%%%%%%%%%%%%%%%% WITHIN_NETWORKS %%%%%%%%%%%%%%%%%%%%%%%%%%%
       %calculate the FC within DN_core_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:3
          for jsub =  isub+1:4
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC_L(i,1)= Sum_temp/iCount_temp;
       %calculate the FC within DN_core_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 5:8
          for jsub =  isub+1:9
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC_R(i,1)= Sum_temp/iCount_temp;

       %calculate the FC within DN_dmPFC_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 10:13
          for jsub =  isub+1:14
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC_L(i,2)= Sum_temp/iCount_temp;
       %calculate the FC within DN_dmPFC_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 15:17
          for jsub =  isub+1:18
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC_R(i,2)= Sum_temp/iCount_temp;
       
       %calculate the FC within DN_MTL_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 19:20
          for jsub =  isub+1:21
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC_L(i,3)= Sum_temp/iCount_temp;
       %calculate the FC within DN_MTL_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 22:23
          for jsub =  isub+1:24
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC_R(i,3)= Sum_temp/iCount_temp;
       
        %calculate the FC within FPCN_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 25:29
          for jsub =  isub+1:30
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC_L(i,4)= Sum_temp/iCount_temp;
       %calculate the FC within FPCN_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 31:34
          for jsub =  isub+1:35
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC_R(i,4) = Sum_temp/iCount_temp;
       
        %calculate the FC within SN_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 36:45
          for jsub =  isub+1:46
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC_L(i,5)= Sum_temp/iCount_temp;
       %calculate the FC within SN_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 47:58
          for jsub =  isub+1:59
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Within_Network_FC_R(i,5) = Sum_temp/iCount_temp;
       
       %the FC within DMN-total: the mean of three within_subsystems' FC
       Within_Network_FC_R(i,6) = (Within_Network_FC_R(i,1) + Within_Network_FC_R(i,2) + Within_Network_FC_R(i,3))/3;
       Within_Network_FC_L(i,6) = (Within_Network_FC_L(i,1) + Within_Network_FC_L(i,2) + Within_Network_FC_L(i,3))/3;
       
       %%%%%%%%%%%%%%%%%%%%%%% BETWEEN_NETWORKS %%%%%%%%%%%%%%%%%%%%%%%%%%%
      
       %%%%%%%%%%%%%%%%%%%%%%% DN_core and others %%%%%%%%%%%%%%%%%%%%%%%%
       
       %calculate FCs between DN_core_L and DN_dmPFC_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:4
          for jsub = 10:14
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(1,2,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_core_R and DN_dmPFC_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 5:9
          for jsub = 15:18
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(1,2,i) = Sum_temp/iCount_temp;
       
       %calculate FCs between DN_core_L and DN_MTL_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:4
          for jsub = 19:21
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(1,3,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_core_R and DN_MTL_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 5:9
          for jsub = 22:24
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(1,3,i) = Sum_temp/iCount_temp;
       
       
       %calculate FCs between DN_core_L and FPCN_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:4
          for jsub = 25:30
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(1,4,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_core_R and FPCN_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 5:9
          for jsub = 31:35
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(1,4,i) = Sum_temp/iCount_temp;
       
       %calculate FCs between DN_core_L and SN_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:4
          for jsub = 36:46
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(1,5,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_core_R and SN_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 5:9
          for jsub = 47:59
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(1,5,i) = Sum_temp/iCount_temp;
       
       
       %For those between DMN_subsystems and DMN_total, I calculated them
       %and the rest ROIs
       %calculate FCs between DN_core_L and DMN_total_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 1:4
          for jsub = [10:14, 19:21]
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(1,6,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_core_R and DMN_total_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 5:9
          for jsub = [15:18, 22:24]
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(1,6,i) = Sum_temp/iCount_temp;
       
       
       %%%%%%%%%%%%%%%%%%%% DN_dmPFC and others %%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate FCs between DN_dmPFC_L and DN_MTL_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 10:14
          for jsub = 19:21
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(2,3,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_dmPFC_R and DN_MTL_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 15:18
          for jsub = 22:24
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(2,3,i) = Sum_temp/iCount_temp;
       
        %calculate FCs between DN_dmPFC_L and FPCN_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 10:14
          for jsub = 25:30
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(2,4,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_dmPFC_R and FPCN_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 15:18
          for jsub = 31:35
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(2,4,i) = Sum_temp/iCount_temp;
       
        %calculate FCs between DN_dmPFC_L and SN_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 10:14
          for jsub = 36:46
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(2,5,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_dmPFC_R and SN_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 15:18
          for jsub = 47:59
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(2,5,i) = Sum_temp/iCount_temp;
       
       %calculate FCs between DN_dmPFC_L and DMN_total_L
        iCount_temp=0;
       Sum_temp=0;
       for isub = 10:14
          for jsub = [1:4, 19:21]
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(2,6,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_dmPFC_R and DMN_total_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 15:18
          for jsub =  [5:9, 22:24]
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(2,6,i) = Sum_temp/iCount_temp;
       
       %%%%%%%%%%%%%%%%%%%%%%DN_MTL and others%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate FCs between DN_MTL_L and FPCN_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 19:21
          for jsub = 25:30
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(3,4,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_MTL_R and FPCN_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 22:24
          for jsub = 31:35
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(3,4,i) = Sum_temp/iCount_temp;
       
        %calculate FCs between DN_MTL_L and SN_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 19:21
          for jsub = 36:46
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(3,5,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_MTL_R and FPCN_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 22:24
          for jsub = 47:59
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(3,5,i) = Sum_temp/iCount_temp;
       
       %calculate FCs between DN_MTL_L and DMN_total_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 19:21
          for jsub = [1:4, 10:13]
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(3,6,i) = Sum_temp/iCount_temp;
       %calculate FCs between DN_MTL_R and DMN_total_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 22:24
          for jsub =  [5:9, 15:18]
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(3,6,i) = Sum_temp/iCount_temp;
       
       %%%%%%%%%%%%%%%%%%%%%%% FPCN and others %%%%%%%%%%%%%%%%%%%%%%%%%%%
       %calculate FCs between FPCN_L and SN_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 25:30
          for jsub = 36:46
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(4,5,i) = Sum_temp/iCount_temp;
       %calculate FCs between FPCN_R and SN_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 31:35
          for jsub =  47:59
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(4,5,i) = Sum_temp/iCount_temp;
       
       %calculate FCs between FPCN_L and DMN_total_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 25:30
          for jsub = [1:4, 10:14, 19:21]
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(4,6,i) = Sum_temp/iCount_temp;
       %calculate FCs between FPCN_R and DMN_total_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 31:35
          for jsub =  [5:9, 15:18, 22:24]
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(4,6,i) = Sum_temp/iCount_temp;
       
       %calculate FCs between SN_L and DMN_total_L
       iCount_temp=0;
       Sum_temp=0;
       for isub = 36:46
          for jsub = [1:4, 10:14, 19:21]
              Sum_temp = Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp = iCount_temp+1;
          end
       end
       Between_Network_FC_L(5,6,i) = Sum_temp/iCount_temp;
       %calculate FCs between SN_R and DMN_total_R
       iCount_temp=0;
       Sum_temp=0;
       for isub = 47:59
          for jsub =  [5:9, 15:18, 22:24]
              Sum_temp=Sum_temp+ROICorrelation_FisherZ(isub,jsub);
              iCount_temp=iCount_temp+1;
          end
       end
       Between_Network_FC_R(5,6,i) = Sum_temp/iCount_temp;
       
       Between_Network_FC(:,:,i) = (Between_Network_FC_L(:,:,i) + Between_Network_FC_R(:,:,i))/2;
    end
    %Calculate the whole brain FCs
    AllWithin_Network_FC_L{iFC} = Within_Network_FC_L;
    AllWithin_Network_FC_R{iFC} = Within_Network_FC_R;
    Within_Network_FC = (Within_Network_FC_L + Within_Network_FC_R)/2;
    AllWithin_Network_FC{iFC} = Within_Network_FC;
    AllBetween_Network_FC_L{iFC} = Between_Network_FC_L;
    AllBetween_Network_FC_R{iFC} = Between_Network_FC_R;
    AllBetween_Network_FC{iFC} = Between_Network_FC;
end


save([Resultdir,'/',SiteSet{iSite},'_Network_FC.mat'],'All*');

end