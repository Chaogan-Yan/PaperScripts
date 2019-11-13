%To investigate the correlation between scale scores (RRS, BDI) and neural
%signals

clear;clc;
load '/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Demographic_data/DemographicInfo_V2.mat';
SiteSet = {'IPCAS','PKUGE','PKUSIMENS'};
Workdir='/MD3860F/RfMRILab/ChenX/Rumination_project/Analysis/Analysis4Publish';
Resultdir = '/MD3860F/RfMRILab/ChenX/Rumination_project/Analysis/Analysis4Publish/ScaleAnalysis/';
Datadir='/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing';
Scaledatadir = '/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Raw/Behavior_data/Scale_Sorted';
SubList = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/Analysis/IPCAS_Sublist.txt');

%Z transform
% Brooding = (Brooding-mean(Brooding))/std(Brooding);
% Depression = (Depression-mean(Depression))/std(Depression);
% Reflection = (Reflection-mean(Reflection))/std(Reflection);
% Rumination = (Rumination-mean(Rumination))/std(Rumination);

%Calculate the correlation between network wise FCs and scale
for iSite = 1:length(SiteSet)
    load([Workdir,'/',SiteSet{iSite},'_Network_FC_across_brain_Smooth.mat']);
    load([Scaledatadir,'/',SiteSet{iSite},'_ScaleData.mat']);
    RumThinkingContent = cell2mat(RumThinkingContent);
    DisThinkingContent = cell2mat(DisThinkingContent);
    ResidualThinkingContent = RumThinkingContent - DisThinkingContent;
    Scaleset = {Rumination,Brooding,Reflection};
    % extract FC under each network, then corr this FC with scales      
    for i = 1:5
        for j = i+1:6 
           %betweennetwork correlation
           BetweenNetworkFC = {};
           ResidualBetweenNeworkFC = [];
           
           for iSub = 1:41
                for iCondition = 1:5
                    BetweenNetworkFC{iSub,iCondition} =  AllBetween_Network_FC{iCondition}(i,j,iSub);
                end
           end
           BetweenNetworkFC = cell2mat(BetweenNetworkFC);
           ResidualBetweenNeworkFC = BetweenNetworkFC(:,4)-BetweenNetworkFC(:,5);
           
           for iScale = 1:3
               for iCondition = 1:5
                   [ScaleBetweenFC_RMap{iCondition,iScale}(i,j),ScaleBetweenFC_PMap{iCondition,iScale}(i,j)] = partialcorr(BetweenNetworkFC(:,iCondition),Scaleset{iScale},Depression);
               end
                   [ScaleBetweenFC_RMap{6,iScale}(i,j),ScaleBetweenFC_PMap{6,iScale}(i,j)] = partialcorr(ResidualBetweenNeworkFC,Scaleset{iScale},Depression);
           end
           
           for iTC = 1:10
               for iCondition = 1:5
                [MRIScaleRum_BetweenFC_RMap{i,j}(iCondition,iTC),MRIScaleRum_BetweenFC_PMap{i,j}(iCondition,iTC)] = corr(RumThinkingContent(:,iTC), BetweenNetworkFC(:,iCondition));
                [MRIScaleResidual_BetweenFC_RMap{i,j}(iCondition,iTC),MRIScaleResidual_BetweenFC_PMap{i,j}(iCondition,iTC)] = corr(ResidualThinkingContent(:,iTC), BetweenNetworkFC(:,iCondition));
               end
                [MRIScaleRum_BetweenFC_RMap{i,j}(6,iTC),MRIScaleRum_BetweenFC_PMap{i,j}(6,iTC)] = corr(RumThinkingContent(:,iTC), ResidualBetweenNeworkFC);
                [MRIScaleResidual_BetweenFC_RMap{i,j}(6,iTC),MRIScaleResidual_BetweenFC_PMap{i,j}(6,iTC)] = corr(ResidualThinkingContent(:,iTC), ResidualBetweenNeworkFC);
           end
        end
    end
    
    
       %within network correlation
       for iCondition = 1:5
           for iNetwork = 1:6
               for iScale = 1:3
                 [ScaleWithinFC_RMap{iCondition,iScale}(iNetwork),ScaleWithinFC_PMap{iCondition,iScale}(iNetwork)] = partialcorr(AllWithin_Network_FC{iCondition}(:,iNetwork),Scaleset{iScale},Depression);
               end
           end
       end

       for iNetwork = 1:6
            ResidualWithinNetwork_FC(:,iNetwork) = AllWithin_Network_FC{4}(:,iNetwork)-AllWithin_Network_FC{5}(:,iNetwork);
            for iScale = 1:3
                [ScaleWithinFC_RMap{6,iScale}(iNetwork),ScaleWithinFC_PMap{6,iScale}(iNetwork)] = partialcorr(ResidualWithinNetwork_FC(:,iNetwork),Scaleset{iScale},Depression);
            end
       end

       for iNetwork = 1:6
           for iTC = 1:10
               for iCondition = 1:5
                [MRIScaleRum_WithinFC_RMap{iNetwork}(iCondition,iTC),MRIScaleRum_WithinFC_PMap{iNetwork}(iCondition,iTC)] = corr(RumThinkingContent(:,iTC), AllWithin_Network_FC{iCondition}(:,iNetwork));
                [MRIScaleResidual_WithinFC_RMap{iNetwork}(iCondition,iTC),MRIScaleResidual_WithinFC_PMap{iNetwork}(iCondition,iTC)] = corr(ResidualThinkingContent(:,iTC), AllWithin_Network_FC{iCondition}(:,iNetwork));
               end
                [MRIScaleRum_WithinFC_RMap{iNetwork}(6,iTC),MRIScaleRum_WithinFC_PMap{iNetwork}(6,iTC)] = corr(RumThinkingContent(:,iTC), ResidualWithinNetwork_FC(:,iNetwork));
                [MRIScaleResidual_WithinFC_RMap{iNetwork}(6,iTC),MRIScaleResidual_WithinFC_PMap{iNetwork}(6,iTC)] = corr(ResidualThinkingContent(:,iTC), ResidualWithinNetwork_FC(:,iNetwork));
           end
       end
    
    save([Resultdir,'/',SiteSet{iSite},'_NetworkFCScaleCorrelation.mat'],'*Map*');
end


%Calculate the correlation between ROI level FCs and scale scores
% for iSite = 1:length(SiteSet)
%     RestFCDir=['/',SiteSet{iSite},'_rest/Results/ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
%     HappyFCDir=['/',SiteSet{iSite},'_task/Results/ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
%     SadFCDir=['/',SiteSet{iSite},'_task/S2_Results/S2_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
%     RumFCDir=['/',SiteSet{iSite},'_task/S3_Results/S3_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
%     DisFCDir=['/',SiteSet{iSite},'_task/S4_Results/S4_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_'];
%     FCDirset = {RestFCDir, HappyFCDir,SadFCDir,RumFCDir,DisFCDir};
%     Scaleset = {Rumination,Brooding,Reflection};
%     AllFCMatrix = zeros(41,5);
%     ROIRMap = {};
%     ROIPMap = {};
%     
%     for iROI = 1:23
%        for jROI = iROI+1:24
%             for iFC = 1:length(FCDirset)
%                 for i = 1:length(SubList)
%                    load([Datadir,FCDirset{iFC},SubList{i},'.mat']);
%                    AllFCMatrix(i,iFC) = ROICorrelation_FisherZ(iROI,jROI);
%                 end 
%                 for iScale = 1:3
%                     % a 6x3 matrix 6: 5 conditions and rum-dis FC
%                     % differences; 3: rum, brooding, reflection
%                     ResidualFC = [];
%                     [ROIRMap{iFC,iScale}(iROI,jROI),ROIPMap{iFC,iScale}(iROI,jROI)] = partialcorr(AllFCMatrix(:,iFC),Scaleset{iScale},Depression);
%                     ResidualFC = AllFCMatrix(:,4) - AllFCMatrix(:,5);
%                     [ROIRMap{6,iScale}(iROI,jROI),ROIPMap{6,iScale}(iROI,jROI)] = partialcorr(ResidualFC,Scaleset{iScale},Depression);
%                 end
%             end
%             
%        end
%     end
%     
%     %fdr correction
%     for iFC = 1:6
%         for iScale = 1:3
%             
%             icount = 0;
%             ROIPMap4FDR = [];
%             ROIPMapAfterFDR = [];
%             for i = 1:23
%                 for j = i+1:24
%                     icount = icount +1;
%                     ROIPMap4FDR(icount,1) = ROIPMap{iFC,iScale}(i,j);
%                 end
%             end
%             ROIPMapAfterFDR = mafdr(ROIPMap4FDR,'BHFDR',true);
%             icount = 0;
%             for i = 1:23
%                 for j = i+1:24
%                     icount = icount+1;
%                     ROIFDR_PMap{iFC,iScale}(i,j) = ROIPMapAfterFDR(icount);
%                 end
%             end
%             
%         end
%     end
%     
%     save([Workdir,'/',SiteSet{iSite},'ROIScaleCorrelation.mat'],'ROIRMap','ROIPMap','ROIFDR_PMap');
% end
