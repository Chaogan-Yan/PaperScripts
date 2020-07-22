%Perform statistical analysis
clear;clc;

%original
addpath('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/functions');
ResultDir='/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/Analysis4Pub';
restsitesuffix = {'_rest','_rest','_rest_newVersion'};
tasksitesuffix = {'_task4Pub','_task4Pub','_task4Pub_oldVersion'};
if ~exist([ResultDir,'/rawFCmats'])
    mkdir([ResultDir,'/rawFCmats']);
end
if ~exist([ResultDir,'/csvfiles'])
    mkdir([ResultDir,'/csvfiles']);
end

% %scrubbing
% ResultDir='/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/stats_scrubbing';
% if ~exist(ResultDir); mkdir(ResultDir); end
% restsitesuffix = '_rest_scrubbing';
% tasksitesuffix = '_task_scrubbing';

%where you store the preprocessed fMRI data
DataDir = '/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing';

SiteSet = {'IPCAS','PKUGE','PKUSIEMENS'};
ThinkingContentNameSet = {'Past','Future','Self','Others','Positive','Negative','Image','Speech','Happy','Sad'};
%please replace this path with the 'IPCAS_Sublist.txt' you download
SubList = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/Analysis/IPCAS_Sublist.txt');
[~,ROIName] = xlsread('/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/DMN_ROI_Name.xlsx', ...
    1,'A1:A24');
%read nodes' location and color flag for plotting in BrainNetviewer
%please replace this path with the 'DMN_nodeV2.xlsx' you download
[~,NodeName] = xlsread('/MD3860F/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/nodeplot/DMN_nodeV2.xlsx', ...
    1,'A1:A24');
NodeCoord_flag = xlsread('/MD3860F/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/nodeplot/DMN_nodeV2.xlsx', ...
    1,'B1:E24');

%please replace this path with the 'DemographicInfo4Pub.mat' you download
load '/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/Analysis4Pub/DemographicInfo4Pub.mat';
ScaleMat = [Rumination,Brooding,Reflection];
% load thinking content ratings
ScaledataDir = '/MD3860F/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/Analysis4Pub';
rum_thinking_content = {};
dis_thinking_content = {};
rest_thinking_content = {};
for iSite = 1:length(SiteSet)
    load([ScaledataDir,'/',SiteSet{iSite},'_QuestionnaireData.mat']);
    rum_thinking_content{iSite} = cell2mat(RumThinkingContent);
    dis_thinking_content{iSite} = cell2mat(DisThinkingContent);
    rest_thinking_content{iSite} = cell2mat(RestThinkingContent);
end 

% %exclude left handed sub010 or sub011
% wantedSubMatrix = ones(length(SubList),1);
% wantedSubMatrix(11,1) = 0;
% WantedSubIndex = find(wantedSubMatrix);
% SubList = SubList(WantedSubIndex);
% ScaleMat = ScaleMat(WantedSubIndex,:);
% for iSite = 1:length(SiteSet)
%     rum_thinking_content{iSite} = rum_thinking_content{iSite}(WantedSubIndex,:);
%     dis_thinking_content{iSite} = dis_thinking_content{iSite}(WantedSubIndex,:);
%     rest_thinking_content{iSite} = rest_thinking_content{iSite}(WantedSubIndex,:);
% end
% ResultDir='/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/Analysis4Pub_noSub011';
% if ~exist(ResultDir) 
%     mkdir(ResultDir); 
% end
% if ~exist([ResultDir,'/csvfiles']) 
%     mkdir([ResultDir,'/csvfiles']); 
% end

%Extract head motion
%1,Rest; 2, Sad; 3, Rum; 4, Dis;
headmotion = {};
for iSite = 1:length(SiteSet)
    fprintf(['\nreading ',SiteSet{iSite},' head motions...\n\n']);
    for i=1:length(SubList)
        temp=load([DataDir,'/',SiteSet{iSite},restsitesuffix{iSite},'/RealignParameter/',SubList{i},'/FD_Jenkinson_',SubList{i},'.txt']);
        headmotion{iSite}(i,1)=mean(temp);
        temp=load([DataDir,'/',SiteSet{iSite},tasksitesuffix{iSite},'/RealignParameter/',SubList{i},'/FD_Jenkinson_',SubList{i},'.txt']);
        headmotion{iSite}(i,2)=mean(temp);
        temp=load([DataDir,'/',SiteSet{iSite},tasksitesuffix{iSite},'/RealignParameter/',SubList{i},'/S2_FD_Jenkinson_',SubList{i},'.txt']);
        headmotion{iSite}(i,3)=mean(temp);
        temp=load([DataDir,'/',SiteSet{iSite},tasksitesuffix{iSite},'/RealignParameter/',SubList{i},'/S3_FD_Jenkinson_',SubList{i},'.txt']);
        headmotion{iSite}(i,4)=mean(temp);
    end
end

%paired t test's regressors
nSub = length(SubList);
Regressors = [ones(nSub,1);-1*ones(nSub,1)];
for i=1:nSub
    SubjectRegressors(:,i) = zeros(nSub*2,1);
    SubjectRegressors(i:nSub:nSub*2,i) = 1;
end
OtherCovariatesMatrix = [];

for iSite = 1:length(SiteSet)
    fprintf(['\nperforming paired t tests ',SiteSet{iSite},'\n\n']);
    %%%Network level analysis
    %extract networkFC
    FCDataDir = [];
    FCDataDir = [DataDir,'/',SiteSet{iSite},restsitesuffix{iSite},'/Results/ROISignals_FunImgARCWFS'];
    RawFCMatsName = [ResultDir,'/rawFCmats/',SiteSet{iSite},'_restFC.mat'];
    [rest_within_networkFC,rest_between_networkFC,within_networkFCName,between_networkFCName] ...
        = c_ExtractDMN_FC(SubList,FCDataDir,RawFCMatsName);
    FCDataDir = [DataDir,'/',SiteSet{iSite},tasksitesuffix{iSite},'/S2_Results/S2_ROISignals_FunImgARCWFS'];
    RawFCMatsName = [ResultDir,'/rawFCmats/',SiteSet{iSite},'_rumFC.mat'];
    [rum_within_networkFC,rum_between_networkFC,~,~] ...
        = c_ExtractDMN_FC(SubList,FCDataDir,RawFCMatsName);
    FCDataDir = [DataDir,'/',SiteSet{iSite},tasksitesuffix{iSite},'/S3_Results/S3_ROISignals_FunImgARCWFS'];
    RawFCMatsName = [ResultDir,'/rawFCmats/',SiteSet{iSite},'_disFC.mat'];
    [dis_within_networkFC,dis_between_networkFC,~,~] ...
        = c_ExtractDMN_FC(SubList,FCDataDir,RawFCMatsName);
    FCDataDir = [DataDir,'/',SiteSet{iSite},tasksitesuffix{iSite},'/Results/ROISignals_FunImgARCWFS'];
    RawFCMatsName = [ResultDir,'/rawFCmats/',SiteSet{iSite},'_sadFC.mat'];
    [sad_within_networkFC,sad_between_networkFC,~,~] ...
        = c_ExtractDMN_FC(SubList,FCDataDir,RawFCMatsName);
    
    rum_dis_within_stats = {'network','t','p','Bon_p','fdr_p'};
    for ii =2:length(within_networkFCName)+1
        rum_dis_within_stats{ii,1} = within_networkFCName{ii-1};
    end
    rum_rest_within_stats = {'network','t','p','Bon_p','fdr_p'};
    for ii =2:length(within_networkFCName)+1
        rum_rest_within_stats{ii,1} = within_networkFCName{ii-1};
    end
    rum_sad_within_stats = {'network','t','p','Bon_p','fdr_p'};
    for ii =2:length(within_networkFCName)+1
        rum_sad_within_stats{ii,1} = within_networkFCName{ii-1};
    end
    sad_dis_within_stats = {'network','t','p','Bon_p','fdr_p'};
    for ii =2:length(within_networkFCName)+1
        sad_dis_within_stats{ii,1} = within_networkFCName{ii-1};
    end
    intersubVarWithin_rum_dis_stats = {'','Increase','decrease','equal'};
    for ii = 2:length(within_networkFCName)+1
        intersubVarWithin_rum_dis_stats{ii,1} = within_networkFCName{ii-1};
    end
    
    % comparing withinnetwork FC
    for iNetwork = 1:size(rum_within_networkFC,2)
        % rum vs. dis
        y = [rum_within_networkFC(:,iNetwork);dis_within_networkFC(:,iNetwork)];
        OtherCovariatesMatrix = [headmotion{iSite}(:,3);headmotion{iSite}(:,4)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        rum_dis_within_stats{iNetwork+1,2} = TF_ForContrast;
        rum_dis_within_stats{iNetwork+1,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        rum_dis_within_stats{iNetwork+1,4} = rum_dis_within_stats{iNetwork+1,3}*3;
        rum_dis_within_p4fdr(iNetwork,1) = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        
        %As the reviewer required, I count the number of increase or
        %decreased
        increaseCount = 0;
        decreaseCount = 0;
        equalCount = 0;
        for ii = 1:nSub
            if y(ii) > y(ii+nSub)
                increaseCount = increaseCount+1;
            elseif y(ii) == y(ii+nSub)
                equalCount = equalCount+1;
            else
                decreaseCount = decreaseCount+1;
            end
        end
        intersubVarWithin_rum_dis_stats{iNetwork+1,2} = increaseCount;
        intersubVarWithin_rum_dis_stats{iNetwork+1,3} = decreaseCount;
        intersubVarWithin_rum_dis_stats{iNetwork+1,4} = equalCount;
        
        %write out csv files for violin plots
        rawFC = y;
        rawFC(1:nSub,2) = 1;
        rawFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_within',within_networkFCName{iNetwork},'_rum_dis_rawFC.csv'], ...
            rawFC,'delimiter',',');
        adjustedFC = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
        adjustedFC(1:nSub,2) = 1;
        adjustedFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_within',within_networkFCName{iNetwork},'_rum_dis_adjustFC.csv'], ...
            adjustedFC,'delimiter',',');
        
        % rum vs. rest
        y = [rum_within_networkFC(:,iNetwork);rest_within_networkFC(:,iNetwork)];
        OtherCovariatesMatrix = [headmotion{iSite}(:,3);headmotion{iSite}(:,1)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        rum_rest_within_stats{iNetwork+1,2} = TF_ForContrast;
        rum_rest_within_stats{iNetwork+1,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        rum_rest_within_stats{iNetwork+1,4} = rum_rest_within_stats{iNetwork+1,3}*3;
        rum_rest_within_p4fdr(iNetwork,1) = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        
        %write out csv files for violin plots
        rawFC = y;
        rawFC(1:nSub,2) = 1;
        rawFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_within',within_networkFCName{iNetwork},'_rum_rest_rawFC.csv'], ...
            rawFC,'delimiter',',');
        adjustedFC = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
        adjustedFC(1:nSub,2) = 1;
        adjustedFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_within',within_networkFCName{iNetwork},'_rum_rest_adjustFC.csv'], ...
            adjustedFC,'delimiter',',');
        
        % rum vs. sad
        y = [rum_within_networkFC(:,iNetwork);sad_within_networkFC(:,iNetwork)];
        OtherCovariatesMatrix = [headmotion{iSite}(:,3);headmotion{iSite}(:,2)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        rum_sad_within_stats{iNetwork+1,2} = TF_ForContrast;
        rum_sad_within_stats{iNetwork+1,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        rum_sad_within_stats{iNetwork+1,4} = rum_sad_within_stats{iNetwork+1,3}*3;
        rum_sad_within_p4fdr(iNetwork,1) = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        
        %write out csv files for violin plots
        rawFC = y;
        rawFC(1:nSub,2) = 1;
        rawFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_within',within_networkFCName{iNetwork},'_rum_sad_rawFC.csv'], ...
            rawFC,'delimiter',',');
        adjustedFC = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
        adjustedFC(1:nSub,2) = 1;
        adjustedFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_within',within_networkFCName{iNetwork},'_rum_sad_adjustFC.csv'], ...
            adjustedFC,'delimiter',',');
        
        % sad vs. dis
        y = [sad_within_networkFC(:,iNetwork);dis_within_networkFC(:,iNetwork)];
        OtherCovariatesMatrix = [headmotion{iSite}(:,2);headmotion{iSite}(:,4)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        sad_dis_within_stats{iNetwork+1,2} = TF_ForContrast;
        sad_dis_within_stats{iNetwork+1,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        sad_dis_within_stats{iNetwork+1,4} = sad_dis_within_stats{iNetwork+1,3}*3;
        sad_dis_within_p4fdr(iNetwork,1) = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        
        %write out csv files for violin plots
        rawFC = y;
        rawFC(1:nSub,2) = 1;
        rawFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_within',within_networkFCName{iNetwork},'_sad_dis_rawFC.csv'], ...
            rawFC,'delimiter',',');
        adjustedFC = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
        adjustedFC(1:nSub,2) = 1;
        adjustedFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_within',within_networkFCName{iNetwork},'_sad_dis_adjustFC.csv'], ...
            adjustedFC,'delimiter',',');
    end
    
    %comparing betweennetwork FC
    rum_dis_between_stats = {'network','t','p','Bon_p','fdr_p'};
    for ii =2:length(between_networkFCName)+1
        rum_dis_between_stats{ii,1} = between_networkFCName{ii-1};
    end
    rum_rest_between_stats = {'network','t','p','Bon_p','fdr_p'};
    for ii =2:length(between_networkFCName)+1
        rum_rest_between_stats{ii,1} = between_networkFCName{ii-1};
    end
    rum_sad_between_stats = {'network','t','p','Bon_p','fdr_p'};
    for ii =2:length(between_networkFCName)+1
        rum_sad_between_stats{ii,1} = between_networkFCName{ii-1};
    end
    sad_dis_between_stats = {'network','t','p','Bon_p','fdr_p'};
    for ii =2:length(between_networkFCName)+1
        sad_dis_between_stats{ii,1} = between_networkFCName{ii-1};
    end
    intersubVarBetween_rum_dis_stats = {'','Increase','decrease','equal'};
    for ii = 2:length(between_networkFCName)+1
        intersubVarBetween_rum_dis_stats{ii,1} = between_networkFCName{ii-1};
    end
    
    for iNetwork = 1:size(rum_between_networkFC,2)
        % rum vs. dis
        y = [rum_between_networkFC(:,iNetwork);dis_between_networkFC(:,iNetwork)];
        OtherCovariatesMatrix = [headmotion{iSite}(:,3);headmotion{iSite}(:,4)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        rum_dis_between_stats{iNetwork+1,2} = TF_ForContrast;
        rum_dis_between_stats{iNetwork+1,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        rum_dis_between_stats{iNetwork+1,4} = rum_dis_between_stats{iNetwork+1,3}*3;
        rum_dis_between_p4fdr(iNetwork,1) = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        
        %As the reviewer required, I count the number of increase or
        %decreased
        increaseCount = 0;
        decreaseCount = 0;
        equalCount = 0;
        for ii = 1:nSub
            if y(ii) > y(ii+nSub)
                increaseCount = increaseCount+1;
            elseif y(ii) == y(ii+nSub)
                equalCount = equalCount+1;
            else
                decreaseCount = decreaseCount+1;
            end
        end
        intersubVarBetween_rum_dis_stats{iNetwork+1,2} = increaseCount;
        intersubVarBetween_rum_dis_stats{iNetwork+1,3} = decreaseCount;
        intersubVarBetween_rum_dis_stats{iNetwork+1,4} = equalCount;
        
        %write out csv files for violin plots
        rawFC = y;
        rawFC(1:nSub,2) = 1;
        rawFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_between',between_networkFCName{iNetwork},'_rum_dis_rawFC.csv'], ...
            rawFC,'delimiter',',');
        adjustedFC = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
        adjustedFC(1:nSub,2) = 1;
        adjustedFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_between',between_networkFCName{iNetwork},'_rum_dis_adjustFC.csv'], ...
            adjustedFC,'delimiter',',');
        
        % rum vs. rest
        y = [rum_between_networkFC(:,iNetwork);rest_between_networkFC(:,iNetwork)];
        OtherCovariatesMatrix = [headmotion{iSite}(:,3);headmotion{iSite}(:,1)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        rum_rest_between_stats{iNetwork+1,2} = TF_ForContrast;
        rum_rest_between_stats{iNetwork+1,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        rum_rest_between_stats{iNetwork+1,4} = rum_rest_between_stats{iNetwork+1,3}*3;
        rum_rest_between_p4fdr(iNetwork,1) = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        
         %write out csv files for violin plots
        rawFC = y;
        rawFC(1:nSub,2) = 1;
        rawFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_between',between_networkFCName{iNetwork},'_rum_rest_rawFC.csv'], ...
            rawFC,'delimiter',',');
        adjustedFC = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
        adjustedFC(1:nSub,2) = 1;
        adjustedFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_between',between_networkFCName{iNetwork},'_rum_rest_adjustFC.csv'], ...
            adjustedFC,'delimiter',',');
        
        % rum vs. sad
        y = [rum_between_networkFC(:,iNetwork);sad_between_networkFC(:,iNetwork)];
        OtherCovariatesMatrix = [headmotion{iSite}(:,3);headmotion{iSite}(:,2)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        rum_sad_between_stats{iNetwork+1,2} = TF_ForContrast;
        rum_sad_between_stats{iNetwork+1,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        rum_sad_between_stats{iNetwork+1,4} = rum_sad_between_stats{iNetwork+1,3}*3;
        rum_sad_between_p4fdr(iNetwork,1) = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        
         %write out csv files for violin plots
        rawFC = y;
        rawFC(1:nSub,2) = 1;
        rawFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_between',between_networkFCName{iNetwork},'_rum_sad_rawFC.csv'], ...
            rawFC,'delimiter',',');
        adjustedFC = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
        adjustedFC(1:nSub,2) = 1;
        adjustedFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_between',between_networkFCName{iNetwork},'_rum_sad_adjustFC.csv'], ...
            adjustedFC,'delimiter',',');
        
        %sad vs. dis
        y = [sad_between_networkFC(:,iNetwork);dis_between_networkFC(:,iNetwork)];
        OtherCovariatesMatrix = [headmotion{iSite}(:,2);headmotion{iSite}(:,4)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        sad_dis_between_stats{iNetwork+1,2} = TF_ForContrast;
        sad_dis_between_stats{iNetwork+1,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        sad_dis_between_stats{iNetwork+1,4} = sad_dis_between_stats{iNetwork+1,3}*3;
        sad_dis_between_p4fdr(iNetwork,1) = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        
        %write out csv files for violin plots
        rawFC = y;
        rawFC(1:nSub,2) = 1;
        rawFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_between',between_networkFCName{iNetwork},'_sad_dis_rawFC.csv'], ...
            rawFC,'delimiter',',');
        adjustedFC = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
        adjustedFC(1:nSub,2) = 1;
        adjustedFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_between',between_networkFCName{iNetwork},'_sad_dis_adjustFC.csv'], ...
            adjustedFC,'delimiter',',');
    end
    %fdr correction
    rum_dis_within_p4fdr = rum_dis_within_p4fdr(1:3,1);
    rum_rest_within_p4fdr = rum_rest_within_p4fdr(1:3,1);
    rum_sad_within_p4fdr = rum_sad_within_p4fdr(1:3,1);
    sad_dis_within_p4fdr = sad_dis_within_p4fdr(1:3,1);
    rum_dis_within_pfdr = mafdr(rum_dis_within_p4fdr,'BHFDR',true);
    rum_rest_within_pfdr = mafdr(rum_rest_within_p4fdr,'BHFDR',true);
    rum_sad_within_pfdr = mafdr(rum_sad_within_p4fdr,'BHFDR',true);
    sad_dis_within_pfdr = mafdr(sad_dis_within_p4fdr,'BHFDR',true);
    for ii = 1:3
        rum_dis_within_stats{ii+1,5} = rum_dis_within_pfdr(ii,1);
        rum_rest_within_stats{ii+1,5} = rum_rest_within_pfdr(ii,1);
        rum_sad_within_stats{ii+1,5} = rum_sad_within_pfdr(ii,1);
        sad_dis_within_stats{ii+1,5} = sad_dis_within_pfdr(ii,1);
    end
    rum_dis_between_p4fdr = rum_dis_between_p4fdr(:,1);
    rum_rest_between_p4fdr = rum_rest_between_p4fdr(:,1);
    rum_sad_between_p4fdr = rum_sad_between_p4fdr(:,1);
    sad_dis_between_p4fdr = sad_dis_between_p4fdr(:,1);
    rum_dis_between_pfdr = mafdr(rum_dis_between_p4fdr,'BHFDR',true);
    rum_rest_between_pfdr = mafdr(rum_rest_between_p4fdr,'BHFDR',true);
    rum_sad_between_pfdr = mafdr(rum_sad_between_p4fdr,'BHFDR',true);
    sad_dis_between_pfdr = mafdr(sad_dis_between_p4fdr,'BHFDR',true);
    iCount = 0;
    for ii = 1:3
        iCount = iCount+1;
        rum_dis_between_stats{ii+1,5} = rum_dis_between_pfdr(iCount,1);
        rum_rest_between_stats{ii+1,5} = rum_rest_between_pfdr(iCount,1);
        rum_sad_between_stats{ii+1,5} = rum_sad_between_pfdr(iCount,1);
        sad_dis_between_stats{ii+1,5} = sad_dis_between_pfdr(iCount,1);
    end
    
    fprintf('\ndoing correlation analysis\n\n');
    %%% correlation to scales and thinking contents
    % rum state FC corr
    subScalemMat = [ScaleMat,rum_thinking_content{iSite}(:,1:10)];
    subScaleNameSet = [{'Rumination','Brooding','Reflection'},ThinkingContentNameSet];
    %within network
    for iNetwork = 1:size(rum_within_networkFC,2)
        rum_within_corr_stats{iNetwork+1,1} = within_networkFCName{iNetwork};
        x = rum_within_networkFC(:,iNetwork);
        for iscale = 1:size(subScalemMat,2)
            y = subScalemMat(:,iscale);
            [r, p] = corr(x,y);
            rum_within_corr_stats{1,2*iscale} = subScaleNameSet{iscale};
            rum_within_corr_stats{1,2*iscale+1} = 'p';
            rum_within_corr_stats{iNetwork+1,2*iscale} = r;
            if p < 0.05
                rum_within_corr_stats{iNetwork+1,2*iscale+1} = p;
            else
                rum_within_corr_stats{iNetwork+1,2*iscale+1} = 0; 
            end
        end
    end
    %between network
    for iNetwork = 1:size(rum_between_networkFC,2)
        rum_between_corr_stats{iNetwork+1,1} = between_networkFCName{iNetwork};
        x = rum_between_networkFC(:,iNetwork);
        for iscale = 1:size(subScalemMat,2)
            y = subScalemMat(:,iscale);
            [r, p] = corr(x,y);
            rum_between_corr_stats{1,2*iscale} = subScaleNameSet{iscale};
            rum_between_corr_stats{1,2*iscale+1} = 'p';
            rum_between_corr_stats{iNetwork+1,2*iscale} = r;
            if p < 0.05
                rum_between_corr_stats{iNetwork+1,2*iscale+1} = p;
            else
                rum_between_corr_stats{iNetwork+1,2*iscale+1} = 0; 
            end
        end
    end
     
    %%%ROI level analysis
    %load FC matrixes
    restROI_FC = [];rumROI_FC = [];disROI_FC = [];
    for isubject = 1:length(SubList) 
        load([DataDir,'/',SiteSet{iSite},...
            restsitesuffix{iSite},'/Results/ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_',SubList{isubject},'.mat']);
        restROI_FC(:,:,isubject) = ROICorrelation_FisherZ;
        load([DataDir,'/',SiteSet{iSite}, ...
            tasksitesuffix{iSite},'/S2_Results/S2_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_',SubList{isubject},'.mat']);
        rumROI_FC(:,:,isubject) = ROICorrelation_FisherZ;
        load([DataDir,'/',SiteSet{iSite}, ...
            tasksitesuffix{iSite},'/S3_Results/S3_ROISignals_FunImgARCWFS/ROICorrelation_FisherZ_',SubList{isubject},'.mat']);
        disROI_FC(:,:,isubject) = ROICorrelation_FisherZ;
    end
    
    rum_dis_ROI_stats = {'ROI','t','p','fdr_p','rep'};
    rum_rest_ROI_stats = {'ROI','t','p','fdr_p','rep'};
    %output for heatmap plotting
    rum_dis_tmap = zeros(24,24);
    rum_rest_tmap = zeros(24,24);
    rum_dis_fdrpmap = ones(24,24);
    rum_rest_fdrpmap = ones(24,24);
    iCount = 0;
    rum_dis_p4fdr = [];
    rum_rest_p4fdr = [];
    for iROI = 1:23
       for jROI = iROI+1:24
           iCount = iCount +1;
           %Rum vs. Dis
           y = [squeeze(rumROI_FC(iROI,jROI,:));squeeze(disROI_FC(iROI,jROI,:))];
           OtherCovariatesMatrix = [headmotion{iSite}(:,3);headmotion{iSite}(:,4)];
           FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
           Contrast = zeros(1,size(FullRegressors,2));
           Contrast(1) = 1;
           TF_Flag = 'T';
           [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
           rum_dis_ROI_stats{iCount+1,1} = [ROIName{iROI},'-',ROIName{jROI}];
           rum_dis_ROI_stats{iCount+1,2} = TF_ForContrast;
           rum_dis_ROI_stats{iCount+1,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
           rum_dis_p4fdr(iCount,1) =  2*tcdf(-abs(TF_ForContrast),nSub-2);
           rum_dis_tmap(iROI,jROI) = TF_ForContrast;
           rum_dis_tmap(jROI,iROI) = TF_ForContrast;
           
           %Rum vs. Rest
           y = [squeeze(rumROI_FC(iROI,jROI,:));squeeze(restROI_FC(iROI,jROI,:))];
           OtherCovariatesMatrix = [headmotion{iSite}(:,3);headmotion{iSite}(:,1)];
           FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
           Contrast = zeros(1,size(FullRegressors,2));
           Contrast(1) = 1;
           TF_Flag = 'T';
           [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
           rum_rest_ROI_stats{iCount+1,1} = [ROIName{iROI},'-',ROIName{jROI}];
           rum_rest_ROI_stats{iCount+1,2} = TF_ForContrast;
           rum_rest_ROI_stats{iCount+1,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
           rum_rest_p4fdr(iCount,1) =  2*tcdf(-abs(TF_ForContrast),nSub-2);
           rum_rest_tmap(iROI,jROI) = TF_ForContrast;
           rum_rest_tmap(jROI,iROI) = TF_ForContrast;
       end
    end
    
    rum_dis_pfdr = mafdr(rum_dis_p4fdr,'BHFDR',true);
    rum_rest_pfdr = mafdr(rum_rest_p4fdr,'BHFDR',true);
    
    iCount = 0;
    for iROI = 1:23
       for jROI = iROI+1:24
           iCount = iCount +1;
           rum_dis_fdrpmap(iROI,jROI) = rum_dis_pfdr(iCount);
           rum_rest_fdrpmap(iROI,jROI) = rum_rest_pfdr(iCount);
           rum_dis_fdrpmap(jROI,iROI) = rum_dis_pfdr(iCount);
           rum_rest_fdrpmap(jROI,iROI) = rum_rest_pfdr(iCount);
       end
    end
    
    rum_dis_pfdr(find(rum_dis_pfdr>0.05)) = 0;
    eval(['rum_dis_ROI_mask',num2str(iSite),'=rum_dis_pfdr;']);
    eval(['rum_dis_ROI_mask',num2str(iSite),'(find(rum_dis_ROI_mask',num2str(iSite),'~=0))=1;']);
    
    rum_rest_pfdr(find(rum_rest_pfdr>0.05)) = 0;
    eval(['rum_rest_ROI_mask',num2str(iSite),'=rum_rest_pfdr;']);
    eval(['rum_rest_ROI_mask',num2str(iSite),'(find(rum_rest_ROI_mask',num2str(iSite),'~=0))=1;']);
    
    iCount = 0;
    for iROI = 1:23
       for jROI = iROI+1:24
           iCount = iCount +1;
           rum_dis_ROI_stats{iCount+1,4} = rum_dis_pfdr(iCount);
           rum_rest_ROI_stats{iCount+1,4} = rum_rest_pfdr(iCount);
       end
    end
    
    rim(iSite,1) = ceil(max(max(rum_dis_tmap)));
    rim(iSite,2) = floor(min(min(rum_dis_tmap)));
%   csvwrite([ResultDir,'/',SiteSet{iSite},'_rum_dis_tmap.csv'],rum_dis_tmap);
%   csvwrite([ResultDir,'/',SiteSet{iSite},'_rum_dis_fdrpmap.csv'],rum_dis_fdrpmap);
%   csvwrite([ResultDir,'/',SiteSet{iSite},'_rum_rest_tmap.csv'],rum_rest_tmap);
%   csvwrite([ResultDir,'/',SiteSet{iSite},'_rum_rest_fdrpmap.csv'],rum_rest_fdrpmap);
    
    dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_rum_dis_tmap.csv'],rum_dis_tmap,'delimiter',',');
    dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_rum_dis_fdrpmap.csv'],rum_dis_fdrpmap,'delimiter',',');
    dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_rum_rest_tmap.csv'],rum_rest_tmap,'delimiter',',');
    dlmwrite([ResultDir,'/csvfiles/',SiteSet{iSite},'_rum_rest_fdrpmap.csv'],rum_rest_fdrpmap,'delimiter',',');
    
    save([ResultDir,'/',SiteSet{iSite},'_Stats.mat'],'*stats');
end

% %get reproducible ROI level results
% rum_dis_ROI_Results = rum_dis_ROI_mask1.*rum_dis_ROI_mask2.*rum_dis_ROI_mask3;
% rum_rest_ROI_Results = rum_rest_ROI_mask1.*rum_rest_ROI_mask2.*rum_rest_ROI_mask3;
% 
% %get reproducible ROI level scale corr results
% rum_ROI_corr_Results = rum_ROIcorr_mask1.*rum_ROIcorr_mask2.*rum_ROIcorr_mask3;
% rum_dis_ROI_corr_Results = rum_dis_ROIcorr_mask1.*rum_dis_ROIcorr_mask2.*rum_dis_ROIcorr_mask3;
% rum_rest_ROI_corr_Results = rum_rest_ROIcorr_mask1.*rum_rest_ROIcorr_mask2.*rum_rest_ROIcorr_mask3;
% rest_ROI_corr_Results = rest_ROIcorr_mask1.*rest_ROIcorr_mask2.*rest_ROIcorr_mask3;

% for iSite = 1:length(SiteSet)
%     load([ResultDir,'/',SiteSet{iSite},'_Stats.mat']);
%     iCount = 0;
%     for iROI = 1:23
%        for jROI = iROI+1:24
%            iCount = iCount +1;
%            rum_dis_ROI_stats{iCount+1,5} = rum_dis_ROI_Results(iCount,1);
%            rum_rest_ROI_stats{iCount+1,5} = rum_rest_ROI_Results(iCount,1);
%            for iscale = 1:size(subScalemMat,2)
%                 rum_ROI_corr_stats{iCount,2*iscale+1} = rum_ROI_corr_stats{iCount,2*iscale+1}.* ...
%                     rum_ROI_corr_Results(iCount,2*iscale+1);
%                 rum_dis_ROI_corr_stats{iCount,2*iscale+1} = rum_dis_ROI_corr_stats{iCount,2*iscale+1}.* ...
%                     rum_dis_ROI_corr_Results(iCount,2*iscale+1);
%                 rum_rest_ROI_corr_stats{iCount,2*iscale+1} = rum_rest_ROI_corr_stats{iCount,2*iscale+1}.* ...
%                     rum_rest_ROI_corr_Results(iCount,2*iscale+1);
%                 rest_ROI_corr_stats{iCount,2*iscale+1} = rest_ROI_corr_stats{iCount,2*iscale+1}.* ...
%                     rest_ROI_corr_Results(iCount,2*iscale+1);
%            end
%        end
%     end
%     save([ResultDir,'/',SiteSet{iSite},'_Stats.mat'],'*stats');   
% end

% %For ploting ROI level reproducible connection matrix
% rum_dis_nodesize_weight = zeros(24,1);
% rum_dis_nodalconnec_matrix = zeros(24,24);
% rum_dis_nodalconnec_matrix_p = zeros(24,24);
% rum_dis_nodalconnec_matrix_n = zeros(24,24);
% iCount = 0;
% for iROI = 1:23
%    for jROI = iROI+1:24
%         iCount = iCount +1;
%         if rum_dis_ROI_Results(iCount) == 1
%             rum_dis_nodesize_weight(iROI) = 1;
%             rum_dis_nodesize_weight(jROI) = 1;
%             if rum_dis_ROI_stats{iCount+1,2} > 0
%                 rum_dis_nodalconnec_matrix(iROI,jROI) = 1;
%                 rum_dis_nodalconnec_matrix_p(iROI,jROI) = 1;
%             else
%                 rum_dis_nodalconnec_matrix(iROI,jROI) = -1;
%                 rum_dis_nodalconnec_matrix_n(iROI,jROI) = -1;
%             end
%         end
%    end
% end
% nodeInfo = [NodeCoord_flag,rum_dis_nodesize_weight];
% nodeInfo = num2cell(nodeInfo);
% nodeInfo = [nodeInfo,NodeName];
% [m,n] = size(nodeInfo);
% fid = fopen([ResultDir,'/csvfiles/rum_dis_RepConnect.node'],'w');
% for i = 1:m
%     fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%s\n',nodeInfo{i,:});   
% end
% fclose(fid);
% dlmwrite([ResultDir,'/csvfiles/rum_dis_RepConnect_p.edge'],rum_dis_nodalconnec_matrix_p,'delimiter','\t');
% dlmwrite([ResultDir,'/csvfiles/rum_dis_RepConnect_n.edge'],rum_dis_nodalconnec_matrix_n,'delimiter','\t');
% allrim(1,1) = max(rim(:,1));
% allrim(1,2) = min(rim(:,2));
% csvwrite([ResultDir,'/csvfiles/rim.csv'],allrim);

% addpath(genpath('/mnt/Data/RfMRILab/ChenX/CX_software/BrainNetViewer_20181031/'));
% surfacefile = '/mnt/Data/RfMRILab/ChenX/CX_software/BrainNetViewer_20181031/Data/SurfTemplate/ ...
%                    BrainMesh_ICBM152_smoothed.nv';
% nodefile = [ResultDir,'/rum_dis_RepConnect.node'];
% edgefile = [ResultDir,'/rum_dis_RepConnect_p.edge'];
% volumefile = '/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/DMN_subsystem.nii';
% settingfile = '/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/OptimizedSetting.mat';
% savepath = [ResultDir,'/rum_dis_RepConnect_p.jpg'];
% 
% BrainNet_MapCfg(surfacefile,nodefile,edgefile,volumefile,settingfile,savepath);
% close all;
