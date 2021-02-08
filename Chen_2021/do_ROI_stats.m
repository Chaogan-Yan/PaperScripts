% The whole surface analysis found several ROIs on the PKUSIEMENS data
%
% this script first call y_ExtractROISignal_Surf to get the values from all
% these ROIs, then conducted paired t tests with these values on the PKUGE
% and IPCAS sites to see whether they can be generalized
%
% written by Xiao Chen on 20201229
% chenxiaochina@hotmail.com

clear; clc;
%initializing
ROI_Path = '/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/PairedT_Standardized_NoLeftHanded/Sig_ROIs';
Metric_Path = '/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Stability_Surf/Organized_Metrics_NoLeftHanded_Standardized';
csv_Dir = '/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/PairedT_Standardized_NoLeftHanded/Sig_ROIs_csvFiles';

ROI_Cell_lh = {'lh_DLPFC','lh_hippocampus'};
ROI_Cell_rh = {'rh_hippocampus', 'rh_MCC', 'rh_MPFC', 'rh_Occipital', 'rh_PostCentral', ...
                'rh_TPJdor', 'rh_TPJros'};
Site_Set = {'IPCAS','PKUGE', 'PKUSIEMENS'};
Scale_Set = {'Rumination', 'Brooding','Reflection'};
Condition_Set = {'rum', 'dis'};

Sub_List = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/Sublist_NoLeftHanded.txt');

%%%%%%%%%%%%%%%%%%% Extract roi signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dealing with the left hemisphere
ROIDef = {};
for iROI = 1:length(ROI_Cell_lh)
    eval([ROI_Cell_lh{iROI}, '= {};']);
end
for iSite = 1:length(Site_Set)
    for iCondition = 1:length(Condition_Set)
        for iROI = 1:length(ROI_Cell_lh)
            ROIDef{1,1} = [ROI_Path,'/PKUSIEMENS_',ROI_Cell_lh{iROI},'_Mask.gii'];
            for iSub = 1:length(Sub_List)
                ROISignals = [];
                Img = [Metric_Path,'/',Site_Set{iSite},'/',Condition_Set{iCondition},'/LH/lh/',Sub_List{iSub},'-LH.gii'];
                eval([ROI_Cell_lh{iROI},'{iSite}(iSub,iCondition) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);']);
            end
        end
    end
end

%dealing with the right hemisphere
ROIDef = {};
for iROI = 1:length(Site_Set)
    eval([ROI_Cell_rh{iROI}, '= {};']);
end
for iSite = 1:length(Site_Set)
    for iCondition = 1:length(Condition_Set)
        for iROI = 1:length(ROI_Cell_rh)
            ROIDef{1,1} = [ROI_Path,'/PKUSIEMENS_',ROI_Cell_rh{iROI},'_Mask.gii'];
            for iSub = 1:length(Sub_List)
                ROISignals = [];
                Img = [Metric_Path,'/',Site_Set{iSite},'/',Condition_Set{iCondition},'/RH/rh/',Sub_List{iSub},'-RH.gii'];
                eval([ROI_Cell_rh{iROI},'{iSite}(iSub,iCondition) = y_ExtractROISignal_Surf(Img, ROIDef, [], [], 0);']);
            end
        end
    end
end

% read in head motion
% 1: rum 2: dis
Head_Motion = {};
for iSite = 1:length(Site_Set)
    for i=1:length(Sub_List)
        Temp=load(['/mnt/Data3/RfMRILab/ChenX/Rumination_surf/',Site_Set{iSite},'_task/RealignParameter/sub-',Sub_List{i},'/S3_FD_Jenkinson_sub-',Sub_List{i},'.txt']);
        Head_Motion{iSite}(i,1)=mean(Temp);
        Temp=load(['/mnt/Data3/RfMRILab/ChenX/Rumination_surf/',Site_Set{iSite},'_task/RealignParameter/sub-',Sub_List{i},'/S4_FD_Jenkinson_sub-',Sub_List{i},'.txt']);
        Head_Motion{iSite}(i,2)=mean(Temp);
    end
end

%Using y_regress_ss to conduct t test with GLM
%build up result cell
Results = {};
Results{1,2} = 'IPCAS';
Results{1,4} = 'PKUGE';
Results{1,6} = 'PKUSIEMENS';
%generating subject regressors 
nSub = length(Sub_List);
Regressors = [ones(nSub,1);-1*ones(nSub,1)];
for i=1:nSub
    SubjectRegressors(:,i) = zeros(nSub*2,1);
    SubjectRegressors(i:nSub:nSub*2,i) = 1;
end

%dealing with the left hemisphere
OtherCovariatesMatrix = [];
for iSite = 1:length(Site_Set)
    for iROI = 1:length(ROI_Cell_lh)
        Results{iROI+1,1} = ROI_Cell_lh{iROI};
        eval(['y = [',ROI_Cell_lh{iROI},'{iSite}(:,1);',ROI_Cell_lh{iROI},'{iSite}(:,2)];']);
        OtherCovariatesMatrix = [];
        OtherCovariatesMatrix = [Head_Motion{iSite}(:,1);Head_Motion{iSite}(:,2)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        Results{iROI+1,iSite*2} = TF_ForContrast;
        Results{iROI+1,iSite*2+1} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        
        %write out csv files for violin plots
        adjustedFC = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
        adjustedFC(1:nSub,2) = 1;
        adjustedFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([csv_Dir,'/',Site_Set{iSite},'_',ROI_Cell_lh{iROI},'_rum_dis_adjustedStability.csv'], ...
            adjustedFC,'delimiter',',');
    end
end

%dealing with the right hemisphere
OtherCovariatesMatrix = [];
for iSite = 1:length(Site_Set)
    for iROI = 1:length(ROI_Cell_rh)
        Results{length(ROI_Cell_lh)+iROI+1,1} = ROI_Cell_rh{iROI};
        eval(['y = [',ROI_Cell_rh{iROI},'{iSite}(:,1);',ROI_Cell_rh{iROI},'{iSite}(:,2)];']);
        OtherCovariatesMatrix = [];
        OtherCovariatesMatrix = [Head_Motion{iSite}(:,1);Head_Motion{iSite}(:,2)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        Results{length(ROI_Cell_lh)+iROI+1,iSite*2} = TF_ForContrast;
        Results{length(ROI_Cell_lh)+iROI+1,iSite*2+1} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
        
        %write out csv files for violin plots
        adjustedFC = y - b(size(FullRegressors,2)).*FullRegressors(:,size(FullRegressors,2));
        adjustedFC(1:nSub,2) = 1;
        adjustedFC(nSub+1:2*nSub,2) = 2;
        dlmwrite([csv_Dir,'/',Site_Set{iSite},'_',ROI_Cell_rh{iROI},'_rum_dis_adjustedStability.csv'], ...
            adjustedFC,'delimiter',',');
    end
end

save('/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/ROI_stats_Standardized_NoLeftHanded.mat','Results');

%%%%%%%%%%%%%%%Do correlation analysis on scale measured scores%%%%%%%%%%%
% read scale data
load '/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Demographic_data/DemographicInfo_V3.mat';
ScaleMat = [Rumination,Brooding,Reflection];
DemographicMat = [Age,Sex];

%exclude left handed sub010
wantedSubMatrix = ones(length(Sublist),1);
%the variable Sublist here is stored in the "DemographicInfo_V2.mat", which contains the Sub010, not
%the Sub_List above!
wantedSubMatrix(10,1) = 0;
WantedSubIndex = find(wantedSubMatrix);
Sublist = Sublist(WantedSubIndex);
%now the Sublist is equal to Sub_List
ScaleMat = ScaleMat(WantedSubIndex,:);
DemographicMat = DemographicMat(WantedSubIndex,:);

% corr between stability during the rumination state and scale scores
%dealing with the left hemisphere
for iSite = 1:length(Site_Set)
    eval(['corr_RumStability_lh_',Site_Set{iSite}, '= {};']);
    eval(['corr_RumStability_lh_',Site_Set{iSite}, '{1,2} = "Rumination";']);
    eval(['corr_RumStability_lh_',Site_Set{iSite}, '{1,4} = "Brooding";']);
    eval(['corr_RumStability_lh_',Site_Set{iSite}, '{1,6} = "Reflection";']);
end

for iSite = 1:length(Site_Set)
    for iROI = 1:length(ROI_Cell_lh)
        eval(['corr_RumStability_lh_',Site_Set{iSite}, '{',num2str(iROI+1),',1} = "',ROI_Cell_lh{iROI},'";']);
        for iScale = 1:3
            eval(['y = ',ROI_Cell_lh{iROI},'{iSite}(:,1);']);% rumination state's stability
            x = ScaleMat(:,iScale);
            [r, p] = corr(x,y);
            %no correction
%             if p < 0.05
%                 eval(['corr_RumStability_lh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2+1),'} = p;']);
%             else
%                 eval(['corr_RumStability_lh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2+1),'} = 0;']);
%             end
            eval(['corr_RumStability_lh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2+1),'} = p;']);
            eval(['corr_RumStability_lh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2),'} = r;']);
        end
    end
end

%dealing with the right hemisphere
for iSite = 1:length(Site_Set)
    eval(['corr_RumStability_rh_',Site_Set{iSite}, '= {};']);
    eval(['corr_RumStability_rh_',Site_Set{iSite}, '{1,2} = "Rumination";']);
    eval(['corr_RumStability_rh_',Site_Set{iSite}, '{1,4} = "Brooding";']);
    eval(['corr_RumStability_rh_',Site_Set{iSite}, '{1,6} = "Reflection";']);
end

for iSite = 1:length(Site_Set)
    for iROI = 1:length(ROI_Cell_rh)
        eval(['corr_RumStability_rh_',Site_Set{iSite}, '{',num2str(iROI+1),',1} = "',ROI_Cell_rh{iROI},'";']);
        for iScale = 1:3
            eval(['y = ',ROI_Cell_rh{iROI},'{iSite}(:,1);']);% rumination state's stability
            x = ScaleMat(:,iScale);
            [r, p] = corr(x,y);
%             if p < 0.05
%                 eval(['corr_RumStability_rh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2+1),'} = p;']);
%             else
%                 eval(['corr_RumStability_rh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2+1),'} = 0;']);
%             end
            eval(['corr_RumStability_rh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2+1),'} = p;']);
            eval(['corr_RumStability_rh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2),'} = r;']);
            % write out csv files with for corr plot
            outputMatrix = [];
            outputMatrix = [x,y];
            dlmwrite([csv_Dir,'/corr_RumStability_rh_',Site_Set{iSite},'_',ROI_Cell_rh{iROI},'_',Scale_Set{iScale},'.csv'], ...
                outputMatrix,'delimiter',',');
        end
    end
end

% fdr correction
rum4fdr = {};
brooding4fdr = {};
ref4fdr = {};
for iSite = 1:length(Site_Set)
    for ii = 1:length(ROI_Cell_lh)
        eval(['rum4fdr{iSite}(ii,1) = corr_RumStability_lh_',Site_Set{iSite}, '{',num2str(ii+1),',3};']);
        eval(['brooding4fdr{iSite}(ii,1) = corr_RumStability_lh_',Site_Set{iSite}, '{',num2str(ii+1),',5};']);
        eval(['ref4fdr{iSite}(ii,1) = corr_RumStability_lh_',Site_Set{iSite}, '{',num2str(ii+1),',7};']);
    end
end

for iSite = 1:length(Site_Set)
    for ii = 1:length(ROI_Cell_rh)
        eval(['rum4fdr{iSite}(ii+length(ROI_Cell_lh),1) = corr_RumStability_rh_',Site_Set{iSite}, '{',num2str(ii+1),',3};']);
        eval(['brooding4fdr{iSite}(ii+length(ROI_Cell_lh),1) = corr_RumStability_rh_',Site_Set{iSite}, '{',num2str(ii+1),',5};']);
        eval(['ref4fdr{iSite}(ii+length(ROI_Cell_lh),1) = corr_RumStability_rh_',Site_Set{iSite}, '{',num2str(ii+1),',7};']);
    end
end

rum_fdr = {};
brooding_fdr = {};
ref_fdr = {};
for iSite = 1:length(Site_Set)
    for ii = 1:length(ROI_Cell_lh)+length(ROI_Cell_rh)
        rum_fdr{iSite} = mafdr(rum4fdr{iSite},'BHFDR',true);
    end
end


% corr between stability of rumination - distraction states and scale scores
%dealing with the left hemisphere
for iSite = 1:length(Site_Set)
    eval(['corr_RumMinusDisStability_lh_',Site_Set{iSite}, '= {};']);
    eval(['corr_RumMinusDisStability_lh_',Site_Set{iSite}, '{1,2} = "rumination";']);
    eval(['corr_RumMinusDisStability_lh_',Site_Set{iSite}, '{1,4} = "Brooding";']);
    eval(['corr_RumMinusDisStability_lh_',Site_Set{iSite}, '{1,6} = "Reflection";']);
end

for iSite = 1:length(Site_Set)
    for iROI = 1:length(ROI_Cell_lh)
        eval(['corr_RumMinusDisStability_lh_',Site_Set{iSite}, '{',num2str(iROI+1),',1} = "',ROI_Cell_lh{iROI},'";']);
        for iScale = 1:3
            % difference between rumination and distraction state's stability
            eval(['y = ',ROI_Cell_lh{iROI},'{iSite}(:,1) - ',ROI_Cell_lh{iROI},'{iSite}(:,2);']);
            x = ScaleMat(:,iScale);
            [r, p] = corr(x,y);
            if p < 0.05
                eval(['corr_RumMinusDisStability_lh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2+1),'} = p;']);
            else
                eval(['corr_RumMinusDisStability_lh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2+1),'} = 0;']);
            end
            eval(['corr_RumMinusDisStability_lh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2),'} = r;']);
        end
    end
end

%dealing with the right hemisphere
for iSite = 1:length(Site_Set)
    eval(['corr_RumMinusDisStability_rh_',Site_Set{iSite}, '= {};']);
    eval(['corr_RumMinusDisStability_rh_',Site_Set{iSite}, '{1,2} = "Rumination";']);
    eval(['corr_RumMinusDisStability_rh_',Site_Set{iSite}, '{1,4} = "Brooding";']);
    eval(['corr_RumMinusDisStability_rh_',Site_Set{iSite}, '{1,6} = "Reflection";']);
end

for iSite = 1:length(Site_Set)
    for iROI = 1:length(ROI_Cell_rh)
        eval(['corr_RumMinusDisStability_rh_',Site_Set{iSite}, '{',num2str(iROI+1),',1} = "',ROI_Cell_rh{iROI},'";']);
        for iScale = 1:3
            % difference between rumination and distraction state's stability
            eval(['y = ',ROI_Cell_rh{iROI},'{iSite}(:,1) - ',ROI_Cell_rh{iROI},'{iSite}(:,2);']);
            x = ScaleMat(:,iScale);
            [r, p] = corr(x,y);
            if p < 0.05
                eval(['corr_RumMinusDisStability_rh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2+1),'} = p;']);
            else
                eval(['corr_RumMinusDisStability_rh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2+1),'} = 0;']);
            end
            eval(['corr_RumMinusDisStability_rh_',Site_Set{iSite}, '{',num2str(iROI+1),',',num2str(iScale*2),'} = r;']);
        end
    end
end

save('/mnt/Data/RfMRILab/ChenX/Rumination_Stability/Analysis/ROI_stats/corrResults.mat','corr*');