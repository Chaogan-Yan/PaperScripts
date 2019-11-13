%To compute the ROI level's ISC (Inter-subject correlation) and ISFC
%First average the time courses of four segments from a mental state
%Then calculate the Pearson's correlation between a subject and the rest
%subjects' time courses

clc;clear all;
DataDir = '/mnt/Data/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing';
SubList = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/Analysis/IPCAS_Sublist.txt');
ResultDir = '/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis4Publish/ISC/ROI_ISC_CorrMaps/';
SiteSet = {'IPCAS','PKUGE','PKUSIMENS'};
%Resting state
%Compute the sum of all ROIs' time courses
AllSumTp = {};
for iSite = 1:length(SiteSet)
    SumTp = {};
    SumTp{1} = zeros(55,66);
    SumTp{2} = zeros(60,66);
    SumTp{3} = zeros(60,66);
    SumTp{4} = zeros(60,66);
    
    for iSub = 1:length(SubList)
        %load ROISignals
        ROISignals = [];
        load([DataDir,'/',SiteSet{iSite},'_rest/Results/ROISignals_FunImgARCWFS/ROISignals_',SubList{iSub},'.mat']);
        %add each segment's time course
        for iStimulus = 1:4
            if iStimulus == 1
                SumTp{iStimulus} = SumTp{iStimulus} + ROISignals(1:55,:);
            else
                SumTp{iStimulus} = SumTp{iStimulus} + ROISignals(iStimulus*60-64:iStimulus*60-5,:);
            end
        end
    end
    AllSumTp{iSite} = SumTp;
end
%Calculate the ISFC map
for iSite = 1:length(SiteSet)
    for iSub = 1:length(SubList)
        %load subjects ROI's signals
        MeanISCorrelationMap = [];
        ROISignals = [];
        load([DataDir,'/',SiteSet{iSite},'_rest/Results/ROISignals_FunImgARCWFS/ROISignals_',SubList{iSub},'.mat']);
        SubjectTp = {};
        RemTp = {};
        for iStimulus = 1:4
            if iStimulus == 1
                SubjectTp{iStimulus} = ROISignals(1:55,:);
                RemTp{iStimulus} = AllSumTp{iSite}{iStimulus} - ROISignals(1:55,:);
            else
                SubjectTp{iStimulus} = ROISignals(iStimulus*60-64:iStimulus*60-5,:);
                RemTp{iStimulus} = AllSumTp{iSite}{iStimulus} - ROISignals(iStimulus*60-64:iStimulus*60-5,:);
            end
        end
        for i = 1:66
            for j = 1:66
                for iStimulus = 1:4
                   ISCorrelation{iStimulus}(i,j) =  corr(SubjectTp{iStimulus}(:,i),RemTp{iStimulus}(:,j));
                end
            end
        end
        %Fisher's r to z
        for iStimulus = 1:4
            ISCorrelation{iStimulus} = 0.5 .* log( (1+ISCorrelation{iStimulus}) ./ (1-ISCorrelation{iStimulus}) );
        end
        MeanISCorrelationMap = (ISCorrelation{1}+ISCorrelation{2}+ISCorrelation{3}+ISCorrelation{4})./4;
        OutputDir = [ResultDir,'/Rest/',SiteSet{iSite}];
        if ~exist(OutputDir, 'dir'); mkdir(OutputDir); end
        save([OutputDir,'/',SubList{iSub},'_ISCMap.mat'],'MeanISCorrelationMap')
    end
end

%Four mental states
SessionSet = {'','S2_', 'S3_', 'S4_'}; % happy, sad, rum, dis
SessionNameSet = {'happy', 'sad', 'rum', 'dis'};
%read stimulus order info generated with Sort_behavioral_data.m:
%StimulusOrder_Site is ordered by 1:IPCAS, 2: PKUGE, 3: PKUSIMENS
load('/mnt/Data/RfMRILab/ChenX/Rumination_project/Data/Raw/Behavior_data/Stimulus_order/StimulusOrder.mat');
StimulusOrderSet = {IPCAS_StimulusOrder,PKUGE_StimulusOrder,PKUSIMENS_StimulusOrder};
%Calculate the sum time series of mental states
AllSumTp = {};
for iSite = 1:length(SiteSet)
    SumTp = cell(4);
    for iSession = 1:length(SessionSet)
        %a 4 x 4 cell: 4 conditions x 4 stimulus cell
        SumTp{iSession,1} = zeros(60,66);
        SumTp{iSession,2} = zeros(60,66);
        SumTp{iSession,3} = zeros(60,66);
        SumTp{iSession,4} = zeros(60,66);
        for iSub = 1:length(SubList)
            ROISignals = [];
            load([DataDir,'/',SiteSet{iSite},'_task/',SessionSet{iSession},'Results/',SessionSet{iSession},'ROISignals_FunImgARCWFS/ROISignals_',SubList{iSub},'.mat']);
            StimulusOrder = StimulusOrderSet{iSite};
            SubjectNum = str2num(SubList{iSub,1}(4:end));
            %iStimulus is the serial number of a stimulus
            %Order is the sequential number
            for iStimulus = 1:4
                Order = StimulusOrder{SubjectNum,iSession}(1,iStimulus);
                if Order == 1
                    SubjectTp{iStimulus} = ROISignals(1:55,:);
                else
                    SubjectTp{iStimulus} = ROISignals(Order*60-64:Order*60-5, :);
                end
            end
            for iStimulus = 1:4
                if size(SubjectTp{iStimulus},1) == 55
                    SumTp{iSession,iStimulus}(6:60,:) = SumTp{iSession,iStimulus}(6:60,:) + SubjectTp{iStimulus};
                else 
                    SumTp{iSession,iStimulus} = SumTp{iSession,iStimulus} + SubjectTp{iStimulus};
                end 
            end
        end
    end
    AllSumTp{iSite} = SumTp;
end

%calculate ISC map
for iSite = 1:length(SiteSet)
    for iSession = 1:length(SessionSet)
        for iSub = 1:length(SubList)
            SubjectTp = {};
            RemTp = {};
            ISCorrelation = {};
            ROISignals = [];
            load([DataDir,'/',SiteSet{iSite},'_task/',SessionSet{iSession},'Results/',SessionSet{iSession},'ROISignals_FunImgARCWFS/ROISignals_',SubList{iSub},'.mat']);
            StimulusOrder = StimulusOrderSet{iSite};
            SubjectNum = str2num(SubList{iSub,1}(4:end));
            for iStimulus = 1:4
                Order = StimulusOrder{SubjectNum,iSession}(1,iStimulus);
                if Order == 1
                    SubjectTp{iStimulus} = ROISignals(1:55,:);
                    RemTp{iStimulus} = AllSumTp{iSite}{iSession,iStimulus}(6:60,:) - ROISignals(1:55,:);
                else
                    SubjectTp{iStimulus} = ROISignals(Order*60-64:Order*60-5, :);
                    RemTp{iStimulus} = AllSumTp{iSite}{iSession,iStimulus} - ROISignals(Order*60-64:Order*60-5, :);
                end
            end
            for i = 1:66
                for j = 1:66
                    for iStimulus = 1:4
                        ISCorrelation{iStimulus}(i,j) =  corr(SubjectTp{iStimulus}(:,i),RemTp{iStimulus}(:,j));   
                    end
                end
            end
            %Fisher's r to z
            for iStimulus = 1:4
                ISCorrelation{iStimulus} = 0.5 .* log( (1+ISCorrelation{iStimulus}) ./ (1-ISCorrelation{iStimulus}) );
            end
            MeanISCorrelationMap = (ISCorrelation{1}+ISCorrelation{2}+ISCorrelation{3}+ISCorrelation{4})./4;
            OutputDir = [ResultDir,'/',SessionNameSet{iSession},'/',SiteSet{iSite},];
            if ~exist(OutputDir, 'dir'); mkdir(OutputDir); end
            save([OutputDir,'/',SubList{iSub},'_ISCMap.mat'],'MeanISCorrelationMap')
        end
    end
end