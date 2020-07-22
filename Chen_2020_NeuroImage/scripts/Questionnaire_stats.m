%conducting stats on the questionnaires at the end of each states
clear;clc;
SiteSet = {'IPCAS','PKUGE','PKUSIEMENS'};
ThinkingContentNameSet = {'Past','Future','Self','Others','Positive','Negative','Image','Speech','Happy','Sad'};
WorkDir = '/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/Analysis4Pub';
ResultDir='/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/Analysis4Pub/Scale_stats';
if ~exist(ResultDir); mkdir(ResultDir); end

EmotionFstats = {'','F','p'};
EmotionTstats = {'','sad-rest','p','rum-rest','p','dis-rest','p'};
for iSite = 1:length(SiteSet)
    load([WorkDir,'/',SiteSet{iSite},'_QuestionnaireData.mat']);
    %Emotion: before rest; after rest; sad; rum; dis
    Emotion = [];
    Emotion = cell2mat(EmotionScore);
    Emotion(:,3) = Emotion(:,4);
    Emotion(:,4) = Emotion(:,5);
    Emotion(:,5) = Emotion(:,6);
    Emotion(:,6) = [];

    t = table(Emotion(:,1), Emotion(:,2), Emotion(:,3), Emotion(:,4), Emotion(:,5),...
            'VariableNames',{'meas1','meas2','meas3','meas4','meas5'});
    Meas = table([1 2 3 4 5]', 'VariableNames', {'Measurements'});
    rm = fitrm(t, 'meas1-meas5~1','WithinDesign',Meas);
    ranovabl = ranova(rm);
    EmotionFstats{2,1} = SiteSet{iSite};
    EmotionFstats{2,2} = ranovabl.F(1);
    EmotionFstats{2,3} = ranovabl.pValue(1);

    EmotionTstats{2,1} = SiteSet{iSite};
    [H,P,CI,STATS] = ttest(Emotion(:,3),Emotion(:,2));
    EmotionTstats{2,3} = P;
    EmotionTstats{2,2} = STATS.tstat;
    [H,P,CI,STATS] = ttest(Emotion(:,4),Emotion(:,2));
    EmotionTstats{2,5} = P;
    EmotionTstats{2,4} = STATS.tstat;
    [H,P,CI,STATS] = ttest(Emotion(:,5),Emotion(:,2));
    EmotionTstats{2,7} = P;
    EmotionTstats{2,6} = STATS.tstat;
    
    %thinking content
    RumThinkingContent = cell2mat(RumThinkingContent);
    DisThinkingContent = cell2mat(DisThinkingContent);
    RestThinkingContent = cell2mat(RestThinkingContent);
    ThinkCotFstats = {'content','F','p'};
    ThinkCotTstats = {'content','rum-dis','corrected_p','rum-rest','corrected_p'};
    for i = 1:10
        t = table(RumThinkingContent(:,i), DisThinkingContent(:,i), RestThinkingContent(:,i),...
            'VariableNames',{'meas1','meas2','meas3'});
        Meas = table([1 2 3]', 'VariableNames', {'Measurements'});
        rm = fitrm(t, 'meas1-meas3~1','WithinDesign',Meas);
        ranovabl = ranova(rm);
        ThinkCotFstats{i+1,1} = ThinkingContentNameSet{i};
        ThinkCotFstats{i+1,2} = ranovabl.F(1);
        ThinkCotFstats{i+1,3} = ranovabl.pValue(1)*10;
        
        ThinkCotTstats{i+1,1} = ThinkingContentNameSet{i};
        [H,P,CI,STATS] = ttest(RumThinkingContent(:,i),DisThinkingContent(:,i));
        P = P*3;
        if P < 0.05
            ThinkCotTstats{i+1,3} = P;
        else
            ThinkCotTstats{i+1,3} = 0;
        end
        ThinkCotTstats{i+1,2} = STATS.tstat;
    
        [H, P,CI,STATS] = ttest(RumThinkingContent(:,i),RestThinkingContent(:,i));
        P = P*3;
        if P < 0.05
            ThinkCotTstats{i+1,5} = P;
        else
            ThinkCotTstats{i+1,5} = 0;
        end
        ThinkCotTstats{i+1,4} = STATS.tstat;
    end
    save([ResultDir,'/',SiteSet{iSite},'_Stats.mat'],'*stats');
end