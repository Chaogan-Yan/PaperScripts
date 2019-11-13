%Conducting a paired t-test (rum vs. rest and rum vs. dis) and one-way ANOVA to test
%whether FCs within or between networks are significantly different

%Written by ChenXiao 20180925
% 1 Rest 2 Happy 3 Sad 4 Rum 5 Dis
% 1, DN_core; 2, DN_dmPFC; 3, DN_MTL; 4, FPCN; 5, SN; 6, DMN_total

%modified in 20190213
%also compare dis and rest

%Prepare for publish
%Using GLM (y_regress) to conduct paired t-tests, including head motion as
%covarites

%ChenXiao 20190320

clear; clc;

Resultdir='/MD3860F/RfMRILab/ChenX/Rumination_project/Analysis/Analysis4Publish';
SiteSet = {'IPCAS','PKUGE','PKUSIMENS'};
SubList = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/Analysis/IPCAS_Sublist.txt');

%Extract head motion
%1,Rest; 2, Happy; 3, Sad; 4, Rum; 5, Dis;
HeadMotion = {};
for iSite = 1:length(SiteSet)
    for i=1:length(SubList)
        Temp=load(['/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/',SiteSet{iSite},'_rest/RealignParameter/',SubList{i},'/FD_Jenkinson_',SubList{i},'.txt']);
        HeadMotion{iSite}(i,1)=mean(Temp);
        Temp=load(['/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/',SiteSet{iSite},'_task/RealignParameter/',SubList{i},'/FD_Jenkinson_',SubList{i},'.txt']);
        HeadMotion{iSite}(i,2)=mean(Temp);
        Temp=load(['/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/',SiteSet{iSite},'_task/RealignParameter/',SubList{i},'/S2_FD_Jenkinson_',SubList{i},'.txt']);
        HeadMotion{iSite}(i,3)=mean(Temp);
        Temp=load(['/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/',SiteSet{iSite},'_task/RealignParameter/',SubList{i},'/S3_FD_Jenkinson_',SubList{i},'.txt']);
        HeadMotion{iSite}(i,4)=mean(Temp);
        Temp=load(['/MD3860F/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/',SiteSet{iSite},'_task/RealignParameter/',SubList{i},'/S4_FD_Jenkinson_',SubList{i},'.txt']);
        HeadMotion{iSite}(i,5)=mean(Temp);
    end
end


FDR_Within_rum_rest_Mask = {};
FDR_Within_rum_dis_Mask = {};
Between_rum_rest_FDRmask = {};
Between_rum_dis_FDRmask = {};

for iSite = 1:length(SiteSet)
    load([Resultdir,'/',SiteSet{iSite},'_Network_FC_across_brain_Smooth.mat']);

    %Extract FCs from different networks across all 5 mental states
    WithinNetworks_DN_core = zeros(41,5);
    WithinNetworks_DN_dmPFC = zeros(41,5);
    WithinNetworks_DN_MTL = zeros(41,5);
    WithinNetworks_FPCN = zeros(41,5);
    WithinNetworks_SN = zeros(41,5);
    WithinNetworks_DMN_total = zeros(41,5);

    for iCondition = 1:5
        WithinNetworks_DN_core(:,iCondition) = (AllWithin_Network_FC{iCondition}(:,1));
        WithinNetworks_DN_dmPFC(:,iCondition) = (AllWithin_Network_FC{iCondition}(:,2));
        WithinNetworks_DN_MTL(:,iCondition) = (AllWithin_Network_FC{iCondition}(:,3));
        WithinNetworks_FPCN(:,iCondition) = (AllWithin_Network_FC{iCondition}(:,4));
        WithinNetworks_SN(:,iCondition) = (AllWithin_Network_FC{iCondition}(:,5));
        WithinNetworks_DMN_total(:,iCondition) = (AllWithin_Network_FC{iCondition}(:,6));
    end
    Networkset = {WithinNetworks_DN_core, WithinNetworks_DN_dmPFC, WithinNetworks_DN_MTL, WithinNetworks_FPCN, WithinNetworks_SN, WithinNetworks_DMN_total};
    % Using y_regress_ss to conduct t test with GLM
    nSub = length(SubList);
    Regressors = [ones(nSub,1);-1*ones(nSub,1)];
    for i=1:nSub
        SubjectRegressors(:,i) = zeros(nSub*2,1);
        SubjectRegressors(i:nSub:nSub*2,i) = 1;
    end
    OtherCovariatesMatrix = [];
    for iNetwork = 1:6
        % rum vs. dis
        y = [Networkset{iNetwork}(:,4);Networkset{iNetwork}(:,5)];
        OtherCovariatesMatrix = [HeadMotion{iSite}(:,4);HeadMotion{iSite}(:,5)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        Within_rum_dis_tMap(1,iNetwork) = TF_ForContrast;
        Within_rum_dis_pMap(1,iNetwork) = 2*tcdf(-abs(TF_ForContrast),nSub-1);
        
        %rum vs. rest
        y = [Networkset{iNetwork}(:,4);Networkset{iNetwork}(:,1)];
        OtherCovariatesMatrix = [HeadMotion{iSite}(:,4);HeadMotion{iSite}(:,1)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        Within_rum_rest_tMap(1,iNetwork) = TF_ForContrast;
        Within_rum_rest_pMap(1,iNetwork) = 2*tcdf(-abs(TF_ForContrast),nSub-1);
        
        %dis vs. rest
        y = [Networkset{iNetwork}(:,5);Networkset{iNetwork}(:,1)];
        OtherCovariatesMatrix = [HeadMotion{iSite}(:,5);HeadMotion{iSite}(:,1)];
        FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
        Contrast = zeros(1,size(FullRegressors,2));
        Contrast(1) = 1;
        TF_Flag = 'T';
        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
        Within_dis_rest_tMap(1,iNetwork) = TF_ForContrast;
        Within_dis_rest_pMap(1,iNetwork) = 2*tcdf(-abs(TF_ForContrast),nSub-1);
    end

    FDR_Within_rum_rest_pMap = Within_rum_rest_pMap';
    FDR_Within_rum_rest_pMap = mafdr(FDR_Within_rum_rest_pMap,'BHFDR',true);
    FDR_Within_rum_rest_pMap = FDR_Within_rum_rest_pMap';
    %[Within_rum_rest_threshold,b] = FDR(FDR_Within_rum_rest_pMap,0.05);

    FDR_Within_rum_dis_pMap = Within_rum_dis_pMap';
    FDR_Within_rum_dis_pMap = mafdr(FDR_Within_rum_dis_pMap,'BHFDR',true);
    FDR_Within_rum_dis_pMap = FDR_Within_rum_dis_pMap';
    %[Within_rum_dis_threshold,b] = FDR(FDR_Within_rum_dis_pMap,0.05);
    
    FDR_Within_dis_rest_pMap = Within_dis_rest_pMap';
    FDR_Within_dis_rest_pMap = mafdr(FDR_Within_dis_rest_pMap,'BHFDR',true);
    FDR_Within_dis_rest_pMap = FDR_Within_dis_rest_pMap';

    Current_FDR_rum_rest_mask = zeros(1,6);
    Current_FDR_rum_rest_mask(find(FDR_Within_rum_rest_pMap < 0.05)) = 1;
%     if ~isempty(Within_rum_rest_threshold)
%         Current_FDR_rum_rest_mask(find(Within_rum_rest_pMap < Within_rum_rest_threshold)) = 1;
%     end

    Current_FDR_rum_dis_mask = zeros(1,6);
    Current_FDR_rum_dis_mask(find(FDR_Within_rum_dis_pMap < 0.05)) = 1;
%     if ~isempty(Within_rum_dis_threshold)
%         Current_FDR_rum_dis_mask(find(Within_rum_dis_pMap < Within_rum_dis_threshold)) = 1;
%     end


    FDR_Within_rum_rest_Mask{iSite} = Current_FDR_rum_rest_mask;
    FDR_Within_rum_dis_Mask{iSite} = Current_FDR_rum_dis_mask;

    %%%%%%%%%%%%%%%%%%%%%%%%Between Networks%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    BetweenNetworkFC = {};
    Between_rum_rest_PMap = [];
    Between_rum_dis_PMap = [];
    Between_rum_rest_TMap = [];
    Between_rum_dis_TMap = [];
    
    
    nSub = length(SubList);
    Regressors = [ones(nSub,1);-1*ones(nSub,1)];
    for i=1:nSub
        SubjectRegressors(:,i) = zeros(nSub*2,1);
        SubjectRegressors(i:nSub:nSub*2,i) = 1;
    end
    OtherCovariatesMatrix = [];
    
    for i = 1:5 %six networks
       for j = i+1:6 %six networks
           for iCondition = 1:5
               for iSub = 1:41
                    BetweenNetworkFC{iSub,iCondition} =  AllBetween_Network_FC{iCondition}(i,j,iSub);
               end
           end
           BetweenNetworkFC = cell2mat(BetweenNetworkFC);
           
            % rum vs. dis
            y = [BetweenNetworkFC(:,4);BetweenNetworkFC(:,5)];
            OtherCovariatesMatrix = [HeadMotion{iSite}(:,4);HeadMotion{iSite}(:,5)];
            FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
            Contrast = zeros(1,size(FullRegressors,2));
            Contrast(1) = 1;
            TF_Flag = 'T';
            [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
            Between_rum_dis_TMap(i,j) = TF_ForContrast;
            Between_rum_dis_PMap(i,j) = 2*tcdf(-abs(TF_ForContrast),nSub-1);

            %rum vs. rest
            y = [BetweenNetworkFC(:,4);BetweenNetworkFC(:,1)];
            OtherCovariatesMatrix = [HeadMotion{iSite}(:,4);HeadMotion{iSite}(:,1)];
            FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
            Contrast = zeros(1,size(FullRegressors,2));
            Contrast(1) = 1;
            TF_Flag = 'T';
            [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
            Between_rum_rest_TMap(i,j) = TF_ForContrast;
            Between_rum_rest_PMap(i,j) = 2*tcdf(-abs(TF_ForContrast),nSub-1);

            %dis vs. rest
            y = [BetweenNetworkFC(:,5);BetweenNetworkFC(:,1)];
            OtherCovariatesMatrix = [HeadMotion{iSite}(:,5);HeadMotion{iSite}(:,1)];
            FullRegressors = [Regressors,SubjectRegressors,OtherCovariatesMatrix];
            Contrast = zeros(1,size(FullRegressors,2));
            Contrast(1) = 1;
            TF_Flag = 'T';
            [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,FullRegressors,Contrast,TF_Flag);
            Between_dis_rest_TMap(i,j) = TF_ForContrast;
            Between_dis_rest_PMap(i,j) = 2*tcdf(-abs(TF_ForContrast),nSub-1);
           BetweenNetworkFC = {};
       end
    end

    icount = 0;
    for i = 1:5
        for j = i+1:6
            icount = icount +1;
            Between_rum_rest_p(icount,1) = Between_rum_rest_PMap(i,j);
            Between_rum_dis_p(icount,1) = Between_rum_dis_PMap(i,j);
            Between_dis_rest_p(icount,1) = Between_dis_rest_PMap(i,j);
        end
    end

        FDR_Between_rum_rest = mafdr(Between_rum_rest_p,'BHFDR',true);
        FDR_Between_rum_dis = mafdr(Between_rum_dis_p,'BHFDR',true);
        FDR_Between_dis_rest = mafdr(Between_dis_rest_p,'BHFDR',true);
%     [Between_rum_rest_threshold,b] = FDR(Between_rum_rest_p,0.05);
%     [Between_rum_dis_threshold,b] = FDR(Between_rum_dis_p,0.05);

        icount = 0;
        for i = 1:5
            for j = i+1:6
                icount = icount +1;
                Between_rum_rest_fdrcorrected_pMap(i,j) = FDR_Between_rum_rest(icount);
                Between_rum_dis_fdrcorrected_pMap(i,j) = FDR_Between_rum_dis(icount);
                Between_dis_rest_fdrcorrected_pMap(i,j) = FDR_Between_dis_rest(icount);
            end
        end


        id_rum_rest = find(Between_rum_rest_fdrcorrected_pMap < 0.05);
        id_rum_dis = find(Between_rum_dis_fdrcorrected_pMap < 0.05);

         mask = zeros(5,6);
%          if ~isempty(Between_rum_rest_threshold)
%             id_rum_rest = find(Between_rum_rest_PMap < Between_rum_rest_threshold);
%             mask(id_rum_rest) = 1;
%          end
         mask(id_rum_rest) = 1;
         Between_rum_rest_FDRmask{iSite} = mask;


         mask = zeros(5,6);
%          if ~isempty(Between_rum_dis_threshold)
%             id_rum_dis = find(Between_rum_dis_PMap < Between_rum_dis_threshold);
%             mask(id_rum_dis) = 1;
%          end
        mask(id_rum_dis) = 1;
        Between_rum_dis_FDRmask{iSite} = mask;


    %comparing Networks within each hemisphere

    % BetweenNetworkFC_L = {};
    % BetweenPMap_L = [];
    % BetweenFMap_L = [];
    % for i = 1:5
    %    for j = i+1:6
    %        for iCondition = 1:5
    %            for iSub = 1:41
    %                 BetweenNetworkFC_L{iSub,iCondition} =  AllBetween_Network_FC_L{iCondition}(i,j,iSub);
    %            end
    %        end
    %        BetweenNetworkFC_L = cell2mat(BetweenNetworkFC_L);
    %        t = table(BetweenNetworkFC_L(:,1), BetweenNetworkFC_L(:,2), BetweenNetworkFC_L(:,3), BetweenNetworkFC_L(:,4), BetweenNetworkFC_L(:,5),...
    %         'VariableNames',{'meas1','meas2','meas3','meas4','meas5'});
    %        Meas = table([1 2 3 4 5]', 'VariableNames', {'Measurements'});
    %        rm = fitrm(t, 'meas1-meas5~1','WithinDesign',Meas);
    %        ranovabl = ranova(rm);
    %        BetweenPMap_L(i,j) = ranovabl.pValue(1);
    %        BetweenFMap_L(i,j) = ranovabl.F(1);
    %        BetweenNetworkFC_L = {};
    %    end
    % end
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % BetweenNetworkFC_R = {};
    % BetweenPMap_R = [];
    % BetweenFMap_R = [];
    % for i = 1:5
    %    for j = i+1:6
    %        for iCondition = 1:5
    %            for iSub = 1:41
    %                 BetweenNetworkFC_R{iSub,iCondition} =  AllBetween_Network_FC_R{iCondition}(i,j,iSub);
    %            end
    %        end
    %        BetweenNetworkFC_R = cell2mat(BetweenNetworkFC_R);
    %        t = table(BetweenNetworkFC_R(:,1), BetweenNetworkFC_R(:,2), BetweenNetworkFC_R(:,3), BetweenNetworkFC_R(:,4), BetweenNetworkFC_R(:,5),...
    %         'VariableNames',{'meas1','meas2','meas3','meas4','meas5'});
    %        Meas = table([1 2 3 4 5]', 'VariableNames', {'Measurements'});
    %        rm = fitrm(t, 'meas1-meas5~1','WithinDesign',Meas);
    %        ranovabl = ranova(rm);
    %        BetweenPMap_R(i,j) = ranovabl.pValue(1);
    %        BetweenFMap_R(i,j) = ranovabl.F(1);
    %        BetweenNetworkFC_R = {};
    %    end
    % end

    save([Resultdir,'/',SiteSet{iSite},'_Stats_test_Smooth.mat'],'*Map*');

end

AllWithin_rum_rest_fdrmask = FDR_Within_rum_rest_Mask{1}.*FDR_Within_rum_rest_Mask{2};
AllWithin_rum_rest_fdrmask = AllWithin_rum_rest_fdrmask.*FDR_Within_rum_rest_Mask{3};

AllWithin_rum_dis_fdrmask = FDR_Within_rum_dis_Mask{1}.*FDR_Within_rum_dis_Mask{2};
AllWithin_rum_dis_fdrmask = AllWithin_rum_dis_fdrmask.*FDR_Within_rum_dis_Mask{3};


AllBetween_rum_rest_fdrmask = Between_rum_rest_FDRmask{1}.*Between_rum_rest_FDRmask{2};
AllBetween_rum_rest_fdrmask = AllBetween_rum_rest_fdrmask.*Between_rum_rest_FDRmask{3};

AllBetween_rum_dis_fdrmask = Between_rum_dis_FDRmask{1}.*Between_rum_dis_FDRmask{2};
AllBetween_rum_dis_fdrmask = AllBetween_rum_dis_fdrmask.*Between_rum_dis_FDRmask{3};

AllBetween_rum_dis_TMap = {};
AllWithin_rum_dis_TMap = {};
AllBetween_rum_rest_TMap = {};
AllWithin_rum_rest_TMap = {};
for iSite = 1:3
    load([Resultdir,'/',SiteSet{iSite},'_Stats_test_Smooth.mat']);
    AllBetween_rum_dis_TMap{iSite} = Between_rum_dis_TMap;
    AllWithin_rum_dis_TMap{iSite} = Within_rum_dis_tMap;
    AllBetween_rum_rest_TMap{iSite} = Between_rum_rest_TMap;
    AllWithin_rum_rest_TMap{iSite} = Within_rum_rest_tMap;
end

MeanBetween_rum_dis_TMap = (AllBetween_rum_dis_TMap{1}+AllBetween_rum_dis_TMap{2}+AllBetween_rum_dis_TMap{3})/3;
MeanWithin_rum_dis_TMap = (AllWithin_rum_dis_TMap{1}+AllWithin_rum_dis_TMap{2}+AllWithin_rum_dis_TMap{3})/3;

MeanBetween_rum_rest_TMap = (AllBetween_rum_rest_TMap{1}+AllBetween_rum_rest_TMap{2}+AllBetween_rum_rest_TMap{3})/3;
MeanWithin_rum_rest_TMap = (AllWithin_rum_rest_TMap{1}+AllWithin_rum_rest_TMap{2}+AllWithin_rum_rest_TMap{3})/3;

Reproducible_Between_rum_dis_TMap = MeanBetween_rum_dis_TMap.*AllBetween_rum_dis_fdrmask;
Reproducible_Between_rum_rest_TMap = MeanBetween_rum_rest_TMap.*AllBetween_rum_rest_fdrmask;

Reproducible_Within_rum_dis_TMap = MeanWithin_rum_dis_TMap.*AllWithin_rum_dis_fdrmask;
Reproducible_Within_rum_rest_TMap = MeanWithin_rum_rest_TMap.*AllWithin_rum_rest_fdrmask;

save([Resultdir,'/ReproducibleNetwork_maps_Smooth_FDR.mat'],'Reproducible*');


