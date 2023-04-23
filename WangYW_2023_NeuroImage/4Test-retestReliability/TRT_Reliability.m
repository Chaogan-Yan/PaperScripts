clear;clc;
% 420
load /mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/SubInfo/SubInfo_420.mat;

%% mask
MaskFile ='/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Mask/Overlap38810.nii';
[MaskData,MaskHeader] = y_Read(MaskFile);
MaskIndex = find(MaskData);

SiteUnique=unique(Site);
SiteRegressor=zeros(size(Site,1),length(SiteUnique));
for iSite=1:length(SiteUnique)
    SiteRegressor(find(Site==SiteUnique(iSite)),iSite)=1;
end
SiteRegressor=SiteRegressor(:,1:end-1); % refer the last site as the base

% PALMSettings.nPerm = 5000;
% PALMSettings.ClusterInference=1;
% PALMSettings.ClusterFormingThreshold=2.3;
% PALMSettings.TFCE=1;
% PALMSettings.FDR=0;
% PALMSettings.TwoTailed=1;
% PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'

% dir name
%%
%ResultsSet={'Results','Results_S2'};
ResultsSet = {'Results','S2_Results'};
IndexName ={'ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','ReHo_FunImgARCWF','FC_D142'};
Harmo_methods = {'raw','reg','adj','lmm','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','ICVAE','SMA'};
status  = {'without_ScannerRegressor'};


do_stat = 1;
do_harmonize_btw=0;
do_fdr_btw=1;
do_GRF_btw=1;
caldice_FDR =1;
caldice_GRF =1;

%GRF
VoxelPThreshold=0.001;
IsTwoTailed=1;
ClusterPThreshold=0.05;
%%
%parpool(5);
for i_ResultSet = 1:length(ResultsSet)
    for i_method = 10:numel(Harmo_methods)
        for i_Index = 1:length(IndexName)
            statOutDir=[ '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/', ResultsSet{i_ResultSet},'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{1}];
            mkdir(statOutDir);
            fprintf(' T without site as a regressor \n');
            Cov=[Sex,Age,Motion(:,i_ResultSet),ones(length(SubID),1)];
            
            Contrast=zeros(1,size(Cov,2));
            Contrast(1)=1;
            disp('stat between gender  \n');
            if ~strcmp(IndexName{i_Index},'FC_D142')
                DataDir = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/Niis/',ResultsSet{i_ResultSet},'/',Harmo_methods{i_method}];
                switch IndexName{i_Index}
                    case 'ReHo_FunImgARCWF'
                        MeasurePrefix = 'szReHoMap_'
                    case 'ALFF_FunImgARCW'
                        MeasurePrefix = 'szALFFMap_'
                    case 'fALFF_FunImgARCW'
                        MeasurePrefix = 'szfALFFMap_'
                    case 'DegreeCentrality_FunImgARCWF'
                        MeasurePrefix =  'szDegreeCentrality_PositiveWeightedSumBrainMap_'
                end
                
                FileList=[];
                OutputName=[statOutDir,'/MaleVsFemaleT.nii'];
                
                for iSub = 1:length(SubID)
                    FileList{iSub,1}=[DataDir,'/',IndexName{i_Index},'/',MeasurePrefix,SubID{iSub},'.nii'];
                end
                [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0);
                
            else
                datapath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',ResultsSet{i_ResultSet},'/',IndexName{i_Index},'_',Harmo_methods{i_method} ,'.mat'];
                dat = importdata(datapath);
                if isstruct(dat)
                    dat = struct2array(dat);
                end

                tMap = zeros(142);
                pMap = zeros(142);
                p_vector = [];
                t_vector = [];
                TF_Flag ="T";
                
                OutputName=[statOutDir,'/MaleVsFemaleT.mat'];
                
                for n_fc = 1:(142*141/2)
                    [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(dat(:,n_fc),Cov,Contrast,TF_Flag);
                    t_vector= [ t_vector;TF_ForContrast];
                    p_vector= [p_vector;2*(1-tcdf(abs(TF_ForContrast),length(SubID)-4))]; 
                end
                
            end
            %% FDR, GRF
            if do_fdr_btw==1
                if strcmp(IndexName{i_Index},'FC_D142')
                    Q = 0.05;
                    SortP=sort(p_vector);
                    V=length(SortP);
                    I=(1:V)';
                    cVID = 1;
                    cVN  = sum(1./(1:V));
                    P   = SortP(find(SortP <= I/V*Q/cVID, 1, 'last' ));
                    
                    Thresholded=zeros(size(p_vector));
                    if ~isempty(P)
                        Thresholded(find(p_vector<=P))=1;
                    end
                    BinarizedDir=['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/',ResultsSet{i_ResultSet},'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{1},'/FDR'];
                    mkdir(BinarizedDir);
                    fdrOutputName =[BinarizedDir,'/MaleVsFemaleT_FDRBinarized.mat'];
                    save(fdrOutputName,'Thresholded');
                else
                    OutDir=[statOutDir,'/FDR'];
                    mkdir(OutDir)
                    cd(OutDir)
                    
                    BinarizedDir=['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/',ResultsSet{i_ResultSet},'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{1},'/FDR'];
                    mkdir(BinarizedDir);
                    
                    %get header from statistic result .nii
                    [Data,Header] = y_Read(strcat(statOutDir,'/MaleVsFemaleT.nii'));
                    SMap = Data(MaskIndex);
                    
                    % FDR
                    Q = 0.05;
                    Header0 = Header;
                    Header = w_ReadDF(Header);
                    
                    switch upper(Header.TestFlag)
                        case 'T'
                            PMap=2*(1-tcdf(abs(SMap), Header.Df));
                        case 'R'
                            PMap=2*(1-tcdf(abs(SMap).*sqrt((Header.Df)./(1-SMap.*SMap)), Header.Df));
                        case 'F'
                            PMap=(1-fcdf(SMap, Header.Df, Header.Df2));
                        case 'Z'
                            PMap=2*(1-normcdf(abs(SMap)));
                    end
                    SortP=sort(PMap); % sort ï¼?arrange order)
                    V=length(SortP);
                    I=(1:V)';
                    cVID = 1;
                    cVN  = sum(1./(1:V));
                    P   = SortP(find(SortP <= I/V*Q/cVID, 1, 'last' ));
                    
                    Thresholded=zeros(size(PMap));
                    if ~isempty(P)
                        Thresholded(find(PMap<=P))=1;
                    end
                    
                    AllBrain = zeros(61,73,61);
                    AllBrain = reshape(AllBrain,1,[]);
                    AllBrain(MaskIndex) = Thresholded;
                    AllBrain = reshape(AllBrain,61,73,61);
                    y_Write(Data.*AllBrain,Header,'MaleVsFemaleT');
                    
                    y_Write(AllBrain,Header,[BinarizedDir,'/MaleVsFemaleT']);
                end
            end
            if do_GRF_btw ==1 && ~strcmp(IndexName{i_Index},'FC_D142') % only for voxel-wise metrics
                OutDir=[statOutDir,'/GRF/ClusterCorrected001_05/'];
                mkdir(OutDir)
                cd(OutDir)
                
                % GRF CORRECTION
                FileName = [statOutDir,'/MaleVsFemaleT.nii'];
                OutName=[OutDir,'/MaleVsFemaleT'];
                y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutName,MaskFile);
                
                % READ GRF_AFTER RESULT and binary it
                GRF_DataDir=[OutDir,'/ClusterThresholded_MaleVsFemaleT'];
                BinarizedDir=['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/',ResultsSet{i_ResultSet},'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{1},'/GRF/ClusterCorrected001_05/'];
                mkdir(BinarizedDir);
                
                %get header from statistic result .nii
                [Data,Header] = y_Read(GRF_DataDir);
                ZMap = reshape(Data,1,[]);
                
                Thresholded=zeros(size(ZMap));
                if ~isempty(ZMap)
                    Thresholded(find(ZMap))=1;
                end
                
                AllBrain = reshape(Thresholded,61,73,61);
                
                y_Write(AllBrain,Header,[BinarizedDir,'/MaleVsFemaleT']);
            end
        end
    end
end
%% call dice
        DataDir ='/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/Results';
        DataDir2 ='/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/S2_Results';
for i_Index = 1:length(IndexName)
    
    if caldice_FDR == 1
        Dice= cell(length(Harmo_methods),1);
        num= cell(length(Harmo_methods),1);

        if ~strcmp(IndexName{i_Index},'FC_D142')
            for i_method = 1:length(Harmo_methods)   
                File=[DataDir,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{1},'/FDR/MaleVsFemaleT.nii'];
                Data=y_Read(File);
                File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{1},'/FDR/MaleVsFemaleT.nii'];
                Data2= y_Read(File);
                Data=Data+Data2;
                FDR_mask = Data;
                
                save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Reproducibility/',IndexName{i_Index},'_FDR_mask'],'FDR_mask');
                num{i_method}= length(find(Data==2));
                Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
            end
        else
  
            for i_method = 1:length(Harmo_methods)
                File=[DataDir,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{1},'/FDR/MaleVsFemaleT_FDRBinarized.mat'];
                Data=importdata(File);
                File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{1},'/FDR/MaleVsFemaleT_FDRBinarized.mat'];
                Data2= importdata(File);
                Data=Data+Data2;
                
                num{i_method}= length(find(Data==2));
                Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
            end
        end
        FDR_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
        csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Reproducibility/FDR_DICE.csv'],FDR_RESULT);
        %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Reproducibility/',IndexName{i_Index},'_FDR_DICE'],'FDR_RESULT');
    end
    
    if caldice_GRF == 1 && ~strcmp(IndexName{i_Index},'FC_D142')
        Dice= cell(length(Harmo_methods),1);
        num= cell(length(Harmo_methods),1);

        for i_method = 1:length(Harmo_methods)   
            File=[DataDir,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{1},'/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data=y_Read(File);
            File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{1},'/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data2= y_Read(File);
            Data=Data+Data2;
            GRF_mask = Data;
            
            save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Reproducibility/',IndexName{i_Index},'_GRF_mask'],'GRF_mask');
            num{i_method}= length(find(Data==2));
            Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
        end
        GRF_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
        csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Reproducibility/GRF_DICE.csv'],GRF_RESULT);
    end
end

