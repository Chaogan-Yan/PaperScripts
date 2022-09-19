%%
clear;clc;
load /mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/SubInfo/SubInfo_420.mat;
%mask
MaskFile ='/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Mask/Overlap38810.nii';
[MaskData,MaskHeader] = y_Read(MaskFile);
MaskIndex = find(MaskData);

SiteUnique=unique(Site);
SiteUnique=unique(Site);
SiteRegressor=zeros(size(Site,1),length(SiteUnique));
for iSite=1:length(SiteUnique)
    SiteRegressor(find(Site==SiteUnique(iSite)),iSite)=1;
end
SiteRegressor=SiteRegressor(:,1:end-1); % refer the last site as the base

PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=1;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
%%
ResultsSet = {'Results','S2_Results'};
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
MeasurePrefixSet={'szReHoMap_','szALFFMap_','szfALFFMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_'};
status  = {'without_ScannerRegressor','with_ScannerRegressor'};

do_stat =1 ;
do_fdr_btw=1;
do_GRF_btw=1;
caldice_FDR =0;
caldice_GRF =0;

%GRF
VoxelPThreshold=0.001;
IsTwoTailed=1;
ClusterPThreshold=0.05;
StatDir = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/TargetSiteChoice';
for i_Index =1:length(IndexName)
    for i_ResultSet = 1:length(ResultsSet) % 4 organize data &hamonize
        for i_site =2:length(SiteUnique)
            for i_status = 1:1%length(status)
                if do_stat==1 & SiteUnique(i_site)~=7
                    Datadir = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/TargetSiteChoice/CoRR/',ResultsSet{i_ResultSet},'/Site',num2str(SiteUnique(i_site)),'_asT'];
                    %put path together? in order to change easily
                    statOutDir=[StatDir,'/',ResultsSet{i_ResultSet},'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/',status{i_status}];
                    mkdir(statOutDir);
                    % load demographic
                    load([Datadir,'/phenotypic_info.mat']);
                    % load data
                    SMA = importdata([Datadir,'/',IndexName{i_Index},'_SMA.mat']);
                    
                    if i_status == 1
                        fprintf(' T without site as a regressor \n');
                        Cov=[Sex,Age,Motion,ones(length(Subid),1)];
                    else
                        fprintf(' T with site as a regressor \n');
                        Cov=[Sex,Age,Motion,ones(length(Subid),1),SiteRegressor];
                    end
                    sprintf('stat between gender  /n');
                    TF_Flag ="T";
                    Contrast=zeros(1,size(Cov,2));
                    Contrast(1)=1;
                    if  strcmp(IndexName{i_Index},'FC_D142')
                        OutputName=[statOutDir,'/MaleVsFemaleT.mat'];
                       
                        p_vector = zeros(142*141/2,1);
                        t_vector = zeros(142*141/2,1);            
                        for n_fc = 1:(142*141/2)
                            [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(SMA(:,n_fc),Cov,Contrast,TF_Flag);
                            t_vector(n_fc,1)=TF_ForContrast;
                            p_vector(n_fc,1)=2*(1-tcdf(abs(TF_ForContrast),length(Subid)-4));
                            %fprintf('%d has been stat \n' , n);
                        end
                        index = find(tril(ones(142),-1));
                        tMap = zeros(142);
                        pMap = zeros(142);
                        tMap(index)  = t_vector; %need to change
                        pMap(index) = p_vector;
                        save(OutputName,'*Map');
                    else
                        if ~exist([Datadir,'/Niis/',IndexName{i_Index}],'dir') || ~exist([statOutDir,'/MaleVsFemaleT.nii'],'file')
                            
                            fprintf('write mask for nii \n')
                            FileList = [];
                            OutList = [];
                            mkdir([Datadir,'/Niis/',IndexName{i_Index},'/']);
                            for iSub=1:length(Subid)
                                FileList{iSub,1}=['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/Niis/',ResultsSet{i_ResultSet},'/raw/',IndexName{i_Index},'/',MeasurePrefixSet{i_Index},Subid{iSub},'.nii'];
                                [~,~,~,header] = y_ReadAll(FileList{iSub,1});
                                OutList{iSub,1}=[Datadir,'/Niis/',IndexName{i_Index},'/',MeasurePrefixSet{i_Index},Subid{iSub},'.nii'];
                                
                                header.fname = OutList{iSub,1};
                                blank = zeros(size(MaskData));
                                blank(MaskIndex) = SMA(iSub,:);
                                %eval(['blank(x) =  para_adj_', Harmo_methods{i_method},'(:,iSub);']);
                                %harmonized = reshape(blank,mask_size);
                                y_Write(blank,header,OutList{iSub,1});
                            end
                            FileList=[];
                            OutputName=[statOutDir,'/MaleVsFemaleT.nii'];
                            for iSub = 1:length(Subid)
                                FileList{iSub,1}=[Datadir,'/Niis/',IndexName{i_Index},'/',MeasurePrefixSet{i_Index},Subid{iSub},'.nii'];
                            end
                            [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'T',0);
                        end
                        if exist([Datadir,'/',IndexName{i_Index},'/NII'],'dir')
                            rmdir([Datadir,'/',IndexName{i_Index},'/NII']);
                            rmdir([Datadir,'/',IndexName{i_Index}]);
                        end
                    end
                end
                    %% FDR, GRF
                    if do_fdr_btw==1
                        if  strcmp(IndexName{i_Index},'FC_D142')
                            
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
                            fdrOutputName =[statOutDir,'/MaleVsFemaleT_FDRbinarized.mat'];
                            save(fdrOutputName,'Thresholded');
                        else
                            OutDir=[statOutDir,'/FDR'];
                            mkdir(OutDir)
                            cd(OutDir)
                            
                            BinarizedDir=[statOutDir,'/SignificantBinarized/FDR']
                            mkdir(BinarizedDir);
                            
                            %get header from statistic result .nii
                            [Data,Header] = y_Read(strcat(statOutDir,'/MaleVsFemaleT.nii'));
                            SMap = Data(MaskIndex);
                            
                            % FDR
                            Q = 0.05;
                            Header0 = Header;
                            Header = w_ReadDF(Header);
                            
                            % P-lization result
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
                            % Following  FDR.m	1.3 Tom Nichols 02/01/18
                            SortP=sort(PMap); % sort ?arrange order)
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
                    if do_GRF_btw ==1  % only for voxel-wise metrics
                        if strcmp(IndexName{i_Index},'FC_D142')
                            break;
                        else
                            OutDir=[statOutDir,'/GRF/ClusterCorrected001_05/'];
                            mkdir(OutDir)
                            cd(OutDir)
                            
                            % GRF CORRECTION
                            FileName = [statOutDir,'/MaleVsFemaleT.nii'];
                            OutName=[OutDir,'/MaleVsFemaleT'];
                            y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutName,MaskFile);
                            
                            % READ GRF_AFTER RESULT and binary it
                            GRF_DataDir=[OutDir,'/ClusterThresholded_MaleVsFemaleT'];
                            BinarizedDir=[statOutDir,'/SignificantBinarized/GRF/ClusterCorrected001_05/'];
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
    end
    
    
    %% call dice
    if caldice_FDR == 1
        Reproducibility=cell(length(SiteUnique),2);
        Jaccard= cell(length(SiteUnique),2);
        Dice= cell(length(SiteUnique),2);
        num= cell(length(SiteUnique),2);
        for i_site =2:length(SiteUnique)
            if SiteUnique(i_site)~=7
                if strcmp(IndexName{i_Index},'FC')
                    DataDir = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/ExperimentsOnTargetSiteSet/Results/Site',num2str(SiteUnique(i_site)),'_asT'];
                    DataDir2 = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/ExperimentsOnTargetSiteSet/Results_S2/Site',num2str(SiteUnique(i_site)),'_asT'];
                    for i_status = 1:2
                        File=[DataDir,'/',IndexName{i_Index},'/',status{i_status},'/MaleVsFemaleT_FDRbinarized.mat'];
                        Data=importdata(File);
                        File=[DataDir2,'/',IndexName{i_Index},'/',status{i_status},'/MaleVsFemaleT_FDRbinarized.mat'];
                        Data2= importdata(File);
                        Data=Data+Data2;
                        num{i_site,i_status}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
                        Reproducibility{i_site,i_status}=mean(Data(find(Data>0)));
                        Jaccard{i_site,i_status}=length(find(Data==2))/length(find(Data>=1));
                        Dice{i_site,i_status}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                    end
                else
                    DataDir = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/ExperimentsOnTargetSiteSet/ResultsS/Site',num2str(SiteUnique(i_site)),'_asT'];
                    DataDir2 = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/ExperimentsOnTargetSiteSet/S2_ResultsS/Site',num2str(SiteUnique(i_site)),'_asT'];
                    for i_status = 1:2
                        File=[DataDir,'/',IndexName{i_Index},'/',status{i_status},'/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                        Data=y_Read(File);
                        File=[DataDir2,'/',IndexName{i_Index},'/',status{i_status},'/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                        Data2= y_Read(File);
                        Data=Data+Data2;
                        num{i_site,i_status}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
                        Reproducibility{i_site,i_status}=mean(Data(find(Data>0)));
                        Jaccard{i_site,i_status}=length(find(Data==2))/length(find(Data>=1));
                        Dice{i_site,i_status}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                    end
                end
            end
        end
        FDR_RESULT = {Dice,num};
        save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/ExperimentsOnTargetSiteSet/Dice/',IndexName{i_Index},'_FDR_DICE'],'FDR_RESULT');
    end
    
    if caldice_GRF == 1
        Reproducibility=cell(length(SiteUnique),2);
        Jaccard= cell(length(SiteUnique),2);
        Dice= cell(length(SiteUnique),2);
        num= cell(length(SiteUnique),2);
        for i_site = 2:length(SiteUnique)
            if SiteUnique(i_site)~=7
                DataDir = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/ExperimentsOnTargetSiteSet/ResultsS/Site',num2str(SiteUnique(i_site)),'_asT'];
                DataDir2 = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/ExperimentsOnTargetSiteSet/S2_ResultsS/Site',num2str(SiteUnique(i_site)),'_asT'];
                for i_status = 1:2
                    File=[DataDir,'/',IndexName{i_Index},'/',status{i_status},'/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data=y_Read(File);
                    File=[DataDir2,'/',IndexName{i_Index},'/',status{i_status},'/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data2= y_Read(File);
                    Data=Data+Data2;
                    num{i_site,i_status}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
                    Reproducibility{i_site,i_status}=mean(Data(find(Data>0)));
                    Jaccard{i_site,i_status}=length(find(Data==2))/length(find(Data>=1));
                    Dice{i_site,i_status}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                end
            end
        end
        GRF_RESULT = {Dice,num};
        save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/ExperimentsOnTargetSiteSet/Dice/',IndexName{i_Index},'_GRF_DICE'],'GRF_RESULT');
    end

end


