% FCP site effect
clear;clc;

load /mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/FCP/demographic.mat;
SiteUnique=unique(Site);
SiteRegressor=zeros(size(Site,1),length(SiteUnique));
for iSite=1:length(SiteUnique)
    SiteRegressor(find(Site==SiteUnique(iSite)),iSite)=1;
end
SiteRegressor=SiteRegressor(:,1:end-1); % refer the last site as the base

MaskFile ='/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Mask/Overlap38810.nii';
[MaskData,MaskHeader] = y_Read(MaskFile);
MaskIndex = find(MaskData);

IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
MeasurePrefixSet={'szReHoMap_','szALFFMap_','szfALFFMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_'};
Harmo_methods = {'raw','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','SMA'};
status  = {'without_covariates','with_covariates'};

%GRF
VoxelPThreshold=0.001;
IsTwoTailed=0;
ClusterPThreshold=0.05;
OUTDIR = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/FCP/Site';
Dir =  '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/FCP/';
ResultsSet={'Results','S2_Results'};
do_fdr_btw=1;
do_GRF_btw=1;
%%
for ses = 1:2
    for i_Index = 1:length(IndexName)
        PValue = [];
        for i_method = 1:length(Harmo_methods)        
            for i_status = 2:length(status) %4 stat
                if i_status == 1
                    fprintf(' F without covariates \n');                    
                else
                    fprintf('  F with  covariates \n');
                    Cov=[ones(length(Subid),1),Sex,Age,Motion,SiteRegressor];
                end
                Contrast=zeros(1,size(Cov,2));
                Contrast(end-size(SiteRegressor,2)+1:end)=1;
                if strcmp(IndexName{i_Index},'FC_D142')
                    DataDir = fullfile(Dir,ResultsSet{ses});
                    load('/mnt/Data3/RfMRILab/Wangyw/software/DPABI_V6.0_ForCamp/Templates/Dosenbach_Science_160ROIs_Info.mat');
                    
                    %% Network/voxel  statistics
                    
                    statOutDir=[ OUTDIR,'/',ResultsSet{ses},'/', Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{i_status}];
                    mkdir(statOutDir);
                    datapath = [DataDir,'/',IndexName{i_Index},'_',Harmo_methods{i_method} ,'.mat'];
                    dat = importdata(datapath);
                    if isstruct(dat)
                        dat = struct2array(dat);
                    end
                    if size(dat,1)==600 
                        dat = dat(index587,:);
                    elseif size(dat,1)==1020
                        dat = dat(421:1020,:);
                        dat = dat(index587,:);
                    end

                    tMap = zeros(142);
                    pMap = zeros(142);
                    p_vector = [];
                    t_vector = [];
                    TF_Flag ="F";
                    
                    OutputName=[statOutDir,'/SiteEffect.mat'];
                    
                    sprintf('stat between SITE /n');
                    for n_fc = 1:(142*141/2)
                        [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(dat(:,n_fc),Cov,Contrast,TF_Flag);
                        
                        t_vector= [ t_vector;TF_ForContrast];
                        p_vector= [p_vector;1-fcdf(TF_ForContrast,length(SiteUnique)-1,numel(Subid)-size(Cov,2))];
                    end
                    index = find(tril(ones(142),-1));
                    pMap(index) = p_vector;
                    tMap(index) = t_vector;
                    save(OutputName,'*Map');
                    PValue(:,i_method)  = p_vector;
                    % fdr
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
                    if isempty(find(Thresholded))
                        fprintf(Harmo_methods{i_method});
                    else
                        temp_map = zeros(142);
                        temp_map(index) = Thresholded;
                        temp_map = temp_map+temp_map';
                        largescale_network = zeros(7);
                        for j =1:7
                            indexj = find(Dos_ExcludeCerebellum_142_YeoNetwork==j);
                            for k = 1:7
                                indexk =  find(Dos_ExcludeCerebellum_142_YeoNetwork==k);
                                largescale_network(j,k) = length(find(temp_map(indexj,indexk)));
                            end
                        end
                        fdrOutputName =[statOutDir,'/SiteEffect_FDRBinarized.mat'];
                        save(fdrOutputName,'Thresholded','largescale_network');
                    end
                    
                else
                    DataDir = fullfile(Dir,'Niis',ResultsSet{ses});
                    FileList=[];
                    for iSub = 1:length(Subid)
                        FileList{iSub,1}=[DataDir,'/',SiteName{iSub},'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',MeasurePrefixSet{i_Index},Subid{iSub},'.nii'];
                    end
                    
                    statOutDir=[ OUTDIR,'/',ResultsSet{ses},'/', Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{i_status}];
                    mkdir(statOutDir);
                    
                    
                    OutputName=[statOutDir,'/scannerF.nii'];
                    
                    sprintf('stat between scanners /n');
                    
                    [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'F',0);
                    
                    if do_fdr_btw==1
                        OutDir=[statOutDir,'/FDR'];
                        mkdir(OutDir)
                        cd(OutDir)
                        
                        BinarizedDir=[ OUTDIR,'/',ResultsSet{ses},'/SignificantBinarized/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{i_status},'/FDR'];
                        mkdir(BinarizedDir);
                        
                        %get header from statistic result .nii
                        [Data,Header] = y_Read(strcat(statOutDir,'/scannerF.nii'));
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
                        PValue(:,i_method) = PMap;
                        % Following  FDR.m	1.3 Tom Nichols 02/01/18
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
                        y_Write(Data.*AllBrain,Header,'scannerF');
                        
                        y_Write(AllBrain,Header,[BinarizedDir,'/scannerF']);
                    end
                    if do_GRF_btw ==1  % only for voxel-wise metrics
                        OutDir=[statOutDir,'/GRF/ClusterCorrected001_05/'];
                        mkdir(OutDir)
                        cd(OutDir)
                        
                        % GRF CORRECTION
                        FileName = [statOutDir,'/scannerF.nii'];
                        OutName=[OutDir,'/scannerF'];
                        y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutName,MaskFile);
                        
                        % READ GRF_AFTER RESULT and binary it
                        GRF_DataDir=[OutDir,'/ClusterThresholded_scannerF'];
                        BinarizedDir=[ OUTDIR,'/',ResultsSet{ses},'/SignificantBinarized/',Harmo_methods{i_method},'/',IndexName{i_Index},'/',status{i_status},'/GRF/ClusterCorrected001_05/'];
                        mkdir(BinarizedDir);
                        
                        %get header from statistic result .nii
                        [Data,Header] = y_Read(GRF_DataDir);
                        ZMap = reshape(Data,1,[]);
                        
                        Thresholded=zeros(size(ZMap));
                        if ~isempty(ZMap)
                            Thresholded(find(ZMap))=1;
                        end
                        
                        AllBrain = reshape(Thresholded,61,73,61);
                        
                        y_Write(AllBrain,Header,[BinarizedDir,'/scannerF']);
                    end
                end
            end
        end
       dlmwrite([OUTDIR,'/',ResultsSet{ses},'/',IndexName{i_Index},'_sitePvalue.csv'],PValue);

    end
end
