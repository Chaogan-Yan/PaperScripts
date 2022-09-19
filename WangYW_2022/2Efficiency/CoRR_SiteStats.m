%%
clear;clc;
addpath('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/codes4pub/CORR/utils');
IndexName ={'fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','ReHo_FunImgARCWF','ALFF_FunImgARCW'};
path = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR';
%Harmo_methods = {'raw','reg','adj','lmm','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','ICVAE','SMA'};
Harmo_methods = {'SMA'};

outputpath = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/Niis';
mkdir(outputpath);
rawniidir='/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/ResultsfromPaper/BeSess_SinRest';
MaskFile ='/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Mask/Overlap38810.nii'; 
load('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/SubInfo/SubInfo_420.mat');
ResultsSet={'Results','S2_Results'};
%parpool(2);
parfor ses =1:numel(ResultsSet)
    ses_path = [path,'/',ResultsSet{ses},'/'];
    ses_Rawdir = [rawniidir,'/',ResultsSet{ses},'S'];
    for i = 1:length(IndexName)
        outputdir = [outputpath, '/',ResultsSet{ses}];
        if ~exist(outputdir,'dir')
            mkdir(outputdir);
        end
        write_nii(ses_path,IndexName{i},Harmo_methods,outputdir,MaskFile,ses_Rawdir);
    end
end
%%
clear;clc;
load /mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/SubInfo/SubInfo_420.mat;

MaskFile ='/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Mask/Overlap38810.nii'; 
[MaskData,MaskHeader] = y_Read(MaskFile);
MaskIndex = find(MaskData);

SiteUnique=unique(Site);
SiteRegressor=zeros(size(Site,1),length(SiteUnique));
for iSite=1:length(SiteUnique)
    SiteRegressor(find(Site==SiteUnique(iSite)),iSite)=1;
end
SiteRegressor=SiteRegressor(:,1:end-1); % refer the last site as the base

%1 datadir
dir = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/';
%2 Results
ResultsSet = {'Results','S2_Results'};
%3 Method
Harmo_methods = {'raw','reg','adj','lmm','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','ICVAE','SMA'};
%Harmo_methods = {'SMA'};
%4 Metric
metrics = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
%5 prefix in switch function
statdir = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/Site';

do_stat = 1;
%GRF
VoxelPThreshold=0.001;
IsTwoTailed=0;
ClusterPThreshold=0.05;

%parpool(10);
for i_Index = 1:numel(metrics)
    if ~strcmp(metrics{i_Index},'FC_D142')
        PValue = zeros(38810,numel(Harmo_methods));
    else
        PValue = zeros(10011,numel(Harmo_methods));
    end
    for i_ResultSet = 1:numel(ResultsSet) % 4 organize data &hamonize
        Cov=[ones(length(SubID),1),Sex,Age,Motion(:,i_ResultSet),SiteRegressor];
        Contrast=zeros(1,size(Cov,2));
        Contrast(end-size(SiteRegressor,2)+1:end)=1;
        % according to index name to work
        VoxelDir =  [dir,'/Niis/',ResultsSet{i_ResultSet}];
        NetworkDir= [ dir,'/', ResultsSet{i_ResultSet}];
        for i_method = 1:10%numel(Harmo_methods)
            if strcmp(metrics{i_Index},'FC_D142')
                DataDir = NetworkDir;
                load('/mnt/Data3/RfMRILab/Wangyw/software/DPABI_V6.0_ForCamp/Templates/Dosenbach_Science_160ROIs_Info.mat');

             %% Network/voxel  statistics
                
                statOutDir=fullfile(statdir,ResultsSet{i_ResultSet},metrics{i_Index},Harmo_methods{i_method});
                mkdir(statOutDir);
                datapath = [NetworkDir,'/',metrics{i_Index},'_',Harmo_methods{i_method} ,'.mat'];
                dat = importdata(datapath);
                if isstruct(dat)
                    dat = struct2array(dat);
                end
                FileList=[];
                tMap = zeros(142);
                pMap = zeros(142);
                %                     tMap_420 = zeros(142);
                %                     pMap_420=zeros(142);
                p_vector = [];
                t_vector = [];
                TF_Flag ="F";
                
                OutputName=[statOutDir,'/SiteEffect.mat'];
                
                sprintf('stat between SITE /n');
                for n_fc = 1:(142*141/2)

                    [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(dat(:,n_fc),Cov,Contrast,TF_Flag);
                  
                    t_vector= [ t_vector;TF_ForContrast];
                    p_vector= [p_vector;1-fcdf(TF_ForContrast,length(SiteUnique)-1,numel(SubID)-size(Cov,2))];
                    
                    %fprintf('%d has been stat \n' , n);
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

            DataDir = [VoxelDir,'/',Harmo_methods{i_method}];

                switch metrics{i_Index}
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
                statOutDir=fullfile(statdir,ResultsSet{i_ResultSet},metrics{i_Index},Harmo_methods{i_method});
                mkdir(statOutDir);
                
                OutputName=[statOutDir,'/SiteEffect.nii'];
                sprintf('stat between sites /n');
                for iSub = 1:length(SubID)
                    FileList{iSub,1}=[DataDir,'/',metrics{i_Index},'/',MeasurePrefix,SubID{iSub},'.nii'];
                end
                [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'F',0);
                
                OutDir=statOutDir;
                cd(OutDir);
                %get header from statistic result .nii
                [Data,Header] = y_Read(OutputName);
                SMap = Data(MaskIndex);
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
                %if do_GRF_btw ==1  % only for voxel-wise metrics
                OutDir=[statOutDir,'/GRF/ClusterCorrected001_05/'];
                mkdir(OutDir)
                cd(OutDir)
                
                % GRF CORRECTION
                FileName = [statOutDir,'/SiteEffect.nii'];
                OutName=[OutDir,'/SiteEffect'];
                y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutName,MaskFile);
                
            end
        end
        dlmwrite([statdir,'/',ResultsSet{i_ResultSet},'/',metrics{i_Index},'_sitePvalue.csv'],PValue);
    end
end

%% visulization
clear;clc;
%1 datadir 
dir = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/Site';
%2 Results
ResultsSet = {'Results','S2_Results'};
%ResultsSet_FC = {'Results','Results_S2'};
%3 Method
methods = {'raw','reg','adj','lmm','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','SMA','ICVAE'};
%4 Metric
metrics = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF'};
%5  GRF
GRFfile = '/GRF/ClusterCorrected001_05/ClusterThresholded_SiteEffect.nii';
addpath('/mnt/Data3/RfMRILab/Wangyw/codes/matlab/functions/mymaps');
for i_ResultSet =1:numel(ResultsSet)
    PicDir =[ '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Pics/CoRR/',ResultsSet{i_ResultSet}];
    mkdir(PicDir);
    for i_metric =1:numel(metrics)
            for i_method = 1:numel(methods)
                [Data,~,Header] = y_ReadRPI(fullfile(dir,ResultsSet{i_ResultSet},metrics{i_metric},methods{i_method},GRFfile));
                maxP4eachmethod(i_method) = max(reshape(Data,[],1)); 
                if  maxP4eachmethod(i_method)~=0
                    minP4eachmethod(i_method) = min(reshape(Data(Data~=0),[],1));
                else
                    minP4eachmethod(i_method) = nan;
                end
            end
            NoSigMethodIndex{i_metric} = find(isnan(minP4eachmethod));
            log10_maxP{i_metric} = log10(max(maxP4eachmethod));% for following set the max/min
            log10_minP{i_metric} = log10(min(minP4eachmethod));
    end
    PMax = max(cell2mat(log10_maxP));
    PMin = min(cell2mat(log10_minP));
    NMax = 0;
    NMin = 0;
    for i_metric = 1:numel(metrics)

        SigScannerMethodIndex = num2cell(setdiff(1:10,NoSigMethodIndex{i_metric}));

        for i_method_index = 1:numel(SigScannerMethodIndex)
            [BrainVol,~,BrainHeader]  = y_ReadRPI(fullfile(dir,ResultsSet{i_ResultSet},metrics{i_metric},methods{SigScannerMethodIndex{i_method_index}},GRFfile));
            BrainVol_binary = BrainVol~=0;
            BrainVol  = log10(BrainVol).*double(BrainVol_binary);
            BrainVol( isnan(BrainVol)) = 0;

            ConnectivityCriterion = 18;
            ClusterSize = 0;
            UnderlayFileName='/mnt/Data3/RfMRILab/Wangyw/software/DPABI_V6.0_ForCamp/Templates/ch2.nii';
            viewtype='MediumView';

            [BrainNetViewerPath, fileN, extn] = fileparts(which('BrainNet.m'));
            SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];
            
            if SigScannerMethodIndex{i_method_index}==1
                ColorMap = mymap('Greens');%raw
            elseif SigScannerMethodIndex{i_method_index}==9
                ColorMap = mymap('Blues');%SMA
            elseif SigScannerMethodIndex{i_method_index}==10
                ColorMap = mymap('Purples');%ICVAE
            else
                ColorMap = y_AFNI_ColorMap(12);
                disp('there are other methods sig, please check!');
            end
        %

            H_BrainNet = y_CallBrainNetViewer(BrainVol,NMin,PMin,ClusterSize,ConnectivityCriterion,SurfFileName,viewtype,ColorMap,NMax,PMax,BrainHeader);

            pics_path =[PicDir,'/', metrics{i_metric},'/',methods{SigScannerMethodIndex{i_method_index}}];
            mkdir(pics_path);
            JPGFile=[pics_path,'/ScannerEffect_GRF.jpg'];
            eval(['print -r300 -djpeg -noui ''',JPGFile,''';']);
        end
    end
end
