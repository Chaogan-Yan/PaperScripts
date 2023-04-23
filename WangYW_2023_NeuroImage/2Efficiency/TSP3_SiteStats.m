clear;clc;
%% nii file name combination
%1 datadir
dir = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/TSP3/Niis';
%2 Metric
metrics = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
%3 Method
methods = {'raw','reg','adj','lmm','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','MMD_pku_ge','ICVAE'};
%4 site name
sites = {'IPCAS','PKUGE','PKUSIEMENS'};
%5 prefix in switch function

%%  Analysis other Inputs
%1 mask
MaskFile = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Mask/Overlap38810.nii';
Mask = y_Read(MaskFile);
MaskIndex = find(Mask);
%2 predictors
SubID = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/TST/TSTafterDpabi/info/sublist.txt');
Age = xlsread('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/TST/TSTafterDpabi/info/subinfo.xlsx', 'F2:F42');
Gender = xlsread('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/TST/TSTafterDpabi/info/subinfo.xlsx', 'D2:D42');
Gender = Gender==1;
Motion = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/TST/TSTafterDpabi/HeadMotion.txt');
batchnum =3;
subnum=[41,41,41];
batch = [];
for i = 1:batchnum
    %         subnum(1,i) = input('subject number for sites\n ')
    batch = [batch;kron(i,ones(subnum(1,i),1))];
end
Site = batch;
SiteUnique=unique(Site);
SiteCov = dummyvar(Site);
SiteCov = [SiteCov(:,1),SiteCov(:,3)];
clear batchnum subnum batch Site;
Cov=[repmat(ones(length(SubID),1),3,1),repmat(Gender,3,1),repmat(Age,3,1),Motion,SiteCov];
%3 contrast
Contrast=zeros(1,size(Cov,2));
Contrast(end-size(SiteCov,2)+1:end)=1;
%4 outputdir
outputpath = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/TSP3';
% %% write nii
% addpath('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/codes4pub/TST/utils');
% for mc = 5:numel(metrics)
%     for md = 1:numel(methods)
%         datapath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/TSP3/ResultsS/',metrics{mc},'_',methods{md},'.mat']
%         data = importdata(datapath);
%         if isstruct(data)
%             data = struct2array(data);
%         end
%         for is = 1:numel(sites)
%             rawnii_path = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/TST/TSTafterDpabi/',sites{is},'/ResultsS'];
%             outputpath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/TSP3/Niis/',metrics{mc},'/',methods{md},'/',sites{is}];
%             mkdir(outputpath);
%             ind = [(is-1)*41+1:is*41];
%             write_nii(data(ind,:),SubID,outputpath,metrics{mc},MaskFile,rawnii_path);
%         end
%     end
% end
%% statistical settings
%GRF
VoxelPThreshold=0.001;
IsTwoTailed=0;
ClusterPThreshold=0.05;

do_savep_btw =1;
do_GRF_btw = 1;
%parpool(4);
for i_Index = 1:numel(metrics)
    for i_method =  1:length(methods)
        if strcmp(metrics{i_Index},'FC_D142')
            load('/mnt/Data3/RfMRILab/Wangyw/software/DPABI_V6.0_ForCamp/Templates/Dosenbach_Science_160ROIs_Info.mat');

            statOutDir = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/TSP3/',metrics{i_Index}];
            mkdir(statOutDir);
            DataDir = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/TSP3/ResultsS';       
            statOutDir=fullfile('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/TSP3/',metrics{i_Index},'/',methods{i_method});
            mkdir(statOutDir);
            datapath = [DataDir,'/',metrics{i_Index},'_',methods{i_method} ,'.mat'];
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
                fprintf(methods{i_method});
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
            for is = 1:numel(sites)
                for iSub = 1:numel(SubID)
                    temp = iSub+41*(is-1);
                    FileList{temp,1}=[dir,'/',metrics{i_Index},'/',methods{i_method},'/',sites{is},'/',MeasurePrefix,SubID{iSub},'.nii'];
                end
            end

            statOutDir=[ '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/TSP3/',metrics{i_Index},'/',methods{i_method},'/'];
            mkdir(statOutDir);
            
            OutputName=[statOutDir,'/ScannerEffect.nii'];
            sprintf('stat between scanners /n');
            
            [b_OLS_brain, t_OLS_brain, TF_ForContrast_brain, r_OLS_brain, Header, SSE_OLS_brain] = y_GroupAnalysis_Image(FileList,Cov,OutputName,MaskFile,[],Contrast,'F',0);
            
            if do_savep_btw==1
                OutDir=statOutDir;
                cd(OutDir);
                %get header from statistic result .nii
                [Data,Header] = y_Read(strcat(statOutDir,'/ScannerEffect.nii'));
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
                AllBrain = zeros(61,73,61);
                AllBrain(MaskIndex) = PMap;
                y_Write(Data.*AllBrain,Header,'ScannerEffect_p.nii');
                
                
            end
            
            if do_GRF_btw ==1  % only for voxel-wise metrics
                OutDir=[statOutDir,'/GRF/ClusterCorrected001_05/'];
                mkdir(OutDir)
                cd(OutDir)
                
                % GRF CORRECTION
                FileName = [statOutDir,'/ScannerEffect.nii'];
                OutName=[OutDir,'/ScannerEffect'];
                y_GRF_Threshold(FileName,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutName,MaskFile);
                
                % READ GRF_AFTER RESULT and binary it
                GRF_DataDir=[OutDir,'/ClusterThresholded_scannerF'];
                %BinarizedDir=['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/FCP_Organized/stat/scanner/WithCORR/SignificantBinarized/',methods{i_method},'/',metrics{i_Index},'/GRF/ClusterCorrected001_05/'];
                %mkdir(BinarizedDir);
                
                %get header from statistic result .nii
                %                 [Data,Header] = y_Read(GRF_DataDir);
                %                 ZMap = reshape(Data,1,[]);
                %
                %                 Thresholded=zeros(size(ZMap));
                %                 if ~isempty(ZMap)
                %                     Thresholded(find(ZMap))=1;
                %                 end
                
                %                 AllBrain = reshape(Thresholded,61,73,61);
                %
                %                 y_Write(AllBrain,Header,[BinarizedDir,'/ScannerEffect']);
            end
        end
    end
    dlmwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/TSP3/',metrics{i_Index},'_sitePvalue.csv'],PValue);
end

