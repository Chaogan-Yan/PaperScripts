%% SMA
clear;
clc;
info = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/CoRR/SubInfo/SubInfo_420.mat');
site = info.Site;
age = info.Age;
sex = info.Sex;
sex(sex==-1)=0;
motion1 = info.Motion(:,1);
motion2 = info.Motion(:,2);
subid = info.SubID;
SiteSWU_ind = find(site==8);
ResultsSet = {'Results','S2_Results'};


SexTarget = sex(SiteSWU_ind);
AgeTarget = age(SiteSWU_ind);
AgeTarget = double(AgeTarget<22)+1;
IndexTarget = AgeTarget;% + 2*double(SexTarget) ;
files = {'Results','S2_Results'};
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};

info = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/FCP_Organized/SubInfo/600_subinfo.mat');  % 1-female 2-male
site = info.SiteID;
age = double(info.Age);
sex = info.Sex-1; % 0-female 1-male
motion = info.MeanFDJ;
eye = info.eyeclosed;
sitename = info.SiteName;
subid = info.Subid;

BEIJING_female_inds = find(strcmp(sitename,"Beijing")& sex==0);
Cambridge_male_inds = find(strcmp(sitename,"Cambridge")& sex==1);
refer_site = {BEIJING_female_inds,Cambridge_male_inds};
inds = cell2mat(refer_site');
refer_sitename = sitename(inds);
    
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};

% overlap mask
Cov = [sex(inds),double(age(inds)),motion(inds),ones(180,1)];
Datadir = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Revision/Data/SexinSeparateSites';
subid = info.Subid(inds);

%parpool(15);
for ses =1:2
    savepath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Revision/Data/SexinSeparateSites/',files{ses}];
    mkdir(savepath);
    for i_Index =  2:2% length(IndexName)
        file = [savepath,'/',IndexName{i_Index},'_SMA_SWU_beijing_cambridge.mat'];
        if ~exist(file,'file')
            if strcmp(IndexName{i_Index},'FC')
                datapath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',files{ses},'/',IndexName{i_Index},'_raw.mat'];
            else
                datapath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',files{ses},'/',IndexName{i_Index},'_raw.mat'];
            end
            data = importdata(datapath);
            Xtarget = data(SiteSWU_ind,:);
            datapath =['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/FCP/Results/',IndexName{i_Index},'_raw.mat'];
            FCP_data = importdata(datapath);
            
            %% define target and source
            for i =1:length(refer_site)
                AgeSource =  age(refer_site{i});
                SexSource = sex(refer_site{i});
                Xsource =  FCP_data(refer_site{i},:);
                
                temp = zeros(size(Xsource));
                AgeSource = double(AgeSource<22)+1;
                IndexSource = AgeSource + 2*double(SexSource) ;
                %MatrixTran{i-1}= zeros(4,3);
                for k = 1:size(temp,2)
                    MatrixTran=[];
                    Source = Xsource(:,k);
                    Target = Xtarget(:,k);
                    [slope,intercept,pvalue] = subsamplingMMD(Source,Target,IndexSource,IndexTarget,100);
                    %[slope,intercept,pvalue] = fitMMD(Source,Target,0);
                    %MatrixTran(k,:) = [slope,intercept,pvalue];
                    temp(:,k) = slope*Source+ones(size(Source))*intercept;
                    %Here, pvalue is for subsamples from first iteration
                end
                new_site_data{i} = temp;
                fprintf('site %d fitting completed \n', i);
            end
            
            
            SMA = cell2mat(new_site_data');
            save([savepath,'/',IndexName{i_Index},'_SMA_SWU_beijing_cambridge.mat'],'SMA');
        else 
            load(file);
        end
        RawNiiPath = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/FCP_Organized/AllMaps/';
        MaskFile = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Mask/Overlap38810.nii';
        
        NII_outputpath = fullfile(Datadir,'/Niis/',ResultsSet{ses},'/',IndexName{i_Index});
        mkdir(NII_outputpath);
        StatOutDir =  fullfile(Datadir,'/Stat/SWU_beijing_cambridge/',ResultsSet{ses},'/',IndexName{i_Index});
        mkdir(StatOutDir);
        
        outputname = [StatOutDir,'/MaleVSFemaleT.nii'];
        
        GRFOutDir=[StatOutDir,'/GRF/ClusterCorrected001_05/'];
        mkdir(GRFOutDir);
        GRFOutName=[GRFOutDir,'/MaleVsFemaleT'];
        
                
        if ~exist(outputname,'file')
            % write nii
            StatList = writeSiteNii(SMA,IndexName{i_Index},refer_sitename,subid,NII_outputpath,MaskFile,RawNiiPath);
            
            % TWO-SAMPLE-TEST
            y_GroupAnalysis_Image(StatList,Cov,outputname,MaskFile,[],[1,0,0,0],'T',0);
        end
        % GRF correction
        if length(dir(GRFOutDir))==2 % dir is empty
            [Data_corrected,~,Header]= y_GRF_Threshold(outputname,0.001,1,0.05,GRFOutName,MaskFile);
        else
            Data_corrected = y_Read([GRFOutDir,'/ClusterThresholded_MaleVsFemaleT.nii']);
        end
        Data_corrected =Data_corrected(Data_corrected~=0);
        NMax = min(Data_corrected(Data_corrected<0));
        NMin = max(Data_corrected(Data_corrected<0));
        PMin = min(Data_corrected(Data_corrected>0));
        PMax = max(Data_corrected(Data_corrected>0));
        ConnectivityCriterion = 18;
        ClusterSize = 0;
        UnderlayFileName='/mnt/Data3/RfMRILab/Wangyw/software/DPABI_V6.0_ForCamp/Templates/ch2.nii';
        GRFfile = [GRFOutDir,'/ClusterThresholded_MaleVsFemaleT.nii']
        ColorMap = y_AFNI_ColorMap(12);
        [BrainNetViewerPath, fileN, extn] = fileparts(which('BrainNet.m'));
        SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];
        viewtype='MediumView';
        pics_path = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Revision/Data/SexinSeparateSites/Pics/',ResultsSet{ses}];
        
        mkdir(pics_path);
         H_BrainNet = y_CallBrainNetViewer(GRFfile,NMin,PMin,ClusterSize,ConnectivityCriterion,SurfFileName,viewtype,ColorMap,NMax,PMax );
        JPGFile=[pics_path,'/',IndexName{i_Index},'_BeijingF_CambM.jpeg'];
        eval(['print -r300 -djpeg -noui ''',JPGFile,''';']);
            
    end
end

%% ComBat/CovBat
clear;
clc;
info = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/CoRR/SubInfo/SubInfo_420.mat');
site = info.Site;
age_corr = info.Age;
sex = info.Sex;
sex(sex==-1)=0;
sex_corr=sex;
motion1 = info.Motion(:,1);
motion2 = info.Motion(:,2);
subid = info.SubID;
SiteSWU_ind = find(site==8);
file = {'Results','S2_Results'};
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};

info = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/FCP_Organized/SubInfo/600_subinfo.mat');  % 1-female 2-male
site = info.SiteID;
age = double(info.Age);
sex = info.Sex-1; % 0-female 1-male
motion = info.MeanFDJ;
eye = info.eyeclosed;
sitename = info.SiteName;
subid = info.Subid;

BEIJING_female_inds = find(strcmp(sitename,"Beijing")& sex==0);
Cambridge_male_inds = find(strcmp(sitename,"Cambridge")& sex==1);
refer_site = {BEIJING_female_inds,Cambridge_male_inds};
inds = cell2mat(refer_site');
refer_sitename = sitename(inds);
batch = repelem([1,2,3]',[221,106,74]);
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
mod = [[sex_corr(SiteSWU_ind);sex(inds)],[age_corr(SiteSWU_ind);double(age(inds))]];
% overlap mask
Cov = [sex(inds),double(age(inds)),motion(inds),ones(180,1)];
Datadir = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Revision/Data/SexinSeparateSites';
subid = info.Subid(inds);
parpool(4);
parfor ses =1:2
    savepath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Revision/Data/SexinSeparateSites/',file{ses}];
    mkdir(savepath);
    for i_Index =  2:2% length(IndexName)
        if strcmp(IndexName{i_Index},'FC')
            datapath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',file{ses},'/',IndexName{i_Index},'_raw.mat'];
        else
            datapath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',file{ses},'/',IndexName{i_Index},'_raw.mat'];
        end
        data = importdata(datapath);
        Xtarget = data(SiteSWU_ind,:);
        datapath =['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/FCP/Results/',IndexName{i_Index},'_raw.mat'];
        FCP_data = importdata(datapath);
        
        Xsource = FCP_data(inds,:);
        d = [Xtarget;Xsource];
      
        sp = [savepath,'/',IndexName{i_Index},'_SWU_beijing_cambridge']; 
        cb(d,batch,mod,sp);
    end
end


clear;
clc;
info = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/CoRR/SubInfo/SubInfo_420.mat');
site = info.Site;
age_corr = info.Age;
sex = info.Sex;
sex(sex==-1)=0;
sex_corr=sex;
motion1 = info.Motion(:,1);
motion2 = info.Motion(:,2);
subid = info.SubID;
SiteSWU_ind = find(site==8);
ResultsSet = {'Results','S2_Results'};
file = {'Results','S2_Results'};
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};

info = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/FCP_Organized/SubInfo/600_subinfo.mat');  % 1-female 2-male
site = info.SiteID;
age = double(info.Age);
sex = info.Sex-1; % 0-female 1-male
motion = info.MeanFDJ;
eye = info.eyeclosed;
sitename = info.SiteName;
subid = info.Subid;

BEIJING_female_inds = find(strcmp(sitename,"Beijing")& sex==0);
Cambridge_male_inds = find(strcmp(sitename,"Cambridge")& sex==1);
refer_site = {BEIJING_female_inds,Cambridge_male_inds};
inds = cell2mat(refer_site');
refer_sitename = sitename(inds);
batch = repelem([1,2,3]',[221,106,74]);
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
mod = [[sex_corr(SiteSWU_ind);sex(inds)],[age_corr(SiteSWU_ind);double(age(inds))]];
% overlap mask
Cov = [sex(inds),double(age(inds)),motion(inds),ones(180,1)];
Datadir = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Revision/Data/SexinSeparateSites';
subid = info.Subid(inds);
methods = {'SWU_beijing_cambridge_para_adj_combat',...
    'SWU_beijing_cambridge_nonpara_adj_combat',...
    'SWU_beijing_cambridge_adjust_covbat',...
    'SWU_beijing_cambridge_para_unadj_combat',... 
    'SWU_beijing_cambridge_nonpara_unadj_combat',...
    'SWU_beijing_cambridge_unadj_covbat'
    };

for ses =1:2
    savepath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Revision/Data/SexinSeparateSites/',file{ses}];
    mkdir(savepath);
    for i_method = 1:numel(methods)
        for i_Index =  2:2% length(IndexName)

            if strcmp(methods{i_method},'raw')
                datapath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',ResultsSet{ses},'/',IndexName{i_Index},'_',methods{i_method},'.mat'];
            else
                datapath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Revision/Data/SexinSeparateSites/',ResultsSet{ses},'/',IndexName{i_Index},'_',methods{i_method},'.mat'];
            end
            data = importdata(datapath);
            if ~strcmp(methods{i_method},'raw')
                d = data(222:end,:);
            end
            %load([savepath,'/',IndexName{i_Index},'_para_adj_combat_SWU_beijing_cambridge.mat']);
            %load([savepath,'/',IndexName{i_Index},'_adjust_covbat_SWU_beijing_cambridge.mat']);
            %covbat = covbat(222:end,:);
            NII_outputpath = fullfile(Datadir,'/Niis/FCP/',ResultsSet{ses},IndexName{i_Index},methods{i_method});
            mkdir(NII_outputpath);
            StatOutDir =  fullfile(Datadir,'/Stat/SWU_beijing_cambridge/FCP/',ResultsSet{ses},IndexName{i_Index},methods{i_method});
            mkdir(StatOutDir);
            
            outputname = [StatOutDir,'/MaleVSFemaleT.nii'];
            GRFOutDir=[StatOutDir,'/GRF/ClusterCorrected001_05/'];
            mkdir(GRFOutDir);
            GRFOutName=[GRFOutDir,'/MaleVsFemaleT'];
            
            
            RawNiiPath = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/FCP_Organized/AllMaps/';
            MaskFile = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Mask/Overlap38810.nii';
            
            %if ~exist(outputname,'file')
            % write nii
            StatList = writeSiteNii(d,IndexName{i_Index},refer_sitename,subid,NII_outputpath,MaskFile,RawNiiPath);
            
            % TWO-SAMPLE-TEST
            y_GroupAnalysis_Image(StatList,Cov,outputname,MaskFile,[],[1,0,0,0],'T',0);
            %end
            % GRF correction
            [Data_corrected,~,Header]= y_GRF_Threshold(outputname,0.001,1,0.05,GRFOutName,MaskFile);

            Data_corrected =Data_corrected(Data_corrected~=0);
            NMax = min(Data_corrected(Data_corrected<0));
            NMin = max(Data_corrected(Data_corrected<0));
            PMin = min(Data_corrected(Data_corrected>0));
            PMax = max(Data_corrected(Data_corrected>0));
            ConnectivityCriterion = 18;
            ClusterSize = 0;
            UnderlayFileName='/mnt/Data3/RfMRILab/Wangyw/software/DPABI_V6.0_ForCamp/Templates/ch2.nii';
            GRFfile = [GRFOutDir,'/ClusterThresholded_MaleVsFemaleT.nii'];
            ColorMap = y_AFNI_ColorMap(12);
            [BrainNetViewerPath, fileN, extn] = fileparts(which('BrainNet.m'));
            SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];
            viewtype='MediumView';
            pics_path = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Revision/Data/SexinSeparateSites/Pics/',ResultsSet{ses},'/',IndexName{i_Index}];
            
            mkdir(pics_path);
            H_BrainNet = y_CallBrainNetViewer(GRFfile,NMin,PMin,ClusterSize,ConnectivityCriterion,SurfFileName,viewtype,ColorMap,NMax,PMax );
            JPGFile=[pics_path,'/',methods{i_method},'.jpeg'];            
            eval(['print -r300 -djpeg -noui ''',JPGFile,''';']);
            
        end
    end
end

function cb(data,batch,mod,savepath)

para_adj_combat = combat(data',batch,mod,1)';
save([savepath,'_para_adj_combat.mat'],'para_adj_combat');
para_unadj_combat = combat(data',batch,[],1)';
save([savepath,'_para_unadj_combat.mat'],'para_unadj_combat');
nonpara_adj_combat = combat(data',batch,mod,0)';
save([savepath,'_nonpara_adj_combat.mat'],'nonpara_adj_combat');
nonpara_unadj_combat = combat(data',batch,[],0)';
save([savepath,'_nonpara_unadj_combat.mat'],'nonpara_unadj_combat');
end