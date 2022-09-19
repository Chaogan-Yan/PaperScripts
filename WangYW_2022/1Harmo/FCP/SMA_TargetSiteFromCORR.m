% pooled data seperation by site index
clear;
clc;
% for that the group division problem , we need to delete site 4 and 11
% from our site list
corrinfo = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/CoRR/SubInfo/SubInfo_420.mat');
corrsite = corrinfo.Site;
corrage = corrinfo.Age;
corrsex = corrinfo.Sex;
corrsex(corrsex==-1)=0; % 0 - female 1-male 
corrsubid = corrinfo.SubID;
corrus = unique(corrsite);

Ind = find(corrsite == 1);
target_site_age = corrage(Ind,:);
target_site_sex =  corrsex(Ind,:);
clear corr*;
%%
info = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/old/FCP_Organized/SubInfo/600_subinfo.mat');  % 1-female 2-male
site = info.SiteID;
age = double(info.Age);
sex = info.Sex-1; % 0-female 1-male
motion = info.MeanFDJ;
eye = info.eyeclosed;
sitename = info.SiteName;
subid = info.Subid;
%IndexName = {'FC'}
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
unique_site = unique(site);
target_site =unique(site);

% overlap mask

%%
for i = 1: length(unique_site)
    index{i} = find(site == unique_site(i,1));
end

%%
parpool(13);
%parpool(4);
Rdir = {'Results','S2_Results'};
for ses = 1:2
    for i_Index =  2: 2%length(IndexName)
        fprintf('%s , %s is starting......',IndexName{i_Index},Rdir{ses});
        % take session1 as target data
        corr_datapath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',Rdir{ses},'/',IndexName{i_Index},'_raw.mat'];
        corr_data = importdata(corr_datapath);
        target_data = corr_data(Ind,:);
        
        datapath =[ '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/FCP/',Rdir{ses},'/',IndexName{i_Index},'_raw.mat'];
        data = importdata(datapath);
        for i = 1: length(unique_site)
            site_data{i} = data(index{i},:);
            site_age{i} = age(index{i},:);
            site_sex{i} = sex(index{i},:);
            site_subid{i} = subid(index{i},:);
            site_motion{i} = motion(index{i},:);
%             site_eye{i} = eye(index{i},:);
        end
        %% define target and source
        % now target site is beyond the dataset, but from another dataset
        %MatrixTran = [];

        %new_site_eye = cell(1,length(unique_site));
        unique_qualified_site = find(unique_site~=3 & unique_site~=4);
        new_site_data = cell(1,length(unique_qualified_site));
        for iu = 1:length(unique_qualified_site)
            new_site_data{iu} = zeros(size(site_data{unique_qualified_site(iu)}));
        end     
        %new_site_data{15} = site_data{1};
        new_site_age = cell(1,length(unique_qualified_site));
        %new_site_age{15} = site_age{1};
        new_site_sex = cell(1,length(unique_qualified_site));
        %new_site_sex{15} = site_sex{1};
        new_site_motion = cell(1,length(unique_qualified_site));
        %new_site_motion{15} = site_motion{1};
        parfor i =1:length(unique_qualified_site)
            addpath(genpath('/mnt/Data3/RfMRILab/Wangyw/software/PNAS2018-master'));
            % if i ~=3 && i ~=4
            tic;
           
            Xtarget = target_data;
            AgeTarget  = target_site_age;
            SexTarget  =target_site_sex;
            AgeSource =  site_age{unique_qualified_site(i)};
            SexSource = site_sex{unique_qualified_site(i)};
            Xsource = site_data{unique_qualified_site(i)};
            
            MotionSource = site_motion{unique_qualified_site(i)};
            SubidSource = site_subid{unique_qualified_site(i)};
            AgeSource_final = AgeSource;
            SexSource_final = SexSource;
            XS = Xsource;
            MotionSource_final = MotionSource;
            SubidSource_final = SubidSource;
            XT = Xtarget;
            % ORIGINAL
            %             AgeLow = max(min(AgeSource), min(AgeTarget));
            %             AgeHigh = min(max(AgeSource), max(AgeTarget));
            %             index_source = find(AgeSource >= AgeLow ); %& AgeSource <= AgeHigh
            %             index_target = find(AgeTarget >= AgeLow); % & AgeTarget <= AgeHigh 
            %
            %             AgeTarget = AgeTarget(index_target);
            %             SexTarget = SexTarget(index_target);
            %             XT = Xtarget(index_target,:);
            %
            %             AgeSource(AgeSource<=23)=0;
            %             AgeSource(AgeSource>23)=1;
            %             AgeTarget(AgeTarget<=23)=0;
            %             AgeTarget(AgeTarget>23)=1;
            %             IndexSource = AgeSource.*2 + double(SexSource) ;
            %             IndexTarget = AgeTarget.*2 + double(SexTarget) ;
%             AgeSource = floor(double(AgeSource - 5)./double(10));
%             AgeTarget = floor(double(AgeTarget - 5)./double(10));
%             IndexSource = AgeSource + double(SexSource) ;
%             IndexTarget = AgeTarget + double(SexTarget) ;
            AgeSource = double(AgeSource<22)+1;  % <22 - 2 >=22 1
            AgeTarget = double(AgeTarget<22)+1;
            IndexSource = AgeSource + 2*double(SexSource) ; 
            %tabulate(IndexSource)   
            IndexTarget = AgeTarget + 2*double(SexTarget) ;
            %tabulate(IndexTarget)
            source_expel = [];
            if ~isequal(unique(IndexSource),unique(IndexTarget))
                l1 = unique(IndexSource);
                l2 = unique(IndexTarget);
                cross = intersect(l1,l2);
                source_expel = setdiff(l1,cross);
                target_expel = setdiff(l2,cross);
                if ~isempty(source_expel)
                    fprintf('there are %d sub expeled for source %d \n', size(source_expel,1), unique_site(unique_qualified_site(i)) );
                    for j = 1:length(source_expel)
                        XS(find(IndexSource(:,1)==source_expel(j,1)),:)=[];
                        AgeSource_final(find(IndexSource(:,1)==source_expel(j,1)))=[];
                        SexSource_final(find(IndexSource(:,1)==source_expel(j,1)))=[];
                        MotionSource_final(find(IndexSource(:,1)==source_expel(j,1)))=[];
                        %EyeSource_final(find(IndexSource(:,1)==source_expel(j,1)))=[];
                        SubidSource_final(find(IndexSource(:,1)==source_expel(j,1)))=[];
                        IndexSource(find(IndexSource(:,1)==source_expel(j,1)))=[];
                    end
                else
                    fprintf('there are no sub expeled for source %d \n',unique_site(unique_qualified_site(i)));
                end
                if ~isempty(target_expel)
                    for j = 1:length(target_expel)
                        XT(find(IndexTarget(:,1)==target_expel(j,1)),:)=[];
                        IndexTarget(find(IndexTarget==target_expel(j,1)))=[];
                    end
                end
            else
                fprintf('there are no sub expeled for source %d \n',unique_site(unique_qualified_site(i)));
            end
            
            temp = zeros(size(XS));
            for k = 1:size(temp,2)
                %MatrixTran={};
                str = ['coming soon ',num2str(k/size(temp,2)),'%'];
                Source = XS(:,k);
                Target = XT(:,k);
                [slope,intercept,pvalue] = subsamplingMMD(Source,Target,IndexSource,IndexTarget,100);
                %MatrixTran = [slope,intercept,pvalue];
                temp(:,k) = slope*Source+ones(size(Source))*intercept;
                if sum(isnan(temp(:,k)))
                    warning('there are %d nan in site %d for feature %d \n', sum(isnan(temp(:,k))),unique_site(unique_qualified_site(i)),k);
                end
            end
            new_site_data{i} = temp;
%                         new_site_age{i} = AgeSource_final;
%                         new_site_sex{i} = SexSource_final;
%                         new_site_motion{i} = MotionSource_final;
%                         new_site_subid{i} = SubidSource_final;
                        %new_site_eye{i} = EyeSource_final;
            fprintf('site %d fitting completed \n', unique_site(unique_qualified_site(i)));
            toc;
        end
%             Age = cell2mat(new_site_age');
%             Sex = cell2mat(new_site_sex');
%             Motion = cell2mat(new_site_motion');
%             Subid =cat(1,new_site_subid{:});
%             %Eyestate  = cell2mat(new_site_eye);
%             save('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/FCP/demographic.mat','Age','Sex','Motion','Subid');
        SMA = cell2mat(new_site_data');        
        save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/FCP/',Rdir{ses},'/',IndexName{i_Index},'_SMA.mat'],'SMA');
    end
end


