%% subsamplingMMD for CORR dataset
% different reference site
% calculate rmse and ICC

% pooled data seperation by site index
clear;
clc;
%parpool(5);
info = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/SubInfo/SubInfo_420.mat');
site = info.Site;
age = info.Age;
sex = info.Sex;
sex(sex==-1)=0;
motion1 = info.Motion(:,1);
motion2 = info.Motion(:,2);
subid = info.SubID;
us = unique(site);
unique_site = unique(site);

%%
for i = 1: length(us)
    index{i} = find(site == us(i,1));
end
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
%%
%file = {'ResultsS','S2_ResultsS'}; %voxle
file = {'Results','S2_Results'}; % fc
%parpool(5);
for i_index = 1:1%length(IndexName)
    for ses = 1:2
         fprintf('%s , %s is starting......',IndexName{i_index},file{ses});
        cd(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',file{ses}]);
        datapath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',file{ses},'/',IndexName{i_index},'_raw.mat'];
        data = importdata(datapath);
        eval(['motion = motion' num2str(ses) ';']);
        for i = 1: length(us)
            site_data{i} = data(index{i},:);
            site_age{i} = age(index{i},:);
            site_sex{i} = sex(index{i},:);
            site_subid{i} = subid(index{i},:);
            site_motion{i} = motion(index{i},:);
        end
        
        %% define target and source
        %MatrixTran = cell(1,length(unique_site)-1);
        new_site_data = cell(1,length(unique_site));
        new_site_data{1} = site_data{1};
        new_site_age = cell(1,length(unique_site));
        new_site_age{1} = site_age{1};
        new_site_sex = cell(1,length(unique_site));
        new_site_sex{1} = site_sex{1};
        new_site_motion = cell(1,length(unique_site));
        new_site_motion{1} = site_motion{1};
        new_site_subid = cell(1,length(unique_site));
        new_site_subid{1} = site_subid{1};
        
        Xtarget = site_data{1};
        AgeT  = site_age{1};
        SexTarget  = site_sex{1};
        parfor i =2:length(us)
            AgeTarget = AgeT;
            AgeSource =  site_age{i};
            SexSource = site_sex{i};
            Xsource = site_data{i};            
            MotionSource = site_motion{i};
            SubidSource = site_subid{i};
            
            AgeSource_final = AgeSource;
            SexSource_final = SexSource;
            XS = Xsource;
            MotionSource_final = MotionSource;
            SubidSource_final = SubidSource;
            XT = Xtarget;
            % subgroup division critirier
%             AgeSource = floor(double(AgeSource - 5)./double(10));
%             AgeTarget = floor(double(AgeTarget - 5)./double(10));
            AgeSource = double(AgeSource<22)+1;
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
                    fprintf('there are %d sub expeled for source %d \n', size(source_expel,1), us(i) );
                    for j = 1:length(source_expel)
                        XS(find(IndexSource(:,1)==source_expel(j,1)),:)=[];
                        AgeSource_final(find(IndexSource(:,1)==source_expel(j,1)))=[];
                        SexSource_final(find(IndexSource(:,1)==source_expel(j,1)))=[];
                        MotionSource_final(find(IndexSource(:,1)==source_expel(j,1)))=[];
                        SubidSource_final(find(IndexSource(:,1)==source_expel(j,1)))=[];
                        IndexSource(find(IndexSource(:,1)==source_expel(j,1)))=[];
                    end
                else
                    fprintf('there are no sub expeled for source %d \n',i);
                end
                if ~isempty(target_expel)
                    for j = 1:length(target_expel)
                        XT(find(IndexTarget(:,1)==target_expel(j,1)),:)=[];
                        IndexTarget(find(IndexTarget==target_expel(j,1)))=[];
                    end
                end
            end
            tic;
            temp = zeros(size(XS));
            for k = 1:size(temp,2)
                %MatrixTran={};
                Source = XS(:,k);
                Target = XT(:,k);
                [slope,intercept,~] = subsamplingMMD(Source,Target,IndexSource,IndexTarget,100);
                %MatrixTran{i-1}(k,:) = [slope,intercept];
                temp(:,k) = slope*Source+ones(size(Source))*intercept;
                %Here, pvalue is for subsamples from first iteration
            end
            new_site_data{i} = temp;
%             new_site_age{i} = AgeSource_final;
%             new_site_sex{i} = SexSource_final;
%             new_site_motion{i} = MotionSource_final;
%             new_site_subid{i} = SubidSource_final;
            fprintf('site %d fitting completed \n', us(i));
            toc;
        end
        SMA = cell2mat(new_site_data');
        save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',file{ses},'/',IndexName{i_index},'_SMA.mat'],'SMA');
    end
end
