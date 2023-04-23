%% bootstrapping for hypothsis test
%  H0.a  sample size does not affect result when distribution unchanged
%  H0.b  distribution does not affect result when sample size unchanged

%% H0.a sample size does not affect result when distribution unchanged
%   target site: beijing
%   source sites： leidein2200 satitlouis
%   bootstrap settings：
%                       sample size conditions:
%                       1. 7, 4, 4, 3
%                       2. 27, 16, 15, 12
%                       3. 49, 30, 29, 22
%                       bootstrapping times: 50

clear;clc;
addpath(genpath('/mnt/Data3/RfMRILab/Wangyw/software/PNAS2018-master'));
sizes = [7, 4, 4, 3 ;27, 16, 15, 12; 49, 30, 29, 22];
bootTime = 50;
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};

%parpool(20);
for i_Metric = 1: numel(IndexName)
    load(['/mnt/Data6/RfMRILab/Wangyuwei/bootstrapping/Beijing/',IndexName{i_Metric},'.mat']);
    BeijingSubgroup = subgroup;
    
    leiden_2200 = load(['/mnt/Data6/RfMRILab/Wangyuwei/bootstrapping/Leiden_2200/',IndexName{i_Metric},'.mat']);
    sl =  load(['/mnt/Data6/RfMRILab/Wangyuwei/bootstrapping/SaintLouis/',IndexName{i_Metric},'.mat']);
    
    for i_size = 1:size(sizes,1)
        for i_boot  = 1:bootTime
            fprintf('%d / 50 ...',i_boot);
            
            tic;
            outputdir = ['/mnt/Data6/RfMRILab/Wangyuwei/bootstrapping/ExperimentResults/TargetSiteSampleSize/Size',int2str(i_size),'/',int2str(i_boot)];
            mkdir(outputdir);
            [SampleInds,targetlabel] = subgroup_subsamping(BeijingSubgroup,sizes(i_size,:),0);
            TargetData = sdata(SampleInds,:);
            
            sourcedata = {leiden_2200.sdata;sl.sdata};
            sourcelabel = {leiden_2200.subgroup;sl.subgroup};

            sma = [];
            for i_source =  1:numel(sourcedata)
                SD = sourcedata{i_source};
                SL = sourcelabel{i_source};
                temp = zeros(size(SD));
                parfor k = 1:size(temp,2)
                    source = SD(:,k);
                    target = TargetData(:,k);
                    [slope,intercept,pvalue] = subsamplingMMD(source,target,SL,targetlabel,50);
                    temp(:,k) = slope*source+ones(size(source))*intercept;
                    if sum(isnan(temp(:,k)))
                        warning('there are nans');
                    end
                end
                sma{i_source} = temp;
            end
            %sma = cell2mat(sma);
            save([outputdir,'/',IndexName{i_Metric},'_SMA.mat'],'sma','SampleInds'); 
            toc;
        end
    end
end
%% H0.b  distribution does not affect result when sample size unchanged
%   target site: beijing
%   source sites： leidein2200 satitlouis
%   bootstrapp settings：
%                       sample distribution conditions:
%                       1. left bias 30, 5, 30 , 5
%                       2. right bias 5, 30, 5, 30
%                       bootstrapping times: 50

%%
function [Subsample_ind,Labels] = subgroup_subsamping(groups,subsize,isrepeat)
% 1 get indexs for each sub group
unique_group_labels = unique(groups);
for i = 1:numel(unique_group_labels)
    subgroup_index{i} = find(groups==unique_group_labels(i));
end

% 2 get assigned size of samples for each group
if numel(subsize)~=numel(unique_group_labels)
    error('The number of sizes of subgroups does not match the number of unique groups.');
elseif ~isrepeat
    for i = 1:numel(unique_group_labels)
        Subgroup_ind_choice{i} = subgroup_index{i}(randperm(numel(subgroup_index{i}),subsize(i)));
    end
else
    for i = 1:numel(unique_group_labels)
        Subgroup_ind_choice{i} = subgroup_index{i}(randi(numel(subgroup_index{i}),subsize(i)));
    end
end
Subsample_ind=cell2mat(Subgroup_ind_choice');
Labels = repelem(unique_group_labels,subsize);
end



