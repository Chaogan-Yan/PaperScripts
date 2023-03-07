%% fitMMD for TST dataset
% different reference site
% calculate rmse and ICC
% IndexName ={'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
% % MeasurePrefixSet={'szReHoMap_','szALFFMap_','szfALFFMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_'};
% datapath = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/TSP3/ResultsS';
% refer_site = 'pku_ge';
% 
% outputpath = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/TSP3/ResultsS';
% parpool(5);
% parfor i = 1:numel(IndexName)
%     dataname = [datapath,'/',IndexName{i},'_raw.mat'];
%     outputname = [outputpath,'/',IndexName{i},'_SMA.mat'];
%     fitMMD4TST(dataname,refer_site,IndexName{i},outputname);
% end

function fitMMD4TST(datapath,refer_site,feature_name,outputpath)

    %datapath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/TST/TSTafterHarmonization/3SIte41Sub/ResultsS/',IndexName{i},'_raw.mat'];
    data = importdata(datapath);

    ipcas_ge_index = [1:41];
    pku_ge_index = [42:82];
    pku_siemens_index = [83:123];

    if size(data,1) ~= 123
        data= data';
    end    
    ipcas_ge = data(ipcas_ge_index,:);
    pku_ge = data(pku_ge_index,:);
    pku_siemens = data(pku_siemens_index,:);

    site_names = {'ipcas_ge','pku_ge','pku_siemens'};
    transform_mat =zeros(size(data,2),length(site_names)-1);
   
    %for target_index =  1:length(site_names)
    target_site = refer_site; %pku_ge
    eval(['target_data  = ',target_site,';']);
    
    source_names = site_names(find(~ismember(site_names,target_site)));
    trans_data = zeros(size(data));
    
    transformation.b = transform_mat;
    transformation.w =  transform_mat;
    for source_index = 1:length(source_names)
        source_site = source_names{source_index};
        eval(['source_data = ' source_site ';']);
        fprintf('we are map %s to %s', source_site, target_site)
        for ele = 1:size(data,2)
            [slope,intercept,~] = fitMMD(source_data(:,ele),target_data(:,ele),0);
            transformation.b(ele,source_index) = intercept;
            transformation.w(ele,source_index) = slope;
            temp(:,ele) = source_data(:,ele).*slope+intercept;
        end
        eval(['trans_data(' source_site '_index,:) = temp;']);
    end
    eval(['trans_data(' target_site '_index,:) = target_data']);
    save([outputpath,'/',feature_name,'_MMD_',target_site,'.mat'],'trans_data');
    %end
end