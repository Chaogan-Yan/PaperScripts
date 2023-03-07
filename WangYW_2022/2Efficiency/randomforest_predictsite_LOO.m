%clear;clc;

IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
methods = {'raw','reg','adj','lmm','para_adj_combat','nonpara_adj_combat','nonpara_unadj_combat','para_unadj_combat','adj_covbat','unadj_covbat','SMA','ICVAE'};

%parpool(12);
tic;
for i_Metric = 5:numel(IndexName)
    for i_Method = 11:numel(methods) %notice 
        Datapath = [IndexName{i_Metric},'_',methods{i_Method},'.mat'];
        alldata=importdata(Datapath);
        if isstruct(alldata)
            alldata = struct2array(alldata);
        end
        %consistency=[];
        
        for i_perm = 1:123
            test_ind = i_perm;
            train_ind = setdiff(1:123,test_ind)';
            
            class_labels = repelem({'IPCASGE','PKUGE','PKUSIEMENS'},41);
            class_labels = char(class_labels);
            % to permute data first
            random_order = randi(122,122,1);
            
            X = alldata(train_ind(random_order,:),:);
            Y = class_labels(train_ind(random_order,1),:);
            
            
            % TREEBAGGER TRY
            Ntrees=600;
            B = TreeBagger(Ntrees,alldata(train_ind(random_order,:),:),class_labels(train_ind(random_order,1),:),...
                 "Method","classification",...
                 "NumPredictorstoSample","all",...
                 "OOBPrediction", "off");
           
            predChar = B.predict(alldata(train_ind,:));
            t = double(cellfun(@strcmp,cellstr(class_labels(train_ind,:)),predChar1));
            training_consistency(i_perm,1) = sum(t)/length(t);
            
            predChar = B.predict(alldata(test_ind,:));
            predchar_list(i_perm,1) =  predChar ;
            t = double(cellfun(@strcmp,cellstr(class_labels(test_ind,:)),predChar));
            testing_consistency(i_perm,1) = t;
          
        end
        accuracy_all{i_Metric,i_Method} = {training_consistency,testing_consistency}; 
        accuracy{i_Metric,i_Method} = cellfun(@mean,accuracy_all{i_Metric,i_Method});
        preadchar{i_Metric,i_Method} = predchar_list;
    end
    save('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Revision/Data/CovBatCovariance/SitePrediction/TSP3_LOO_subject.mat','accuracy','accuracy_all','preadchar');
end
