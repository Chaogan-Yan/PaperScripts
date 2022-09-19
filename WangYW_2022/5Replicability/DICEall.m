%% call dice
clear;clc;
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
MeasurePrefixSet={'szReHoMap_','szALFFMap_','szfALFFMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_'};
Harmo_methods = {'raw','para_adj_combat','nonpara_adj_combat','para_unadj_combat','nonpara_unadj_combat','SMA'};
status  = {'without_covariates','with_covariates'};
caldice_FDR =1;
caldice_GRF =1;
    DataDir1 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/Results';
    DataDir2 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/S2_Results';
    DataDir3 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/FCP/MalevsFemale/Results/SignificantBinarized';
    DataDir4 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/FCP/MalevsFemale/S2_Results/SignificantBinarized';
%% 123
for i_Index = 1:length(IndexName)

    if caldice_FDR == 1
        
        Reproducibility=cell(length(Harmo_methods),1);
        Jaccard= cell(length(Harmo_methods),1);
        Dice= cell(length(Harmo_methods),1);
        num= cell(length(Harmo_methods),1);
        if ~strcmp(IndexName{i_Index},'FC_D142')
          
            for i_method = 1:length(Harmo_methods)
                File=[DataDir1,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT.nii'];
                Data1=y_Read(File);
                File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT.nii'];
                Data2= y_Read(File);
                Data12=(Data1+Data2==2);
                File=[DataDir3,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/FDR/MaleVsFemaleT.nii'];
                Data3 = y_Read(File);
                Data = Data12 +Data3;
                %y_Write(FDR_mask,['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/statsResults/Stat_again/VOXEL/402/DiceMap/',IndexName{i_Index
                
                %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_FDR_mask'],'FDR_mask');
                %y_Write(FDR_mask,['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/statsResults/Stat_again/VOXEL/402/DiceMap/',IndexName{i_Index},'_FDR_mask'])
                %num{i_method}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
                num{i_method}= length(find(Data==2));
                Reproducibility{i_method}=mean(Data(Data>0));
                Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
                Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
            end
        else
            for i_method = 1:length(Harmo_methods)
                File=[DataDir1,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT_FDRBinarized.mat'];
                Data1=importdata(File);
                File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT_FDRBinarized.mat'];
                Data2= importdata(File);
                Data12=(Data1+Data2==2);
                File=[DataDir3,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/FDR/MaleVsFemale_TFDRBinarized.mat'];
                Data3 = importdata(File);
                Data = Data12 +Data3;

                num{i_method}= length(find(Data==2));
                Reproducibility{i_method}=mean(Data(find(Data>0)));
                Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
                Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
            end
        end
        FDR_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
        csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/FDR_DICE31.csv'],FDR_RESULT);
        %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_FDR_DICE'],'FDR_RESULT');
    end
    
    if caldice_GRF == 1 && ~strcmp(IndexName{i_Index},'FC_D142')
        Reproducibility=cell(length(Harmo_methods),1);
        Jaccard= cell(length(Harmo_methods),1);
        Dice= cell(length(Harmo_methods),1);
        num= cell(length(Harmo_methods),1);
        
        for i_method = 1:length(Harmo_methods)
            
            File=[DataDir1,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data1=y_Read(File);
            File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data2= y_Read(File);
            Data12=(Data1+Data2==2);
            File=[DataDir3,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data3 = y_Read(File);
            Data = Data12 +Data3;
            
            %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_GRF_mask'],'GRF_mask');
            %num{i_method}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
            num{i_method}= length(find(Data==2));
            Reproducibility{i_method}=mean(Data(find(Data>0)));
            Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
            Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
        end
        GRF_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
        csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/GRF_DICE31.csv'],GRF_RESULT);
    end
end

%% 124
for i_Index = 1:length(IndexName)
    DataDir1 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/Results';
    DataDir2 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/S2_Results';
    DataDir3 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/FCP/MalevsFemale/Results/SignificantBinarized';
    DataDir4 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/FCP/MalevsFemale/S2_Results/SignificantBinarized';
    if caldice_FDR == 1
        
        Reproducibility=cell(length(Harmo_methods),1);
        Jaccard= cell(length(Harmo_methods),1);
        Dice= cell(length(Harmo_methods),1);
        num= cell(length(Harmo_methods),1);
        if ~strcmp(IndexName{i_Index},'FC_D142')       
            for i_method = 1:length(Harmo_methods)
                File=[DataDir1,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT.nii'];
                Data1=y_Read(File);
                File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT.nii'];
                Data2= y_Read(File);
                Data12=(Data1+Data2==2);
                File=[DataDir4,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/FDR/MaleVsFemaleT.nii'];
                Data3 = y_Read(File);
                Data = Data12 +Data3;
                %y_Write(FDR_mask,['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/statsResults/Stat_again/VOXEL/402/DiceMap/',IndexName{i_Index
                
                %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_FDR_mask'],'FDR_mask');
                %y_Write(FDR_mask,['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/statsResults/Stat_again/VOXEL/402/DiceMap/',IndexName{i_Index},'_FDR_mask'])
                %num{i_method}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
                num{i_method}= length(find(Data==2));
                Reproducibility{i_method}=mean(Data(Data>0));
                Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
                Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
            end
        else
            for i_method = 1:length(Harmo_methods)
                File=[DataDir1,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT_FDRBinarized.mat'];
                Data1=importdata(File);
                File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT_FDRBinarized.mat'];
                Data2= importdata(File);
                Data12=(Data1+Data2==2);
                File=[DataDir4,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/FDR/MaleVsFemaleT_FDRBinarized.mat'];
                Data3 = importdata(File);
                Data = Data12 +Data3;
                
                
                %num{i_method}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
                num{i_method}= length(find(Data==2));
                Reproducibility{i_method}=mean(Data(find(Data>0)));
                Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
                Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
            end
        end
        FDR_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
        csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/FDR_DICE32.csv'],FDR_RESULT);
        %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_FDR_DICE'],'FDR_RESULT');
    end
    
    if caldice_GRF == 1 && ~strcmp(IndexName{i_Index},'FC_D142')
        Reproducibility=cell(length(Harmo_methods),1);
        Jaccard= cell(length(Harmo_methods),1);
        Dice= cell(length(Harmo_methods),1);
        num= cell(length(Harmo_methods),1);
        
        for i_method = 1:length(Harmo_methods)
            
            File=[DataDir1,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data1=y_Read(File);
            File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data2= y_Read(File);
            Data12=(Data1+Data2==2);
            File=[DataDir4,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data3 = y_Read(File);
            Data = Data12 +Data3;
            
            %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_GRF_mask'],'GRF_mask');
            %num{i_method}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
            num{i_method}= length(find(Data==2));
            Reproducibility{i_method}=mean(Data(find(Data>0)));
            Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
            Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
        end
        GRF_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
        csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/GRF_DICE32.csv'],GRF_RESULT);
    end
end
%% 12
for i_Index = 1:length(IndexName)
    DataDir1 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/Results';
    DataDir2 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/S2_Results';
    DataDir3 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/FCP/MalevsFemale/Results/SignificantBinarized';
    DataDir4 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/FCP/MalevsFemale/S2_Results/SignificantBinarized';
    if caldice_FDR == 1
        
        Reproducibility=cell(length(Harmo_methods),1);
        Jaccard= cell(length(Harmo_methods),1);
        Dice= cell(length(Harmo_methods),1);
        num= cell(length(Harmo_methods),1);
        if ~strcmp(IndexName{i_Index},'FC_D142')       
            for i_method = 1:length(Harmo_methods)
                File=[DataDir1,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT.nii'];
                Data1=y_Read(File);
                File=[DataDir3,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/FDR/MaleVsFemaleT.nii'];
                Data3 = y_Read(File);
                Data = Data1 +Data3;
                num{i_method}= length(find(Data==2));
                Reproducibility{i_method}=mean(Data(Data>0));
                Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
                Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
            end
        else
            for i_method = 1:length(Harmo_methods)
                File=[DataDir1,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT_FDRBinarized.mat'];
                Data1=importdata(File);
                File=[DataDir3,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/FDR/MaleVsFemaleT_FDRBinarized.mat'];
                Data3 = importdata(File);
                Data = Data1 +Data3;

                num{i_method}= length(find(Data==2));
                Reproducibility{i_method}=mean(Data(find(Data>0)));
                Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
                Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
            end
        end
        FDR_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
        csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/FDR_DICE12.csv'],FDR_RESULT);
        %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_FDR_DICE'],'FDR_RESULT');
    end
    
    if caldice_GRF == 1 && ~strcmp(IndexName{i_Index},'FC_D142')
        Reproducibility=cell(length(Harmo_methods),1);
        Jaccard= cell(length(Harmo_methods),1);
        Dice= cell(length(Harmo_methods),1);
        num= cell(length(Harmo_methods),1);
        
        for i_method = 1:length(Harmo_methods)
            
            File=[DataDir1,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data1=y_Read(File);
            File=[DataDir3,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data3 = y_Read(File);
            Data = Data1 +Data3;
            
            %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_GRF_mask'],'GRF_mask');
            %num{i_method}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
            num{i_method}= length(find(Data==2));
            Reproducibility{i_method}=mean(Data(find(Data>0)));
            Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
            Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
        end
        GRF_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
        csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/GRF_DICE12.csv'],GRF_RESULT);
    end
end
%% 21
for i_Index = 1:length(IndexName)
    DataDir1 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/Results';
    DataDir2 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/CoRR/MalevsFemale/SignificantBinarized/S2_Results';
    DataDir3 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/FCP/MalevsFemale/Results/SignificantBinarized';
    DataDir4 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/FCP/MalevsFemale/S2_Results/SignificantBinarized';
    if caldice_FDR == 1
        
        Reproducibility=cell(length(Harmo_methods),1);
        Jaccard= cell(length(Harmo_methods),1);
        Dice= cell(length(Harmo_methods),1);
        num= cell(length(Harmo_methods),1);
        if ~strcmp(IndexName{i_Index},'FC_D142')       
            for i_method = 1:length(Harmo_methods)
                File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT.nii'];
                Data1=y_Read(File);
                File=[DataDir4,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/FDR/MaleVsFemaleT.nii'];
                Data3 = y_Read(File);
                Data = Data1 +Data3;
                num{i_method}= length(find(Data==2));
                Reproducibility{i_method}=mean(Data(Data>0));
                Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
                Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
            end
        else
            for i_method = 1:length(Harmo_methods)
                File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/FDR/MaleVsFemaleT_FDRBinarized.mat'];
                Data1=importdata(File);
                File=[DataDir4,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/FDR/MaleVsFemaleT_FDRBinarized.mat'];
                Data3 = importdata(File);
                Data = Data1 +Data3;

                num{i_method}= length(find(Data==2));
                Reproducibility{i_method}=mean(Data(find(Data>0)));
                Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
                Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
            end
        end
        FDR_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
        csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/FDR_DICE21.csv'],FDR_RESULT);
        %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_FDR_DICE'],'FDR_RESULT');
    end
    
    if caldice_GRF == 1 && ~strcmp(IndexName{i_Index},'FC_D142')
        Reproducibility=cell(length(Harmo_methods),1);
        Jaccard= cell(length(Harmo_methods),1);
        Dice= cell(length(Harmo_methods),1);
        num= cell(length(Harmo_methods),1);
        
        for i_method = 1:length(Harmo_methods)
            
            File=[DataDir2,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_ScannerRegressor/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data1=y_Read(File);
            File=[DataDir4,'/',Harmo_methods{i_method},'/',IndexName{i_Index},'/without_covariates/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
            Data3 = y_Read(File);
            Data = Data1 +Data3;
            
            %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_GRF_mask'],'GRF_mask');
            %num{i_method}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
            num{i_method}= length(find(Data==2));
            Reproducibility{i_method}=mean(Data(find(Data>0)));
            Jaccard{i_method}=length(find(Data==2))/length(find(Data>=1));
            Dice{i_method}=2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
        end
        GRF_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
        csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/GRF_DICE21.csv'],GRF_RESULT);
    end
end