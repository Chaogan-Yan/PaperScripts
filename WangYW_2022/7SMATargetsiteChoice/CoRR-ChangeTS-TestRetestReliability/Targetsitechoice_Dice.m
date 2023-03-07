clear;clc;
load /mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/SubInfo/SubInfo_420.mat;
SiteUnique=unique(Site);
SiteUnique=setdiff(SiteUnique,7);
ResultsSet = {'Results','S2_Results'};
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
MeasurePrefixSet={'szReHoMap_','szALFFMap_','szfALFFMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_'};
status  = {'without_ScannerRegressor','with_ScannerRegressor'};
DataDir1 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/TargetSiteChoice/Results';
DataDir2 = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Stats/TargetSiteChoice/S2_Results';

for i_Index = 1:length(IndexName)
    Diceall =zeros(length(SiteUnique));
    Dice_between = zeros(length(SiteUnique));
    for i_site =1:length(SiteUnique)
        for j_site =1:length(SiteUnique)
            
            if ~strcmp(IndexName{i_Index},'FC_D142')
                if i_site == j_site
                    File=[DataDir1,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                    Data1=y_Read(File);
                    File=[DataDir2,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                    Data2= y_Read(File);
                    Data=Data1+Data2;
                    %y_Write(FDR_mask,['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/statsResults/Stat_again/VOXEL/402/DiceMap/',IndexName{i_Index
                    
                    %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_FDR_mask'],'FDR_mask');
                    %y_Write(FDR_mask,['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/statsResults/Stat_again/VOXEL/402/DiceMap/',IndexName{i_Index},'_FDR_mask'])
                    %num{i_method}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
                    %                         num{i_Index}= length(find(Data==2));
                    %                         Reproducibility{i_Index}=mean(Data(Data>0));
                    %                         Jaccard{i_Index}=length(find(Data==2))/length(find(Data>=1));
                    Dice_between(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                elseif i_site < j_site
                    File=[DataDir1,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                    Data1=y_Read(File);
                    File=[DataDir1,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                    Data2= y_Read(File);
                    Data=Data1+Data2;
                    Diceall(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                    
                    File=[DataDir1,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                    Data1=y_Read(File);
                    File=[DataDir2,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                    Data2= y_Read(File);
                    Data=Data1+Data2;
                    Dice_between(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                elseif i_site > j_site
                    File=[DataDir2,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                    Data1=y_Read(File);
                    File=[DataDir2,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                    Data2= y_Read(File);
                    Data=Data1+Data2;
                    Diceall(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                    
                    File=[DataDir2,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                    Data1=y_Read(File);
                    File=[DataDir1,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/FDR/MaleVsFemaleT.nii'];
                    Data2= y_Read(File);
                    Data=Data1+Data2;
                    Dice_between(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                end
                
            else
                if i_site == j_site
                    File=[DataDir1,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/MaleVsFemaleT_FDRbinarized.mat'];
                    Data1=importdata(File);
                    File=[DataDir2,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/MaleVsFemaleT_FDRbinarized.mat'];
                    Data2= importdata(File);
                    Data=Data1+Data2;
                    Dice_between(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                    
                elseif i_site < j_site
                    File=[DataDir1,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/MaleVsFemaleT_FDRbinarized.mat'];
                    Data1=importdata(File);
                    File=[DataDir1,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/MaleVsFemaleT_FDRbinarized.mat'];
                    Data2= importdata(File);
                    Data=Data1+Data2;
                    Diceall(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                    
                    File=[DataDir1,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/MaleVsFemaleT_FDRbinarized.mat'];
                    Data1=importdata(File);
                    File=[DataDir2,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/MaleVsFemaleT_FDRbinarized.mat'];
                    Data2= importdata(File);
                    Data=Data1+Data2;
                    Dice_between(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                elseif i_site > j_site
                    File=[DataDir2,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/MaleVsFemaleT_FDRbinarized.mat'];
                    Data1=importdata(File);
                    File=[DataDir2,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'//without_ScannerRegressor/MaleVsFemaleT_FDRbinarized.mat'];
                    Data2= importdata(File);
                    Data=Data1+Data2;
                    Diceall(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                    
                    File=[DataDir2,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/MaleVsFemaleT_FDRbinarized.mat'];
                    Data1=importdata(File);
                    File=[DataDir1,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/MaleVsFemaleT_FDRbinarized.mat'];
                    Data2= importdata(File);
                    Data=Data1+Data2;
                    Dice_between(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                end
            end
        end
        %FDR_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
        
    end
    csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/TargetSiteChoice/Analysis/Dice/withinSession_Dice_FDR',IndexName{i_Index},'.csv'],Diceall);
    csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/TargetSiteChoice/Analysis/Dice/betweenSession_Dice_FDR',IndexName{i_Index},'.csv'],Dice_between);
end

for i_Index = 1:length(IndexName)
    Diceall =zeros(length(SiteUnique));
    Dice_between = zeros(length(SiteUnique));
    for i_site =1:length(SiteUnique)
        for j_site =1:length(SiteUnique)
            if ~strcmp(IndexName{i_Index},'FC_D142')
                if i_site == j_site
                    File=[DataDir1,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data1=y_Read(File);
                    File=[DataDir2,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data2= y_Read(File);
                    Data=Data1+Data2;
                    %y_Write(FDR_mask,['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/statsResults/Stat_again/VOXEL/402/DiceMap/',IndexName{i_Index
                    
                    %save(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/OtherAnalysis/Dice/Replicability/',IndexName{i_Index},'_FDR_mask'],'FDR_mask');
                    %y_Write(FDR_mask,['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/statsResults/Stat_again/VOXEL/402/DiceMap/',IndexName{i_Index},'_FDR_mask'])
                    %num{i_method}= [length(find(Data==2)),length(find(Data==1)),length(find(Data>0))];
                    %                         num{i_Index}= length(find(Data==2));
                    %                         Reproducibility{i_Index}=mean(Data(Data>0));
                    %                         Jaccard{i_Index}=length(find(Data==2))/length(find(Data>=1));
                    Dice_between(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                elseif i_site < j_site
                    File=[DataDir1,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data1=y_Read(File);
                    File=[DataDir1,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data2= y_Read(File);
                    Data=Data1+Data2;
                    Diceall(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                    
                    File=[DataDir1,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data1=y_Read(File);
                    File=[DataDir2,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data2= y_Read(File);
                    Data=Data1+Data2;
                    Dice_between(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                elseif i_site > j_site
                    File=[DataDir2,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data1=y_Read(File);
                    File=[DataDir2,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data2= y_Read(File);
                    Data=Data1+Data2;
                    Diceall(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                    
                    File=[DataDir2,'/Site',num2str(SiteUnique(i_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data1=y_Read(File);
                    File=[DataDir1,'/Site',num2str(SiteUnique(j_site)),'_asT/',IndexName{i_Index},'/without_ScannerRegressor/SignificantBinarized/GRF/ClusterCorrected001_05/MaleVsFemaleT.nii'];
                    Data2= y_Read(File);
                    Data=Data1+Data2;
                    Dice_between(i_site,j_site)= 2*length(find(Data==2))/(length(find(Data>=1))+length(find(Data==2)));
                end
            end
        end
    end
    %FDR_RESULT(:,2*i_Index-1:2*i_Index) = [Dice,num];
    csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/TargetSiteChoice/Analysis/Dice/withinSession_Dice_GRF_',IndexName{i_Index},'.csv'],Diceall);
    csvwrite(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/TargetSiteChoice/Analysis/Dice/betweenSession_Dice_GRF_',IndexName{i_Index},'.csv'],Dice_between);
end