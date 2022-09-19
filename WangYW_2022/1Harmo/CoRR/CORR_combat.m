
clc;clear;
%% 1 read nii
addpath(fullfile(matlabroot,'toolbox','stats','stats'))

load /mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/SubInfo/SubInfo_420.mat;
    
% mask
MaskFile ='/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Mask/Overlap38810.nii'; 

 if length(unique(Site))~=max(Site)
     s = unique(Site);
     for i = 1:length(s)
         Site(Site==s(i)) =  i ;
     end
 end
 
DataDir='/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/ResultsfromPaper/BeSess_SinRest'
OutDir= '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/';
MeasureSet ={'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','ROISignals_FunImgARCWF'};
MeasurePrefixSet={'szReHoMap_','szALFFMap_','szfALFFMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_','ROICorrelation_FisherZ_'};

Mask = y_Read(MaskFile);
mask_size = size(Mask);
% mapping index
x =  find(reshape(Mask,[],1)==1);
newmap = zeros(size(reshape(Mask,[],1)));
newmap(x,1) = 1; 


%%
SResultsSet={'ResultsS','S2_ResultsS'};
ResultsSet = {'Results','S2_Results'};
%parpool(2);
for isession = 1:numel(ResultsSet)
    for iMeasure = 5:numel(MeasureSet)
        raw = [];
        FileList = [];
        if ~strcmp(MeasureSet{iMeasure},'ROISignals_FunImgARCWF') 
            for iSub=1:length(SubID)
                FileList{iSub,1}=[DataDir,'/',SResultsSet{isession},'/',MeasureSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
                [data,~,~,~] = y_ReadAll(FileList{iSub,1});
                new_data = reshape(data,[],1).* newmap;
                raw(iSub,:) = new_data(new_data~=0) ;
            end
            save([OutDir,'/',ResultsSet{isession},'/',MeasureSet{iMeasure},'_raw.mat'],'raw');
        else
            load('/mnt/Data3/RfMRILab/Wangyw/software/DPABI_V6.0_ForCamp/Templates/Dosenbach_Science_160ROIs_Info.mat')
            for iSub=1:length(SubID)
                if isession==2
                    FileList{iSub,1}=[DataDir,'/',ResultsSet{isession},'/S2_',MeasureSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.mat'];
                else
                    FileList{iSub,1}=[DataDir,'/',ResultsSet{isession},'/',MeasureSet{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.mat'];
                end
                data =importdata(FileList{iSub,1});
                data = data(ROIIndex1409_ExcludeCerebellum_142,ROIIndex1409_ExcludeCerebellum_142);
                trildata = tril(data,-1);
                vec = reshape(trildata(trildata~=0),1,[]);
                raw(iSub,:) = vec;
            end
            save([OutDir,'/',ResultsSet{isession},'/FC_D142_raw.mat'],'raw');
        end
    end
end
%% 2 harmonization
IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
 addpath('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/codes');
 parpool(4);
for ses =1 :2 
        parfor i_Index = 1:length(IndexName)
            raw = [OutDir,'/',ResultsSet{ses},'/',IndexName{i_Index},'_raw.mat'];
            output = [OutDir,'/',ResultsSet{ses}];
            fprintf('harmonize! \n');
            ComBat_CoRR(raw,output,ses,IndexName{i_Index})
        end
        fprintf('harmonize complete! \n');
 end

%% 3 map back to nii
IndexName ={'ALFF_FunImgARCW','DegreeCentrality_FunImgARCWF'};
MeasurePrefixSet={'szALFFMap_','szDegreeCentrality_PositiveWeightedSumBrainMap_'};
for iMeasure=1:length(IndexName) 
    for iSub=1:length(SubID)
       OutList{iSub,1}=[OutDir,'/',IndexName{iMeasure},'/',MeasurePrefixSet{iMeasure},SubID{iSub},'.nii'];
       header.fname = OutList{iSub,1};
       blank = zeros(size(newmap));
       blank(x) = AllData(iSub,:);
       harmonized = reshape(blank,mask_size);
       y_Write(harmonized,header,OutList{iSub,1});   
    end 
end
