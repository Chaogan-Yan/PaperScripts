function [within_networkFC,between_networkFC,within_networkFCName,between_networkFCName] = c_ExtractDMN_FC(SubList,DataDir,OutputName)
%[within_networkFC,between_networkFC,within_networkFCName,between_networkFCName] = c_Extract_NetworkFC66(SubList,DataDir,OutputName)
%Calculate the FCs within the entire DMN, between or within 3 sub-systems of DMN.
%ROIs were defined following Dixon et al. 2017 NeuroImage.
%a ROI set containing 24 ROIs:
%~/DMN_ROI_List.mat
%Inputs: 
%SubList: a n x 1 cell containing the list of subjects
%DataDir: the path FCs stored, should be DPARSFA's <working directory>/Results/***ROISignals_FunImgARCWF* 
%FC should be already calculated with DPARSFA
%OutputName: directory and filename you wish to store the results. If you
%do not need to save a file, please leave a '' here.
%"DN" means "DMN"

%Written by Xiao Chen 20190404

within_networkFC = zeros(length(SubList),4);
between_networkFC = zeros(length(SubList),3);

within_networkFCName = {'DN_core','DN_dmPFC','DN_MTL','DN_total'};
between_networkFCName = {'DN_core-DN_dmPFC','DN_core-DN_MTL','DN_dmPFC-DN_MTL'};
for i = 1:length(SubList)
    %Read in FC matrixes
    load([DataDir,'/ROICorrelation_FisherZ_',SubList{i},'.mat']);
    %%%%%%%%%%%%%%%%%%%%%%% WITHIN_NETWORKS %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate the FC within DN_core
    icount_temp=0;
    sum_temp=0;
    for ii = 1:8
        for jj =  ii+1:9
            sum_temp=sum_temp+ROICorrelation_FisherZ(ii,jj);
            icount_temp=icount_temp+1;
        end
    end
    within_networkFC(i,1)= sum_temp/icount_temp;
    
    %calculate the FC within DN_dmPFC
    icount_temp=0;
    sum_temp=0;
    for ii = 10:17
        for jj =  ii+1:18
            sum_temp=sum_temp+ROICorrelation_FisherZ(ii,jj);
            icount_temp=icount_temp+1;
        end
    end
    within_networkFC(i,2)= sum_temp/icount_temp;
    
    %calculate the FC within DN_MTL
    icount_temp=0;
    sum_temp=0;
    for ii = 19:23
        for jj =  ii+1:24
        sum_temp=sum_temp+ROICorrelation_FisherZ(ii,jj);
        icount_temp=icount_temp+1;
        end
    end
    within_networkFC(i,3)= sum_temp/icount_temp;
    
    %the FC within DN_total: the mean of three within_subsystems' FC
    within_networkFC(i,4) = (within_networkFC(i,1) + within_networkFC(i,2) + within_networkFC(i,3))/3;
    
    %%%%%%%%%%%%%%%%%%%%%%% BETWEEN_NETWORKS %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate FCs between DN_core and DN_dmPFC
    icount_temp=0;
    sum_temp=0;
    for ii = 1:9
        for jj = 10:18
            sum_temp = sum_temp+ROICorrelation_FisherZ(ii,jj);
            icount_temp = icount_temp+1;
        end
    end
    between_networkFC(i,1) = sum_temp/icount_temp;
    
    %calculate FCs between DN_core and DN_MTL
    icount_temp=0;
    sum_temp=0;
    for ii = 1:9
        for jj = 19:24
            sum_temp = sum_temp+ROICorrelation_FisherZ(ii,jj);
            icount_temp = icount_temp+1;
        end
    end
    between_networkFC(i,2) = sum_temp/icount_temp;
    
    %calculate FCs between DN_dmPFC and DN_MTL
    icount_temp=0;
    sum_temp=0;
    for ii = 10:18
       for jj = 19:24
           sum_temp = sum_temp+ROICorrelation_FisherZ(ii,jj);
           icount_temp = icount_temp+1;
       end
    end
    between_networkFC(i,3) = sum_temp/icount_temp;
end
% output a file
if ~isempty(OutputName)
    save(OutputName, '*networkFC','*networkFCName');
end