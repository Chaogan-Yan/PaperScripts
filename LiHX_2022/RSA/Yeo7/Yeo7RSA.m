clear;clc
outpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/WYSWYT/RSAanalysisWYSWYT/RSAnewprepocess/ROIRSAanalysis/1paperYeo7/';
bhvRDMpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/WYSWYT/RSAanalysisWYSWYT/RSAnewprepocess/bhvRDMresult/matfiles/';
bhvRDMSenLenpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/WYSWYT/RSAanalysisWYSWYT/RSAnewprepocess/partialSpearmanAnalysis/bhvRDMsentenceLength/matfiles/';

datapath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/WYSWYT/RSAanalysisWYSWYT/RSAnewprepocess/ROIRSAanalysis/1paper86schaefer400/SpeakEventBeta39';
subid=dir([datapath,'/sub*']);


maskpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/WYSWYT/RSAanalysisWYSWYT/RSAnewprepocess/ROIRSAanalysis/1paperYeo7/Reslice_Yeo2011_7Networks_MNI152.nii';
[MaskData,MaskVox,MaskHead]=y_ReadRPI(maskpath);
MaskROI=reshape(MaskData,1,[]);
Element = unique(MaskROI);
Element(find(isnan(Element))) = [];
Element(find(Element==0)) = [];

for iElement=1:length(Element)
    
    idx_Mask = find(MaskROI==Element(iElement));
    
    SubResults=zeros(size(subid,1), 6);
    for isub=1:size(subid,1)
        subpath=fullfile(datapath,subid(isub).name);
        eventlist=dir([subpath,'/sub*']);
        
        Extractdata = zeros(size(eventlist,1), length(idx_Mask));
        for ie=1:size(eventlist,1)
            if ie<10
                evenname=['beta_000',num2str(ie),'.nii'];
            else
                evenname=['beta_00',num2str(ie),'.nii'];
            end
            eventpath=[subpath,'/',subid(isub).name,'Even',num2str(ie),'_',evenname];
            [eventdata,~]=y_Read(eventpath);
            [dim1, dim2, dim3, tp] = size(eventdata);
            Data = reshape(eventdata, dim1*dim2*dim3, tp)';
            Extractdata(ie,:) = Data(:, idx_Mask);
        end
        %nan--->0
        Extractdata(find(isnan(Extractdata))) = 0;  
        %RDM 1-r
        trannsform_Extractdata=Extractdata';
        D=corrcoef(trannsform_Extractdata);
        NetworkRDM=1-abs(D); 

        b=eye(size(NetworkRDM));
        NetworkRDMxin=NetworkRDM+100*b;
        v1=NetworkRDMxin(tril(true(size(NetworkRDM))));
        v1(find(v1==100))=[];
        
        %%behavioral RDM
        bhvRDM=load([bhvRDMpath,'RDM_',subid(isub).name,'.mat']);
        bhvRDMxin=bhvRDM.Dmatrix+100*b; 
        v2=bhvRDMxin(tril(true(size(bhvRDMxin))));
        v2(find(v2==100))=[];                              
        
        %spearman correlation
        [r,p]=corr(v1,v2,'type','Spearman');
        %r to z
        r_Z = 0.5 .* log( (1+r) ./ (1-r) );
        
        %covariation  sentence length 
        bhvSenLenRDM=load([bhvRDMSenLenpath,'wordNumRDM_',subid(isub).name,'.mat']);
        bhvSenLenRDMxin=bhvSenLenRDM.wordNumchaRDM+10000*b; 
        v3=bhvSenLenRDMxin(tril(true(size(bhvSenLenRDMxin))));
        v3(find(v3==10000))=[];  
        
        %%% partial spearman correlation
        [r_partial,p_partial]=partialcorr(v1,v2,v3,'type','Spearman');
        %r to z
        Zr_partial= 0.5 .* log( (1+r_partial) ./ (1-r_partial) );
        
        SubResults(isub,1)=r;
        SubResults(isub,2)=p; 
        SubResults(isub,3)=r_Z;
        SubResults(isub,4)=r_partial;
        SubResults(isub,5)=p_partial;
        SubResults(isub,6)=Zr_partial;        
    end
    savebpath=[outpath,'E1Yeo',num2str(Element(iElement)),'_SubResultsrp.mat'];
    save(savebpath,'SubResults')
end


