clear;clc
path='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/convergentTwoEMeanBold/1Yeo7/';
bhvRDMpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RASanalyisis/bhvRDMresult/matfiles/';
bhvRDMSenLenpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RSAnewprocess/partialSpearmanAnalysis/bhvRDMsentenceLength/matfiles/';

datapath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/MeanBlodAnalysis/MeanBlodSpeakEvent47';
subid=dir([datapath,'/sub*']);


maskpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RASanalyisis/ROIRSAanalysis/1paperYeo7network/Reslice_Yeo2011_7Networks_MNI152.nii';
[MaskData,MaskVox,MaskHead]=y_ReadRPI(maskpath);
MaskROI=reshape(MaskData,1,[]);
Element = unique(MaskROI);
Element(find(isnan(Element))) = [];% ignore background if encoded as nan.
Element(find(Element==0)) = [];

for iElement=1:length(Element)
    
    idx_Mask = find(MaskROI==Element(iElement));
    
    outpath=fullfile(path,['Yeo',num2str(iElement)]);
    mkdir(outpath)
    SubResults=zeros(size(subid,1), 6);
    for isub=1:size(subid,1)
        subpath=fullfile(datapath,subid(isub).name);
        eventlist=dir([subpath,'/MeanBlod*']);
        
        Extractdata = zeros(size(eventlist,1), length(idx_Mask));
        for ie=1:size(eventlist,1)
            if ie<10
                evenname=['MeanBlod_Event00',num2str(ie),'.nii'];
            else
                evenname=['MeanBlod_Event0',num2str(ie),'.nii'];
            end
            eventpath=[subpath,'/',evenname];
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
        %maybe exist 0, so v2(find(v2==0))=[]; --->100
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
        
        %covariation  sentence length may be very long. 
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
    savebpath=[outpath,'/E2Yeo',num2str(Element(iElement)),'_SubResultsrp.mat'];
    save(savebpath,'SubResults')
end


