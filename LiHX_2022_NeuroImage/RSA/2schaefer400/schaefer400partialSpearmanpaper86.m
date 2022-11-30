clear;clc
outpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/convergentTwoEMeanBold/1schaefer86/';
bhvRDMpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RASanalyisis/bhvRDMresult/matfiles/';
bhvRDMSenLenpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RSAnewprocess/partialSpearmanAnalysis/bhvRDMsentenceLength/matfiles/';

datapath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/MeanBlodAnalysis/MeanBlodSpeakEvent47';
subid=dir([datapath,'/sub*']);

maskpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RASanalyisis/ROIRSAanalysis/schaefer400/schaefer/Reslice_Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii';
[MaskData,MaskVox,MaskHead]=y_ReadRPI(maskpath);
MaskROI=reshape(MaskData,1,[]);
Element = unique(MaskROI);
Element(find(isnan(Element))) = [];% ignore background if encoded as nan.
Element(find(Element==0)) = [];

ROI_MultipleLabel = zeros(size(subid,1),length(Element),3);  %3: r ; p; r-z
for iElement=1:length(Element)
    
    idx_Mask = find(MaskROI==Element(iElement)); 
    
    SubResults=zeros(size(subid,1), 2);
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
            %Set the NaN voxels to 0.
            eventdata(find(isnan(eventdata))) = 0;
            % Convert into 2D
            [dim1, dim2, dim3, tp] = size(eventdata);
            Data = reshape(eventdata,[],tp)';
            Extractdata(ie,:) = Data(:, idx_Mask);
        end  
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
        
        %covariation  sentence length may be very long. 
        bhvSenLenRDM=load([bhvRDMSenLenpath,'wordNumRDM_',subid(isub).name,'.mat']);
        bhvSenLenRDMxin=bhvSenLenRDM.wordNumchaRDM+10000*b; 
        v3=bhvSenLenRDMxin(tril(true(size(bhvSenLenRDMxin))));
        v3(find(v3==10000))=[];
        
        % partial spearman correlation
        [r,p]=partialcorr(v1,v2,v3,'type','Spearman');
        %r to z
        r_Z = 0.5 .* log( (1+r) ./ (1-r) );
        
        ROI_MultipleLabel(isub,iElement,1) = r; 
        ROI_MultipleLabel(isub,iElement,2) = p;
        ROI_MultipleLabel(isub,iElement,3) =r_Z;   
    end
end
savebpath=[outpath,'E2PartialSpearman_schaefer400_SubResultsrp.mat'];
save(savebpath,'ROI_MultipleLabel')      
        
        
        
        
        
        
        
        
        
        
        
        
