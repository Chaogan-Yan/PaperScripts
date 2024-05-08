clear;clc
datapath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020Relationship/MVPA/firstlevel';
outpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020Relationship/MVPA/sortdataAH';
mkdir(outpath)

subid=dir([datapath,'/sub*']);

for isub=1:size(subid,1)
    subpath=fullfile(datapath,subid(isub).name);
    
    [maskdata,maskheader]=y_Read([subpath,'/mask.nii']);
    sortdata=zeros([size(maskdata),8]); 
    
    %beatinfo:
    %1-4: FreeChoice; 5-8:ForceA; 9-12:ForceH; 13-16: Control 
    conindx=[5:8;9:12];
    indx=conindx(:);
    for i=1:size(indx,1)
        ibeat=indx(i);
        if ibeat<10
            [subdata,subheader]=y_Read([subpath,'/beta_000',num2str(ibeat),'.nii']);
        else
            [subdata,subheader]=y_Read([subpath,'/beta_00',num2str(ibeat),'.nii']);
        end
        subdata(find(isnan(subdata))) = 0;
        sortdata(:,:,:,i)=subdata;
    end
    
    outname=[outpath,'/',subid(isub).name,'.nii'];
    y_Write(sortdata,subheader,outname)
end


    
