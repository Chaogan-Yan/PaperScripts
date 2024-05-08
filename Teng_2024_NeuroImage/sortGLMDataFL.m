clear;clc
datapath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020FrameLine/2022MVPA/firstlevel';
outpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020FrameLine/2022MVPA/sortdataReAb';
mkdir(outpath)

subid=dir([datapath,'/sub*']);

for isub=1:size(subid,1)
    subpath=fullfile(datapath,subid(isub).name);
    
    [maskdata,maskheader]=y_Read([subpath,'/mask.nii']);
    sortdata=zeros([size(maskdata),24]);
    
    %beatinfo:1-6: Absolute_Congruent; 7-12:Absolute_Incongruent; 13-18:Relative_Congruent; 19-24:Relative_Incongruent
    
    conindx=[1:12;13:24];
    indx=conindx(:);
    
    for i=1:24
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


    
