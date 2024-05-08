clear;clc
outpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020FrameLine/2022MVPA/convergentwoexperienceMVPA/sortdataFour';
mkdir(outpath)
path='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020FrameLine/2022MVPA/convergentwoexperienceMVPA';
cd(path)
load convergenceSubid.mat

FLpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020FrameLine/2022MVPA/firstlevel';
Reparh='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020Relationship/MVPA/firstlevel';
for isub=1:size(convergenceSubid,1)
    
    subid=convergenceSubid{isub};
    FLsubpath=fullfile(FLpath,subid);
    [maskdata,maskheader]=y_Read([FLsubpath,'/mask.nii']);
    sortdata=zeros([size(maskdata),16]);

    % FL task
    Index=1:2:16;
    %beatinfo:1-6: Absolute_Congruent; 7-12:Absolute_Incongruent; 13-18:Relative_Congruent; 19-24:Relative_Incongruent
    conindx=[1:12;13:24];
    %indx=conindx(:);
    % 1-6;13-18 Congruent; 7-12;19-24 Incongruent random select 2
    indx_1=randperm(6,2);
    indx_2=randperm(6,2)+6;
    Sconindx=conindx(:,[indx_1,indx_2]);
    indxFL=Sconindx(:);
     for i=1:size(indxFL,1)
        ibeat=indxFL(i);
        if ibeat<10
            [subdataFL,subheader]=y_Read([FLsubpath,'/beta_000',num2str(ibeat),'.nii']);
        else
            [subdataFL,subheader]=y_Read([FLsubpath,'/beta_00',num2str(ibeat),'.nii']);
        end
        subdataFL(find(isnan(subdataFL))) = 0;
        sortdata(:,:,:,Index(i))=subdataFL;
     end
    
    % Relationship task
    Index2=2:2:16;
    Resubpath=fullfile(Reparh,subid);
    %beatinfo:
    %1-4: FreeChoice; 5-8:ForceA; 9-12:ForceH; 13-16: Control 
    conindxR=[5:8;9:12];
    indxR=conindxR(:);
    for j=1:size(indxR,1)
        ibeatR=indxR(j);
        if ibeatR<10
            [subdataR,subheader]=y_Read([Resubpath,'/beta_000',num2str(ibeatR),'.nii']);
        else
            [subdataR,subheader]=y_Read([Resubpath,'/beta_00',num2str(ibeatR),'.nii']);
        end
        subdataR(find(isnan(subdataR))) = 0;
        sortdata(:,:,:,Index2(j))=subdataR;
    end
    
    outname=[outpath,'/',subid,'.nii'];
    y_Write(sortdata,subheader,outname)
end
     
    


    