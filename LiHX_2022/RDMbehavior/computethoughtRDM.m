clear;clc
path = '/Users/lihuixian/Documents/1researches/myself/32020experiment/2020Analysis/Report/Textanalysis/Textdata/vectortext';
outpath='/Users/lihuixian/Documents/1researches/myself/32020experiment/2020Analysis/Report/RSAanalysis/bhvRDMresult/';

mkdir(fullfile(outpath,'matfiles'))
mkdir(fullfile(outpath,'csvfiles'))

SubID=dir([path,'/sub*']);

for isub=1:size(SubID,1)
    
    subpath=fullfile(path,SubID(isub).name);
    vectorsdata = readmatrix(subpath);
    
    D=pdist(vectorsdata,'cosine');
    Dmatrix=squareform(D);
    
    savepath=[outpath,'matfiles','/','RDM_',SubID(isub).name(1:6),'.mat'];
    save(savepath,'Dmatrix')
    
    csvsavepath=[outpath,'csvfiles','/','RDM_',SubID(isub).name(1:6),'.csv'];
    csvwrite(csvsavepath,Dmatrix);
end
    