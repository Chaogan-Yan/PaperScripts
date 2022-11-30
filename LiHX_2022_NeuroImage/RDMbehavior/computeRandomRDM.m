clear;clc
path = '/Users/lihuixian/Documents/1researches/myself/32020experiment/2020Analysis/Textanalysis/data/Textdata/vectortext';
outpath='/Users/lihuixian/Documents/1researches/myself/32020experiment/2020Analysis/Report/RSAanalysis/bhvRDMresultRandom/';

mkdir(fullfile(outpath,'matfiles'))
mkdir(fullfile(outpath,'csvfiles'))

SubID=dir([path,'/sub*']);

for isub=1:size(SubID,1)
    
    subpath=fullfile(path,SubID(isub).name);
    vectorsdata = readmatrix(subpath); 
    
    randIndex=randperm(size(vectorsdata,1));
    vectorsdata_new = vectorsdata(randIndex,:);
    
    D=pdist(vectorsdata_new,'cosine');
    Dmatrix=squareform(D);
    
    savepath=[outpath,'matfiles','/','randomRDM_',SubID(isub).name(1:6),'.mat'];
    save(savepath,'Dmatrix')
    
    csvsavepath=[outpath,'csvfiles','/','randomRDM_',SubID(isub).name(1:6),'.csv'];
    csvwrite(csvsavepath,Dmatrix);
end
    