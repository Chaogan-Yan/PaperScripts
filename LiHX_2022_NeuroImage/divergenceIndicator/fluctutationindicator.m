clear;clc
path = '/Users/lihuixian/Documents/1researches/myself/32020experiment/2020Analysis/Textanalysis/data/Textdata/vectortext';
outpath='/Users/lihuixian/Documents/1researches/myself/32020experiment/2020Analysis/Textanalysis/text';

SubID=dir([path,'/sub*']);
divergence_fluctutation=zeros(size(SubID,1),1);
for isub=1:size(SubID,1)
    
    subpath=fullfile(path,SubID(isub).name);
    vectorsdata = readmatrix(subpath);
    
    [row,col] = size(vectorsdata); 
    Hudu = zeros(row,1);
    for im = 1:row   
        Hudu(im,1) = acos(dot(vectorsdata(im,:),mean(vectorsdata))/(norm(vectorsdata(im,:))*norm(mean(vectorsdata))));  
    end
    
    divergence_fluctutation(isub,1)= (std(Hudu)/mean(Hudu));
    
end
save([outpath,'indexsResults.mat'],'Sub_divergence')