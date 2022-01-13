clear; clc;

path = '/Users/lihuixian/Documents/1researches/myself/32020experiment/2020Analysis/Textanalysis/data/Textdata/vectortext';
outpath='/Users/lihuixian/Documents/1researches/myself/32020experiment/2020Analysis/Textanalysis/text';
SubID = dir([path,'/sub*']);

[B5minNum,txt5min,cellB5minNum] = xlsread('/Users/lihuixian/Documents/1researches/myself/32020experiment/2020Analysis/Textanalysis/text/5minCut.xlsx');

huducha_5mincut=zeros(size(SubID,1),2);
for isub = 1:size(SubID,1)
    
    subpath=fullfile(path,SubID(isub).name);
    vectorsdata = readmatrix(subpath);
    
    subid =str2double(SubID(isub).name(isstrprop(SubID(isub).name,'digit')));
    [Row,Colu] = find(B5minNum(:,1)==subid);
    
    B5minMeanVectors = mean(vectorsdata(1:B5minNum(Row,2),:));
    %%sub016 5items, A5:1item
    if subid==16
        A5minMeanVectors=vectorsdata(end,:);
    else
        A5minMeanVectors = mean(vectorsdata(B5minNum(Row,2)+1:end,:));
    end
    
    huducha = acos(dot(B5minMeanVectors,A5minMeanVectors)/(norm(B5minMeanVectors)*norm(A5minMeanVectors)));
    
    huducha_5mincut(isub,1)=subid;
    huducha_5mincut(isub,2)=huducha;
end
save([outpath,'/huducha_5mincut.mat'],'huducha_5mincut');
