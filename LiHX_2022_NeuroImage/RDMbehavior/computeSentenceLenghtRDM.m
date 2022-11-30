clear;clc
path = '/Users/lihuixian/Documents/1researches/myself/32020experiment/2020Analysis/Report/RSAanalysis/bhvRDMsentenceLength/';
outpath='/Users/lihuixian/Documents/1researches/myself/32020experiment/2020Analysis/Report/RSAanalysis/bhvRDMsentenceLength/';

E=['E1';'E2'];

for iE=1:2
    
    outpathmat=fullfile(outpath,E(iE,:),'matfiles');
    mkdir(outpathmat)
    outpathcsv=fullfile(outpath,E(iE,:),'csvfiles');
    mkdir(outpathcsv)
    [data,txtdata,celldata]=xlsread([path,'wordinfo.xlsx'],E(iE,:));
    
    [data_conversion,txtdata_conversion,celldata_conversion]=xlsread([path,'conversioninfo.xlsx'],E(iE,:));
    
    
    SubIndex=zeros(size(data_conversion,1),1);
    for idex=1:size(data_conversion,1)
        
        SubIndex(idex,1)=sum(data_conversion(1:idex));
    end
    
    SubTextIndex= [0;SubIndex];
   
    for isub = 1:size(SubTextIndex,1)-1
        
        SubWord = data(SubTextIndex(isub)+1:SubTextIndex(isub+1),1);

        row=length(SubWord);
        wordNumcha =zeros(row,row);
        for i = 1:row-1
            for j=(i+1):row
                cha = abs(SubWord(i)-SubWord(j));
                wordNumcha(j,i)=cha;
            end
        end
        wordNumchaRDM=wordNumcha+wordNumcha';

        savepath=[outpathmat,'/','wordNumRDM_',celldata{SubTextIndex(isub)+2}(1:6),'.mat'];
        save(savepath,'wordNumchaRDM')
        
        csvsavepath=[outpathcsv,'/','wordNumRDM_',celldata{SubTextIndex(isub)+2}(1:6),'.csv'];
        writematrix(wordNumchaRDM,csvsavepath);
    end
end
    