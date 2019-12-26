clc;clear

load Content_Weibo.mat

for i=1: size(ContentBanl,1)
    
    text = ContentBanl{i,1};
    
    
    A = strfind(text,'[');
    B = strfind(text,']');
    
    C = strfind(text,'?');
    D = strfind(text,'?');
    
    E = strfind(text,'??');
    F = strfind(text,':');
    
    
    G = strfind(text,'/');
    P = strfind(text,'#');
    
    Index1 = [];
    if length(A)==length(B)
        for j = 1:length(A)
            Index1 = [Index1,A(j):B(j)];
        end
    else
        disp('ssss')
    end
    
    Index2 = [];
    if length(C)==length(D)
        for j = 1:length(C)
            Index2 = [Index2,C(j):D(j)];
        end
    else
        disp('ssss')
    end
    
    Index3=E:F;
    
    Index = [G,Index1,Index2,Index3,P];
    
    text(Index) = '';
     
    ContentBanl{i,1}=text;
    Index = [];
end



save('/Users/lihuixian/Documents/1researches/myself/1behaviralexperiment/2019paper/paper-new/NLP-Lu/EmotionPrediction1210-Lu/AnalysisLi/Weibotext/weibosort1.mat','ContentBanl')


%% 
clc;clear

load weibosort1.mat

%delet space begin
for ii=1: size(ContentBanl,1)
    
    text2 = ContentBanl{ii,1};
    ContentBanl{ii,1}=strtrim(text2);
end
    
 for jj=1: size(ContentBanl,1)
     
     text3 = ContentBanl{jj,1};
     
    
     AA = strfind(text3,'@');
     BB = strfind(text3,':');
     
     EE = strfind(text3,'??');
     FF = strfind(text3,':');
    
     AB = AA:BB;
     
     EF = EE:FF;
     Index = [AB,EF];
     text3(Index) = '';
     
     ContentBanl{jj,1}=text3;
     Index = [];        
 end
 
save('/Users/lihuixian/Documents/1researches/myself/1behaviralexperiment/2019paper/paper-new/NLP-Lu/EmotionPrediction1210-Lu/AnalysisLi/Weibotext/weibosort2.mat','ContentBanl')
 
%%
clc;clear 

load weibosort2.mat
    
for k=1: size(ContentBanl,1)
     
     text4 = ContentBanl{k,1};
     
    
     AAA = strfind(text4,'@');
     BBB = strfind(text4,'?');
     
     EEE = strfind(text4,'??');
     FFF = strfind(text4,':');
     
     ABnew = AAA:BBB;
     
     EFnew = EEE:FFF;
     
     Index = [ABnew,EFnew];
     
     text4(Index) = '';
     
     ContentBanl{k,1}=text4;        
end

save('/Users/lihuixian/Documents/1researches/myself/1behaviralexperiment/2019paper/paper-new/NLP-Lu/EmotionPrediction1210-Lu/AnalysisLi/Weibotext/weibosort3.mat','ContentBanl')


%%
clc;clear

load weibosort3.mat

%delet space begin
for kk=1: size(ContentBanl,1)
    
    text5 = ContentBanl{kk,1};
    ContentBanl{kk,1}=strtrim(text5);
end

for g=1: size(ContentBanl,1)
     
     text5 = ContentBanl{g,1};
     
    
     AAAA = strfind(text5,'@');
     BBBB = strfind(text5,' ');
     
     ABnew1 = AAAA:BBBB;
     
     text5(ABnew1) = '';
     
     ContentBanl{g,1}=text5;        
end

save('/Users/lihuixian/Documents/1researches/myself/1behaviralexperiment/2019paper/paper-new/NLP-Lu/EmotionPrediction1210-Lu/AnalysisLi/Weibotext/weibosort4.mat','ContentBanl')

%%
clc;clear
load weibosort4.mat

for h=1: size(ContentBanl,1)
     
     text6 = ContentBanl{h,1};
     
    
     AAAA = strfind(text6,'@');
     
     Index = AAAA:length(text6);
     
     text6(Index) = '';
     
     ContentBanl{h,1}=text6;        
end

save('/Users/lihuixian/Documents/1researches/myself/1behaviralexperiment/2019paper/paper-new/NLP-Lu/EmotionPrediction1210-Lu/AnalysisLi/Weibotext/weibosort5.mat','ContentBanl')










        
            
            
        

    
    
    