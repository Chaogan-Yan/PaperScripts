% This script is for reading Weibo data downloaded from 
%  the website of “The 3rd CFF Conference on Natural Language Processing & Chinese Computing” 
% (http://tcci.ccf.org.cn/conference/2014/pages/page04_ans.html) and
% experiment data.
% Authur: Bin Lu lub@psych.ac.cn

%% material 1

xDoc = xmlread('/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/evtestdata1/Training data for Emotion Classification.xml');
xRoot=xDoc.getDocumentElement();
SentList= xRoot.getElementsByTagName('sentence');
Num_Sent=SentList.getLength();

for i = 1:Num_Sent
    Sent = SentList.item(i-1);
    Opinion{i} = char(Sent.getAttribute('opinionated'));
    if strcmp(Opinion{i},'Y')
        Emotion_1{i} =  char(Sent.getAttribute('emotion-1-type'));
        Emotion_2{i} =  char(Sent.getAttribute('emotion-2-type'));
    else
        Emotion_1{i} =  char('none');
        Emotion_2{i} =  char('none');
    end
    Content{i}=char(Sent.getTextContent());
    disp(num2str(i))
end

EmoUniq_1 = unique(Emotion_1);
% EmoUniq_1 = {'anger','disgust','fear','happiness','like','none','sadness','surprise'}
for i = 1:length(Emotion_1)
    for j = 1:length(EmoUniq_1)
        if strcmp(Emotion_1{i},EmoUniq_1{j})
            Emotion_1_Num(i,1) = j;
        end
    end
end

OutDir = '/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/WorkDir';
mkdir(OutDir);
save([OutDir,filesep,'Emotion_1_Num_1'],'Emotion_1_Num');
save([OutDir,filesep,'Content_1'],'Content');
clear;



%% material 2

xDoc = xmlread('/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/Emotion Analysis in Chinese Weibo Texts/EmotionClassficationTest.xml');
xRoot=xDoc.getDocumentElement();
SentList= xRoot.getElementsByTagName('sentence');
Num_Sent=SentList.getLength();

for i = 1:Num_Sent
    Sent = SentList.item(i-1);
    Opinion{i} = char(Sent.getAttribute('opinionated'));
    if strcmp(Opinion{i},'Y')
        Emotion_1{i} =  char(Sent.getAttribute('emotion-1-type'));
        Emotion_2{i} =  char(Sent.getAttribute('emotion-2-type'));
    else
        Emotion_1{i} =  char('none');
        Emotion_2{i} =  char('none');
    end
    Content{i}=char(Sent.getTextContent());
    disp(num2str(i))
end

EmoUniq_1 = unique(Emotion_1);
% EmoUniq_1 = {'anger','disgust','fear','happiness','like','none','sadness','surprise'}
for i = 1:length(Emotion_1)
    for j = 1:length(EmoUniq_1)
        if strcmp(Emotion_1{i},EmoUniq_1{j})
            Emotion_1_Num(i,1) = j;
        end
    end
end

OutDir = '/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/WorkDir';
mkdir(OutDir);
save([OutDir,filesep,'Emotion_1_Num_2'],'Emotion_1_Num');
save([OutDir,filesep,'Content_2'],'Content');
clear;

%% material 3

xDoc = xmlread('/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/NLPCC2014微博情绪分析样例数据.xml');
xRoot=xDoc.getDocumentElement();
SentList= xRoot.getElementsByTagName('sentence');
Num_Sent=SentList.getLength();

for i = 1:Num_Sent
    Sent = SentList.item(i-1);
    Opinion{i} = char(Sent.getAttribute('opinionated'));
    if strcmp(Opinion{i},'Y')
        Emotion_1{i} =  char(Sent.getAttribute('emotion-1-type'));
        Emotion_2{i} =  char(Sent.getAttribute('emotion-2-type'));
    else
        Emotion_1{i} =  char('none');
        Emotion_2{i} =  char('none');
    end
    Content{i}=char(Sent.getTextContent());
    disp(num2str(i))
end

EmoUniq_1 = unique(Emotion_1);
% EmoUniq_1 = {'anger','disgust','fear','happiness','like','none','sadness','surprise'}
for i = 1:length(Emotion_1)
    for j = 1:length(EmoUniq_1)
        if strcmp(Emotion_1{i},EmoUniq_1{j})
            Emotion_1_Num(i,1) = j;
        end
    end
end

OutDir = '/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/WorkDir';
mkdir(OutDir);
save([OutDir,filesep,'Emotion_1_Num_3'],'Emotion_1_Num');
save([OutDir,filesep,'Content_3'],'Content');
clear;


%% Integration and balance the number of sentences of each emotion type
EmoAll = [];
ContentAll = {};
for i = 1:3
    load(['/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/WorkDir',filesep,'Content_',num2str(i)]);
    load(['/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/WorkDir',filesep,'Emotion_1_Num_',num2str(i)]);
    EmoAll = [EmoAll;Emotion_1_Num];
    ContentAll = [ContentAll;Content'];
end
EmoNum  = zeros(8,1);
for i = 1:length(EmoAll)
    EmoNum(EmoAll(i)) = EmoNum(EmoAll(i))+1;
end
MeanEmoNum = (sum(EmoNum(1:5))+sum(EmoNum(7:8)))/7;
EmoBanl = [];
ContentBanl = {};
for i = 1:8
    if i ~= 6
        EmoBanl = [EmoBanl;EmoAll(find(EmoAll == i))];
        ContentBanl = [ContentBanl;ContentAll(find(EmoAll == i))];
    else
        Index = find(EmoAll == i);
        EmoBanl = [EmoBanl;EmoAll(Index(1:MeanEmoNum))];
        ContentBanl = [ContentBanl;ContentAll(Index(1:MeanEmoNum))];
    end
end

Index = randi(length(EmoBanl),size(EmoBanl));
EmoBanl = EmoBanl(Index);
ContentBanl = ContentBanl(Index);

OutDir = '/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/WorkDir';
mkdir(OutDir);
save([OutDir,filesep,'Emotion_1_Num_Banl'],'EmoBanl');
save([OutDir,filesep,'Content_Banl'],'ContentBanl');   
    
    
%% check acc of a random classification acc
Rand = randi(7,length(Emotion_1_Num),1);
Result = Rand-Emotion_1_Num;
acc = length(find(Result==0))/length(Emotion_1_Num)


%% Experient data
InputDir = {'/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/FinalText/first10.2',...
    '/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/FinalText/first10.3',...
    '/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/FinalText/second10.2',...
    '/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/FinalText/second10.3'};
ContentName = {'Content_Li_1_2','Content_Li_1_3','Content_Li_2_2','Content_Li_2_3'};
IndexName = {'SubIndex_Li_1_2','SubIndex_Li_1_3','SubIndex_Li_2_2','SubIndex_Li_2_3'};

for iPart = 1:length(InputDir)
    SubList = dir([InputDir{iPart},filesep,'sub*']);
    TextAll = {};
    SubIndex = [];
    n=1;
    for iSub = 1:length(SubList)
        fid1 = fopen([InputDir{iPart},filesep,SubList(iSub).name],'r','n','utf-8');
        TextTemp = textscan(fid1,'%s',inf,'Delimiter','\n');
        fclose(fid1);
        Index = n-1+length(TextTemp{1});
        TextAll((n:Index),1) = TextTemp{1}(1:length(TextTemp{1}));
        SubIndex(iSub,1) = Index;
        n = Index +1;
    end
    OutDir = '/Users/Bin/Projects/LiHuiXian/Project_EmotionPredict/WorkDir';
    mkdir(OutDir);
    save([OutDir,filesep,ContentName{iPart}],'TextAll');
    save([OutDir,filesep,IndexName{iPart}],'SubIndex');
end
    



