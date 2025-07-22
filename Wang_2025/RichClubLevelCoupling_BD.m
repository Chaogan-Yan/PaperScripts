%This script is used to caculate the rich club level SFC
%Each subject will be given 3 values representing the degree of SFC in 3 connection types

%%
%Find the rich club nodes and non-rich club nodes
clc;clear;

InputDir = 'D:\BD_SCFC_Coupling\All_Significant_Files_Final\Results\ResultsRichClub\SC_GTA_Results\HC\';%HC
OutputDir = '';%HC
load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\SubID_Reserved\SubID_Reserved_Final_HC.mat');
SubID_Reserved_Final=SubID;
SubjectOrderKey=size(SubID_Reserved_Final,1);

% InputDir = 'D:\BD_SCFC_Coupling\All_Significant_Files_Final\Results\ResultsRichClub\SC_GTA_Results\BD\';%BD
% OutputDir = '';%BD
% load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\SubID_Reserved\SubID_Reserved_Final_BD.mat');
% SubID_Reserved_Final=SubID;
% SubjectOrderKey=size(SubID_Reserved_Final,1);


RichClub_Node=zeros(SubjectOrderKey,78);
for k = 1:SubjectOrderKey
    load([InputDir,'GTA_',SubID_Reserved_Final{k},'.mat']);

    Mean_Degree_AUC = mean(Degree_AUC);
    SD_Degree_AUC= std(Degree_AUC);
    RichClub_Index=find(Degree_AUC>(Mean_Degree_AUC+SD_Degree_AUC));
    RichClub_Num(k)=length(RichClub_Index);

    RichClub_Node(k,1:length(RichClub_Index))=RichClub_Index;
end


% Common_Node = RichClub_Node(1,:);
% for i = 2:size(RichClub_Node, 1)
%     Common_Node = intersect(Common_Node, RichClub_Node(i,:));
% end

%%

clc;clear;

load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\Results\E_ResultsRichClub\RichClub_Nodes\BD\RichClub_Node.mat');
load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_BD\FC_edges.mat');
load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_BD\SC_edges.mat');

% load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\Results\E_ResultsRichClub\RichClub_Nodes\HC\RichClub_Node.mat');
% load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_HC\FC_edges.mat');
% load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_HC\SC_edges.mat');

nSub=size(RichClub_Node,1);
nRegion=454;

SC_matrix=zeros(nRegion,nRegion,nSub);
for k = 1:nSub
    SC_matrix(:, :,k) = squareform(SC_edges(k,:));
end

FC_matrix=zeros(nRegion,nRegion,nSub);
for k = 1:nSub
    FC_matrix(:, :,k) = squareform(FC_edges(k,:));
end

RichClub_Level_Coupling=zeros(nSub,3);%Hub, Feeder, Local

for i=1:nSub

SC_matrix_i=SC_matrix(:,:,i);
FC_matrix_i=FC_matrix(:,:,i);

RichClub_Node_i=RichClub_Node(i,:);
RichClub_Node_i(RichClub_Node_i==0)=[];

Non_RichClub_Node_i=setdiff(1:nRegion, RichClub_Node_i);

SC_matrix_i_Hub=SC_matrix_i(RichClub_Node_i,RichClub_Node_i);%对称矩阵，用squareform
FC_matrix_i_Hub=FC_matrix_i(RichClub_Node_i,RichClub_Node_i);
SC_edges_i_Hub=squareform(SC_matrix_i_Hub);
FC_edges_i_Hub=squareform(FC_matrix_i_Hub);
zero_SC_edges_i_Hub_Index = find(SC_edges_i_Hub==0);
SC_edges_i_Hub(zero_SC_edges_i_Hub_Index)=[];
FC_edges_i_Hub(zero_SC_edges_i_Hub_Index)=[];

SC_edges_i_Hub=zscore(SC_edges_i_Hub,0,2);
[r,p] = corr(SC_edges_i_Hub',FC_edges_i_Hub', 'type', 'Pearson');
RichClub_Level_Coupling(i,1)=r;


SC_matrix_i_Feeder=SC_matrix_i(RichClub_Node_i,Non_RichClub_Node_i);%非对称，用reshape
FC_matrix_i_Feeder=FC_matrix_i(RichClub_Node_i,Non_RichClub_Node_i);
SC_edges_i_Feeder = reshape(SC_matrix_i_Feeder.',1,[]);
FC_edges_i_Feeder = reshape(FC_matrix_i_Feeder.',1,[]);
zero_SC_edges_i_Feeder_Index = find(SC_edges_i_Feeder==0);
SC_edges_i_Feeder(zero_SC_edges_i_Feeder_Index)=[];
FC_edges_i_Feeder(zero_SC_edges_i_Feeder_Index)=[];

SC_edges_i_Feeder=zscore(SC_edges_i_Feeder,0,2);
[r,p] = corr(SC_edges_i_Feeder',FC_edges_i_Feeder', 'type', 'Pearson');
RichClub_Level_Coupling(i,2)=r;


SC_matrix_i_Local=SC_matrix_i(Non_RichClub_Node_i,Non_RichClub_Node_i);%对称矩阵，用squareform
FC_matrix_i_Local=FC_matrix_i(Non_RichClub_Node_i,Non_RichClub_Node_i);
SC_edges_i_Local=squareform(SC_matrix_i_Local);
FC_edges_i_Local=squareform(FC_matrix_i_Local);
zero_SC_edges_i_Local_Index = find(SC_edges_i_Local==0);
SC_edges_i_Local(zero_SC_edges_i_Local_Index)=[];
FC_edges_i_Local(zero_SC_edges_i_Local_Index)=[];

SC_edges_i_Local=zscore(SC_edges_i_Local,0,2);
[r,p] = corr(SC_edges_i_Local',FC_edges_i_Local', 'type', 'Pearson');
RichClub_Level_Coupling(i,3)=r;

end








