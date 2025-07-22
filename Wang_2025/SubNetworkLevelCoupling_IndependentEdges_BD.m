%This script is used to caculate the subnetwork level SC-FC coupling
%Each subject will be given 8 values representing the degree of their brain SC-FC coupling in each subnetwork


clc;clear;

% %Load the BD data file
%  load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_BD\FC_edges.mat');
%  load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_BD\SC_edges.mat');
%  OutDir = 'D:\BD_SCFC_Coupling\SensitivityAnalysis\Results\SubNetworkSplitLR\BD\';
%Load the HC data file
load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_HC\FC_edges.mat');
load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_HC\SC_edges.mat');
OutDir = 'D:\BD_SCFC_Coupling\SensitivityAnalysis\Results\SubNetworkSplitLR\HC\';

%Load the network index and label information
load('D:\DPABI\DPABI_V7.0_ForCamp\DPABISurf\SurfTemplates\DPABISurf_Schaefer2018_400_Tian2020_54_Info.mat');


CorrType='Pearson';

%Converting data from vector form to symmetric matrix form
nreg = 454; % Yeo400+Tian54
nedges = nreg*(nreg-1)/2;
nsub=size(SC_edges,1);

SC_matrix=zeros(nreg,nreg,nsub);
for k = 1:nsub
    SC_matrix(:, :,k) = squareform(SC_edges(k,:));
end

FC_matrix=zeros(nreg,nreg,nsub);
for k = 1:nsub
    FC_matrix(:, :,k) = squareform(FC_edges(k,:));
end

%%Caculate the subnetwork level SC-FC coupling
VN_Index=find(DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork==1);
VN_Index_lh=VN_Index(1:31);%1:ceil(length(VN_Index)/2)
VN_Index_rh=VN_Index(32:61);%ceil(length(VN_Index)/2)+1:end


SMN_Index=find(DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork==2);
SMN_Index_lh=SMN_Index(1:37);
SMN_Index_rh=SMN_Index(38:77);

DAN_Index=find(DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork==3);
DAN_Index_lh=DAN_Index(1:23);
DAN_Index_rh=DAN_Index(24:46);

VAN_Index=find(DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork==4);
VAN_Index_lh=VAN_Index(1:22);
VAN_Index_rh=VAN_Index(23:47);

LN_Index=find(DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork==5);
LN_Index_lh=LN_Index(1:13);
LN_Index_rh=LN_Index(14:26);

FPN_Index=find(DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork==6);
FPN_Index_lh=FPN_Index(1:22);
FPN_Index_rh=FPN_Index(23:52);

DMN_Index=find(DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork==7);
DMN_Index_lh=DMN_Index(1:52);
DMN_Index_rh=DMN_Index(53:91);

SC_Index=find(DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork==8);
SC_Index_lh=SC_Index(28:54);
SC_Index_rh=SC_Index(1:27);

for i=1:nsub

VN_SC_lh=SC_matrix(VN_Index_lh, :,i);
SMN_SC_lh=SC_matrix(SMN_Index_lh, :,i);
DAN_SC_lh=SC_matrix(DAN_Index_lh, :,i);
VAN_SC_lh=SC_matrix(VAN_Index_lh, :,i);
LN_SC_lh=SC_matrix(LN_Index_lh, :,i);
FPN_SC_lh=SC_matrix(FPN_Index_lh, :,i);
DMN_SC_lh=SC_matrix(DMN_Index_lh, :,i);
SC_SC_lh=SC_matrix(SC_Index_lh, :,i);

VN_FC_lh=FC_matrix(VN_Index_lh, :,i);
SMN_FC_lh=FC_matrix(SMN_Index_lh, :,i);
DAN_FC_lh=FC_matrix(DAN_Index_lh, :,i);
VAN_FC_lh=FC_matrix(VAN_Index_lh, :,i);
LN_FC_lh=FC_matrix(LN_Index_lh, :,i);
FPN_FC_lh=FC_matrix(FPN_Index_lh, :,i);
DMN_FC_lh=FC_matrix(DMN_Index_lh, :,i);
SC_FC_lh=FC_matrix(SC_Index_lh, :,i);



VN_SC_rh=SC_matrix(VN_Index_rh, :,i);
SMN_SC_rh=SC_matrix(SMN_Index_rh, :,i);
DAN_SC_rh=SC_matrix(DAN_Index_rh, :,i);
VAN_SC_rh=SC_matrix(VAN_Index_rh, :,i);
LN_SC_rh=SC_matrix(LN_Index_rh, :,i);
FPN_SC_rh=SC_matrix(FPN_Index_rh, :,i);
DMN_SC_rh=SC_matrix(DMN_Index_rh, :,i);
SC_SC_rh=SC_matrix(SC_Index_rh, :,i);

VN_FC_rh=FC_matrix(VN_Index_rh, :,i);
SMN_FC_rh=FC_matrix(SMN_Index_rh, :,i);
DAN_FC_rh=FC_matrix(DAN_Index_rh, :,i);
VAN_FC_rh=FC_matrix(VAN_Index_rh, :,i);
LN_FC_rh=FC_matrix(LN_Index_rh, :,i);
FPN_FC_rh=FC_matrix(FPN_Index_rh, :,i);
DMN_FC_rh=FC_matrix(DMN_Index_rh, :,i);
SC_FC_rh=FC_matrix(SC_Index_rh, :,i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SC
VN_SC_lh_Move=VN_SC_lh;
VN_SC_lh_edges=hl_MatrixToIndependentEdgeVector(VN_SC_lh_Move);

SMN_SC_lh_Move=SMN_SC_lh(:,[32:end,1:31]);
SMN_SC_lh_edges=hl_MatrixToIndependentEdgeVector(SMN_SC_lh_Move);

DAN_SC_lh_Move=DAN_SC_lh(:,[69:end,1:68]);
DAN_SC_lh_edges=hl_MatrixToIndependentEdgeVector(DAN_SC_lh_Move);

VAN_SC_lh_Move=VAN_SC_lh(:,[92:end,1:91]);
VAN_SC_lh_edges=hl_MatrixToIndependentEdgeVector(VAN_SC_lh_Move);

LN_SC_lh_Move=LN_SC_lh(:,[114:end,1:113]);
LN_SC_lh_edges=hl_MatrixToIndependentEdgeVector(LN_SC_lh_Move);

FPN_SC_lh_Move=FPN_SC_lh(:,[127:end,1:126]);
FPN_SC_lh_edges=hl_MatrixToIndependentEdgeVector(FPN_SC_lh_Move);

DMN_SC_lh_Move=DMN_SC_lh(:,[149:end,1:148]);
DMN_SC_lh_edges=hl_MatrixToIndependentEdgeVector(DMN_SC_lh_Move);

SC_SC_lh_Move=SC_SC_lh(:,[428:end,1:427]);
SC_SC_lh_edges=hl_MatrixToIndependentEdgeVector(SC_SC_lh_Move);

%------------------------------------------------------------------------------------------------------------
VN_SC_rh_Move=VN_SC_rh(:,[201:end,1:200]);
VN_SC_rh_edges=hl_MatrixToIndependentEdgeVector(VN_SC_rh_Move);

SMN_SC_rh_Move=SMN_SC_rh(:,[231:end,1:230]);
SMN_SC_rh_edges=hl_MatrixToIndependentEdgeVector(SMN_SC_rh_Move);

DAN_SC_rh_Move=DAN_SC_rh(:,[271:end,1:270]);
DAN_SC_rh_edges=hl_MatrixToIndependentEdgeVector(DAN_SC_rh_Move);

VAN_SC_rh_Move=VAN_SC_rh(:,[294:end,1:293]);
VAN_SC_rh_edges=hl_MatrixToIndependentEdgeVector(VAN_SC_rh_Move);

LN_SC_rh_Move=LN_SC_rh(:,[319:end,1:318]);
LN_SC_rh_edges=hl_MatrixToIndependentEdgeVector(LN_SC_rh_Move);

FPN_SC_rh_Move=FPN_SC_rh(:,[332:end,1:331]);
FPN_SC_rh_edges=hl_MatrixToIndependentEdgeVector(FPN_SC_rh_Move);

DMN_SC_rh_Move=DMN_SC_rh(:,[362:end,1:361]);
DMN_SC_rh_edges=hl_MatrixToIndependentEdgeVector(DMN_SC_rh_Move);

SC_SC_rh_Move=SC_SC_rh(:,[401:end,1:400]);
SC_SC_rh_edges=hl_MatrixToIndependentEdgeVector(SC_SC_rh_Move);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FC
VN_FC_lh_Move=VN_FC_lh;
VN_FC_lh_edges=hl_MatrixToIndependentEdgeVector(VN_FC_lh_Move);

SMN_FC_lh_Move=SMN_FC_lh(:,[32:end,1:31]);
SMN_FC_lh_edges=hl_MatrixToIndependentEdgeVector(SMN_FC_lh_Move);

DAN_FC_lh_Move=DAN_FC_lh(:,[69:end,1:68]);
DAN_FC_lh_edges=hl_MatrixToIndependentEdgeVector(DAN_FC_lh_Move);

VAN_FC_lh_Move=VAN_FC_lh(:,[92:end,1:91]);
VAN_FC_lh_edges=hl_MatrixToIndependentEdgeVector(VAN_FC_lh_Move);

LN_FC_lh_Move=LN_FC_lh(:,[114:end,1:113]);
LN_FC_lh_edges=hl_MatrixToIndependentEdgeVector(LN_FC_lh_Move);

FPN_FC_lh_Move=FPN_FC_lh(:,[127:end,1:126]);
FPN_FC_lh_edges=hl_MatrixToIndependentEdgeVector(FPN_FC_lh_Move);

DMN_FC_lh_Move=DMN_FC_lh(:,[149:end,1:148]);
DMN_FC_lh_edges=hl_MatrixToIndependentEdgeVector(DMN_FC_lh_Move);

SC_FC_lh_Move=SC_FC_lh(:,[428:end,1:427]);
SC_FC_lh_edges=hl_MatrixToIndependentEdgeVector(SC_FC_lh_Move);

%------------------------------------------------------------------------------------------------------------
VN_FC_rh_Move=VN_FC_rh(:,[201:end,1:200]);
VN_FC_rh_edges=hl_MatrixToIndependentEdgeVector(VN_FC_rh_Move);

SMN_FC_rh_Move=SMN_FC_rh(:,[231:end,1:230]);
SMN_FC_rh_edges=hl_MatrixToIndependentEdgeVector(SMN_FC_rh_Move);

DAN_FC_rh_Move=DAN_FC_rh(:,[271:end,1:270]);
DAN_FC_rh_edges=hl_MatrixToIndependentEdgeVector(DAN_FC_rh_Move);

VAN_FC_rh_Move=VAN_FC_rh(:,[294:end,1:293]);
VAN_FC_rh_edges=hl_MatrixToIndependentEdgeVector(VAN_FC_rh_Move);

LN_FC_rh_Move=LN_FC_rh(:,[319:end,1:318]);
LN_FC_rh_edges=hl_MatrixToIndependentEdgeVector(LN_FC_rh_Move);

FPN_FC_rh_Move=FPN_FC_rh(:,[332:end,1:331]);
FPN_FC_rh_edges=hl_MatrixToIndependentEdgeVector(FPN_FC_rh_Move);

DMN_FC_rh_Move=DMN_FC_rh(:,[362:end,1:361]);
DMN_FC_rh_edges=hl_MatrixToIndependentEdgeVector(DMN_FC_rh_Move);

SC_FC_rh_Move=SC_FC_rh(:,[401:end,1:400]);
SC_FC_rh_edges=hl_MatrixToIndependentEdgeVector(SC_FC_rh_Move);

%找到SC零值索引
zero_VN_SC_lh_edges_index = find(VN_SC_lh_edges==0);
zero_VN_SC_rh_edges_index = find(VN_SC_rh_edges==0);

zero_SMN_SC_lh_edges_index = find(SMN_SC_lh_edges==0);
zero_SMN_SC_rh_edges_index = find(SMN_SC_rh_edges==0);

zero_DAN_SC_lh_edges_index = find(DAN_SC_lh_edges==0);
zero_DAN_SC_rh_edges_index = find(DAN_SC_rh_edges==0);

zero_VAN_SC_lh_edges_index = find(VAN_SC_lh_edges==0);
zero_VAN_SC_rh_edges_index = find(VAN_SC_rh_edges==0);

zero_LN_SC_lh_edges_index = find(LN_SC_lh_edges==0);
zero_LN_SC_rh_edges_index = find(LN_SC_rh_edges==0);

zero_FPN_SC_lh_edges_index = find(FPN_SC_lh_edges==0);
zero_FPN_SC_rh_edges_index = find(FPN_SC_rh_edges==0);

zero_DMN_SC_lh_edges_index = find(DMN_SC_lh_edges==0);
zero_DMN_SC_rh_edges_index = find(DMN_SC_rh_edges==0);

zero_SC_SC_lh_edges_index = find(SC_SC_lh_edges==0);
zero_SC_SC_rh_edges_index = find(SC_SC_rh_edges==0);

%按SC零值索引将SC和FC向量赋空值
%lh
VN_SC_lh_edges(zero_VN_SC_lh_edges_index) = [];
VN_SC_lh_edges=zscore(VN_SC_lh_edges,0,2);
VN_FC_lh_edges(zero_VN_SC_lh_edges_index) = [];

SMN_SC_lh_edges(zero_SMN_SC_lh_edges_index) = [];
SMN_SC_lh_edges=zscore(SMN_SC_lh_edges,0,2);
SMN_FC_lh_edges(zero_SMN_SC_lh_edges_index) = [];

DAN_SC_lh_edges(zero_DAN_SC_lh_edges_index) = [];
DAN_SC_lh_edges=zscore(DAN_SC_lh_edges,0,2);
DAN_FC_lh_edges(zero_DAN_SC_lh_edges_index) = [];

VAN_SC_lh_edges(zero_VAN_SC_lh_edges_index) = [];
VAN_SC_lh_edges=zscore(VAN_SC_lh_edges,0,2);
VAN_FC_lh_edges(zero_VAN_SC_lh_edges_index) = [];

LN_SC_lh_edges(zero_LN_SC_lh_edges_index) = [];
LN_SC_lh_edges=zscore(LN_SC_lh_edges,0,2);
LN_FC_lh_edges(zero_LN_SC_lh_edges_index) = [];

FPN_SC_lh_edges(zero_FPN_SC_lh_edges_index) = [];
FPN_SC_lh_edges=zscore(FPN_SC_lh_edges,0,2);
FPN_FC_lh_edges(zero_FPN_SC_lh_edges_index) = [];

DMN_SC_lh_edges(zero_DMN_SC_lh_edges_index) = [];
DMN_SC_lh_edges=zscore(DMN_SC_lh_edges,0,2);
DMN_FC_lh_edges(zero_DMN_SC_lh_edges_index) = [];

SC_SC_lh_edges(zero_SC_SC_lh_edges_index) = [];
SC_SC_lh_edges=zscore(SC_SC_lh_edges,0,2);
SC_FC_lh_edges(zero_SC_SC_lh_edges_index) = [];

%rh
VN_SC_rh_edges(zero_VN_SC_rh_edges_index) = [];
VN_SC_rh_edges=zscore(VN_SC_rh_edges,0,2);
VN_FC_rh_edges(zero_VN_SC_rh_edges_index) = [];

SMN_SC_rh_edges(zero_SMN_SC_rh_edges_index) = [];
SMN_SC_rh_edges=zscore(SMN_SC_rh_edges,0,2);
SMN_FC_rh_edges(zero_SMN_SC_rh_edges_index) = [];

DAN_SC_rh_edges(zero_DAN_SC_rh_edges_index) = [];
DAN_SC_rh_edges=zscore(DAN_SC_rh_edges,0,2);
DAN_FC_rh_edges(zero_DAN_SC_rh_edges_index) = [];

VAN_SC_rh_edges(zero_VAN_SC_rh_edges_index) = [];
VAN_SC_rh_edges=zscore(VAN_SC_rh_edges,0,2);
VAN_FC_rh_edges(zero_VAN_SC_rh_edges_index) = [];

LN_SC_rh_edges(zero_LN_SC_rh_edges_index) = [];
LN_SC_rh_edges=zscore(LN_SC_rh_edges,0,2);
LN_FC_rh_edges(zero_LN_SC_rh_edges_index) = [];

FPN_SC_rh_edges(zero_FPN_SC_rh_edges_index) = [];
FPN_SC_rh_edges=zscore(FPN_SC_rh_edges,0,2);
FPN_FC_rh_edges(zero_FPN_SC_rh_edges_index) = [];

DMN_SC_rh_edges(zero_DMN_SC_rh_edges_index) = [];
DMN_SC_rh_edges=zscore(DMN_SC_rh_edges,0,2);
DMN_FC_rh_edges(zero_DMN_SC_rh_edges_index) = [];

SC_SC_rh_edges(zero_SC_SC_rh_edges_index) = [];
SC_SC_rh_edges=zscore(SC_SC_rh_edges,0,2);
SC_FC_rh_edges(zero_SC_SC_rh_edges_index) = [];

%相关系数计算
%lh
[r,p] = corr(VN_SC_lh_edges',VN_FC_lh_edges', 'type', CorrType);
VN_lh_Network_Coupling(i)=r;

[r,p] = corr(SMN_SC_lh_edges',SMN_FC_lh_edges', 'type', CorrType);
SMN_lh_Network_Coupling(i)=r;

[r,p] = corr(DAN_SC_lh_edges',DAN_FC_lh_edges', 'type', CorrType);
DAN_lh_Network_Coupling(i)=r;

[r,p] = corr(VAN_SC_lh_edges',VAN_FC_lh_edges', 'type', CorrType);
VAN_lh_Network_Coupling(i)=r;

[r,p] = corr(LN_SC_lh_edges',LN_FC_lh_edges', 'type', CorrType);
LN_lh_Network_Coupling(i)=r;

[r,p] = corr(FPN_SC_lh_edges',FPN_FC_lh_edges', 'type', CorrType);
FPN_lh_Network_Coupling(i)=r;

[r,p] = corr(DMN_SC_lh_edges',DMN_FC_lh_edges', 'type', CorrType);
DMN_lh_Network_Coupling(i)=r;

[r,p] = corr(SC_SC_lh_edges',SC_FC_lh_edges', 'type', CorrType);
SC_lh_Network_Coupling(i)=r;

%rh
[r,p] = corr(VN_SC_rh_edges',VN_FC_rh_edges', 'type', CorrType);
VN_rh_Network_Coupling(i)=r;

[r,p] = corr(SMN_SC_rh_edges',SMN_FC_rh_edges', 'type', CorrType);
SMN_rh_Network_Coupling(i)=r;

[r,p] = corr(DAN_SC_rh_edges',DAN_FC_rh_edges', 'type', CorrType);
DAN_rh_Network_Coupling(i)=r;

[r,p] = corr(VAN_SC_rh_edges',VAN_FC_rh_edges', 'type', CorrType);
VAN_rh_Network_Coupling(i)=r;

[r,p] = corr(LN_SC_rh_edges',LN_FC_rh_edges', 'type', CorrType);
LN_rh_Network_Coupling(i)=r;

[r,p] = corr(FPN_SC_rh_edges',FPN_FC_rh_edges', 'type', CorrType);
FPN_rh_Network_Coupling(i)=r;

[r,p] = corr(DMN_SC_rh_edges',DMN_FC_rh_edges', 'type', CorrType);
DMN_rh_Network_Coupling(i)=r;

[r,p] = corr(SC_SC_rh_edges',SC_FC_rh_edges', 'type', CorrType);
SC_rh_Network_Coupling(i)=r;


%左右脑合并的各网络相关系数计算
VN_SC_edges=[VN_SC_lh_edges,VN_SC_rh_edges];
VN_SC_edges=zscore(VN_SC_edges,0,2);
VN_FC_edges=[VN_FC_lh_edges,VN_FC_rh_edges];

SMN_SC_edges=[SMN_SC_lh_edges,SMN_SC_rh_edges];
SMN_SC_edges=zscore(SMN_SC_edges,0,2);
SMN_FC_edges=[SMN_FC_lh_edges,SMN_FC_rh_edges];

DAN_SC_edges=[DAN_SC_lh_edges,DAN_SC_rh_edges];
DAN_SC_edges=zscore(DAN_SC_edges,0,2);
DAN_FC_edges=[DAN_FC_lh_edges,DAN_FC_rh_edges];

VAN_SC_edges=[VAN_SC_lh_edges,VAN_SC_rh_edges];
VAN_SC_edges=zscore(VAN_SC_edges,0,2);
VAN_FC_edges=[VAN_FC_lh_edges,VAN_FC_rh_edges];

LN_SC_edges=[LN_SC_lh_edges,LN_SC_rh_edges];
LN_SC_edges=zscore(LN_SC_edges,0,2);
LN_FC_edges=[LN_FC_lh_edges,LN_FC_rh_edges];

FPN_SC_edges=[FPN_SC_lh_edges,FPN_SC_rh_edges];
FPN_SC_edges=zscore(FPN_SC_edges,0,2);
FPN_FC_edges=[FPN_FC_lh_edges,FPN_FC_rh_edges];

DMN_SC_edges=[DMN_SC_lh_edges,DMN_SC_rh_edges];
DMN_SC_edges=zscore(DMN_SC_edges,0,2);
DMN_FC_edges=[DMN_FC_lh_edges,DMN_FC_rh_edges];

SC_SC_edges=[SC_SC_lh_edges,SC_SC_rh_edges];
SC_SC_edges=zscore(SC_SC_edges,0,2);
SC_FC_edges=[SC_FC_lh_edges,SC_FC_rh_edges];

[r,p] = corr(VN_SC_edges',VN_FC_edges', 'type', CorrType);
VN_Network_Coupling(i)=r;

[r,p] = corr(SMN_SC_edges',SMN_FC_edges', 'type', CorrType);
SMN_Network_Coupling(i)=r;

[r,p] = corr(DAN_SC_edges',DAN_FC_edges', 'type', CorrType);
DAN_Network_Coupling(i)=r;

[r,p] = corr(VAN_SC_edges',VAN_FC_edges', 'type', CorrType);
VAN_Network_Coupling(i)=r;

[r,p] = corr(LN_SC_edges',LN_FC_edges', 'type', CorrType);
LN_Network_Coupling(i)=r;

[r,p] = corr(FPN_SC_edges',FPN_FC_edges', 'type', CorrType);
FPN_Network_Coupling(i)=r;

[r,p] = corr(DMN_SC_edges',DMN_FC_edges', 'type', CorrType);
DMN_Network_Coupling(i)=r;

[r,p] = corr( SC_SC_edges', SC_FC_edges', 'type', CorrType);
SC_Network_Coupling(i)=r;

 end

 %
Subnetwork_String_SplitLR=["VN_lh","SMN_lh","DAN_lh","VAN_lh","LN_lh","FPN_lh","DMN_lh","SC_lh","VN_rh","SMN_rh","DAN_rh","VAN_rh","LN_rh","FPN_rh","DMN_rh","SC_rh"];
CouplingResults_SplitLR=[VN_lh_Network_Coupling',SMN_lh_Network_Coupling', DAN_lh_Network_Coupling',VAN_lh_Network_Coupling',LN_lh_Network_Coupling',FPN_lh_Network_Coupling',DMN_lh_Network_Coupling',SC_lh_Network_Coupling',...
VN_rh_Network_Coupling',SMN_rh_Network_Coupling', DAN_rh_Network_Coupling',VAN_rh_Network_Coupling',LN_rh_Network_Coupling',FPN_rh_Network_Coupling',DMN_rh_Network_Coupling',SC_rh_Network_Coupling'];
Subnetwork_Level_Coupling_SplitLR={Subnetwork_String_SplitLR,CouplingResults_SplitLR};

Subnetwork_String=["VN","SMN","DAN","VAN","LN","FPN","DMN","SC"];
CouplingResults=[VN_Network_Coupling',SMN_Network_Coupling', DAN_Network_Coupling',VAN_Network_Coupling',LN_Network_Coupling',...
FPN_Network_Coupling',DMN_Network_Coupling',SC_Network_Coupling'];
Subnetwork_Level_Coupling={Subnetwork_String,CouplingResults};


% filename=fullfile(OutDir,'Subnetwork_Level_Coupling_SplitLR_BD.mat');
% save(filename,'Subnetwork_Level_Coupling_SplitLR');
% filename=fullfile(OutDir,'Subnetwork_Level_Coupling_IndependentEdgeVector_BD.mat');
% save(filename,'Subnetwork_Level_Coupling');


filename=fullfile(OutDir,'Subnetwork_Level_Coupling_SplitLR_HC.mat');
save(filename,'Subnetwork_Level_Coupling_SplitLR');
filename=fullfile(OutDir,'Subnetwork_Level_Coupling_IndependentEdgeVector_HC.mat');
save(filename,'Subnetwork_Level_Coupling');









