%This script is used to caculate the global level SFC
%Each subject will be given a value representing the degree of their whole brain SFC

clc;clear;

%Load the BD data file

% load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_BD\FC_edges.mat');
% load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_BD\SC_edges.mat');
% OutDir = 'D:\BD_SCFC_Coupling\SensitivityAnalysis\Results\GlobalLevelCoupling\BD\';

% Load the HC data file
load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_HC\FC_edges.mat');
load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_HC\SC_edges.mat');
OutDir = 'D:\BD_SCFC_Coupling\SensitivityAnalysis\Results\GlobalLevelCoupling\HC\';

nsub=size(FC_edges,1);

% Caculate the global level coupling
for i = 1:nsub
    SC_edges_i=SC_edges(i,:);
    FC_edges_i=FC_edges(i,:);
    zero_SC_idx = find(SC_edges_i == 0);

    SC_edges_i(zero_SC_idx ) = [];
    FC_edges_i(zero_SC_idx ) = [];


    s_SC_edges_i=zscore( SC_edges_i,0,2);

    g1_SC_edges_i=s_SC_edges_i*0.1+0.5;

    g2_SC_edges_i=s_SC_edges_i*0.2+0.25;

    [r,p] = corr(SC_edges_i',FC_edges_i', 'type', 'Spearman');
    GlobalLevel_Coupling(i)=r;

    [r,p] = corr( s_SC_edges_i',FC_edges_i', 'type', 'Pearson');
    GlobalLevel_Coupling1(i)=r;

        [r,p] = corr( g1_SC_edges_i',FC_edges_i', 'type', 'Pearson');
    GlobalLevel_Coupling2(i)=r;

            [r,p] = corr( g2_SC_edges_i',FC_edges_i', 'type', 'Pearson');
    GlobalLevel_Coupling3(i)=r;



end

% filename=fullfile(OutDir,'GlobalLevel_Coupling_BD.mat');
% save(filename,'GlobalLevel_Coupling');
filename=fullfile(OutDir,'GlobalLevel_Coupling_HC.mat');
save(filename,'GlobalLevel_Coupling');

