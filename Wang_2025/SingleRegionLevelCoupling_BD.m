%This script is used to caculate the single region level SC-FC coupling
%Each subject will be given 454 values representing the degree of their brain SC-FC coupling in each region defined by Schaefer2018_400_Tian2020_54_atlas
clc;clear;


% %Load the BD data file
% load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_BD\FC_edges.mat');
% load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_BD\SC_edges.mat');
% OutDir = 'D:\BD_SCFC_Coupling\SensitivityAnalysis\Results\C_ResultsSingleRegionLevel\BD\';

% Load the HC data file
load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_HC\FC_edges.mat');
load('D:\BD_SCFC_Coupling\All_Significant_Files_Final\MatrixToMat_HC\SC_edges.mat');
OutDir = 'D:\BD_SCFC_Coupling\SensitivityAnalysis\Results\C_ResultsSingleRegionLevel\HC\';

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

%Caculate the single region level SC-FC coupling
for i=1:nsub
    for j=1:nreg
        SC_matrix_i=SC_matrix(:,:,i);
        FC_matrix_i=FC_matrix(:,:,i);

        SC_matrix_i_j=SC_matrix_i(j,:);
        FC_matrix_i_j=FC_matrix_i(j,:);

        zero_SC_matrix_i_j_index=find(SC_matrix_i_j==0);

        SC_matrix_i_j(zero_SC_matrix_i_j_index)=[];
        FC_matrix_i_j(zero_SC_matrix_i_j_index )=[];

       SC_matrix_i_j=zscore( SC_matrix_i_j,0,2);

        [r,p] = corr( SC_matrix_i_j',FC_matrix_i_j', 'type', 'Pearson');
        SingleRegionLevel_Coupling(i,j)=r;
    end
end

SingleRegionLevel_Coupling_Cortical=SingleRegionLevel_Coupling(:,1:400);
SingleRegionLevel_Coupling_SubCortical=SingleRegionLevel_Coupling(:,401:454);

% filename=fullfile(OutDir,'SingleRegionLevel_Coupling_Cortical_BD.mat');
% save(filename,'SingleRegionLevel_Coupling_Cortical');
% filename=fullfile(OutDir,'SingleRegionLevel_Coupling_SubCortical_BD.mat');
% save(filename,'SingleRegionLevel_Coupling_SubCortical');
% filename=fullfile(OutDir,'SingleRegionLevel_Coupling_Whole_BD.mat');
% save(filename,'SingleRegionLevel_Coupling');

filename=fullfile(OutDir,'SingleRegionLevel_Coupling_Cortical_HC.mat');
save(filename,'SingleRegionLevel_Coupling_Cortical');
filename=fullfile(OutDir,'SingleRegionLevel_Coupling_SubCortical_HC.mat');
save(filename,'SingleRegionLevel_Coupling_SubCortical');
filename=fullfile(OutDir,'SingleRegionLevel_Coupling_Whole_HC.mat');
save(filename,'SingleRegionLevel_Coupling');










