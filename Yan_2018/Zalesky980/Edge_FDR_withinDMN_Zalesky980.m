

load /mnt/Data/RfMRILab/Yan/YAN_Program/Atlas/Zalesky980_YeoNetwork.mat
Temp=load('/mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_848_794/ROIAnalysis/sgACC_Zalesky/FCMatrix/TMatrix_All.mat');

YeoNetwork(120)=7; %Put the 120th ROI, sgACC, to be DMN.

ROI_Remove=Temp.NaNROI(1:end-1); %We still need the 120th ROI
ROIWanted=1:980;
ROIWanted(ROI_Remove)=[];
YeoNetwork(ROI_Remove)=[];
ROIWanted(find(YeoNetwork==0))=[];
YeoNetwork(find(YeoNetwork==0))=[];



DMNIndex=find(YeoNetwork==7);



load /mnt/Data/RfMRILab/Yan/YAN_Program/Atlas/Zalesky980_Center.mat
Zalesky980_Center=Zalesky980_Center(ROIWanted,:);
Zalesky980_Center=Zalesky980_Center(DMNIndex,:);



%Check for networks
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_848_794/Network/Edge/Zalesky980/All.mat

%Restore to symetric
TMatrix=TMatrix+TMatrix';
TriuMat = triu(ones(size(PMatrix)),1);
PMatrix(find(TriuMat))=0;
PMatrix=PMatrix+PMatrix';
PMatrix(find(eye(size(PMatrix))))=1;

DMNTMatrix=TMatrix(DMNIndex,DMNIndex);
DMNPMatrix=PMatrix(DMNIndex,DMNIndex);




%FDR
TriuMat = triu(ones(size(DMNPMatrix)),1)';

PVector = DMNPMatrix(find(TriuMat));

addpath /mnt/Data/RfMRILab/Yan/YAN_Program/gretna

[pID,pN] = FDR(PVector,0.05);

PSig = DMNPMatrix<=pID;









PSurviveP = (DMNPMatrix<=pID).*(DMNTMatrix>0);
PSurviveN = (DMNPMatrix<=pID).*(DMNTMatrix<0);


%Get connection name
Table=[];

PSurviveCountMatrix = PSurviveN;
PSurviveCountMatrixIndex = find(PSurviveCountMatrix);

% for iInd=1:length(PSurviveCountMatrixIndex)
%     [II,JJ] = ind2sub(size(PSurviveCountMatrix),PSurviveCountMatrixIndex(iInd));
%     if ~(JJ>=II)
%         Row={Dos142_WithName{II,4},Dos142_WithName{JJ,4},II,JJ,DMNTMatrix(II,JJ),DMNPMatrix(II,JJ)};
%         Table=[Table;Row];
%     end
% end



NodeDegree=sum(PSurviveCountMatrix);

%Create Node File
fid = fopen(['FullDMN.node'],'w');
for i=1:size(Zalesky980_Center,1)
    Temp='-';
    fprintf(fid,'%g\t%g\t%g\t1\t%g\t%s\n',round(Zalesky980_Center(i,6)),round(Zalesky980_Center(i,7)),round(Zalesky980_Center(i,8)),NodeDegree(i),Temp);
end
fclose(fid);

save(['FullDMN.edge'], 'PSurviveCountMatrix', '-ASCII', '-DOUBLE','-TABS')






