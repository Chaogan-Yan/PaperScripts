

%Check for networks
load /mnt/Data/RfMRILab/Yan/YAN_Program/Atlas/Dos160_WithName.mat
Network=zeros(160,1);
for i=1:length(Dos160_WithName)
    switch Dos160_WithName{i,5}
        case 'occipital'
            Network(i)=1;
        case 'sensorimotor'
            Network(i)=2;
        case 'default'
            Network(i)=3;
        case 'fronto-parietal'
            Network(i)=4;
        case 'cingulo-opercular'
            Network(i)=5;
        case 'cerebellum'
            Network(i)=6;
    end
end

Network142=Network;
Network142(find(Network==6))=[];
Dos142_WithName=Dos160_WithName;
Dos142_WithName(find(Network==6),:)=[];

Network142_Yeo=Dos160_YeoNetwork_YanModified; %Network142_Yeo=Dos160_YeoNetwork;
Network142_Yeo(find(Network==6))=[];

Network142_Yeo(find(Network142_Yeo==10))=5; %Change subcortical to 5



DMNIndex=find(Network142_Yeo==7);






%Check for networks
load /mnt/Data/RfMRILab/Yan/YAN_Work/REST-meta-MDD/Processing/Stats/Stats_MDD_848_794/Network/Edge/Edge_LME.mat

%Restore to symetric
TMatrix=TMatrix+TMatrix';
TriuMat = triu(ones(size(PMatrix)),1);
PMatrix(find(TriuMat))=0;
PMatrix=PMatrix+PMatrix';
PMatrix(find(eye(size(PMatrix))))=1;

DMNTMatrix=TMatrix(DMNIndex,DMNIndex);
DMNPMatrix=PMatrix(DMNIndex,DMNIndex);

Dos142_WithName=Dos142_WithName(DMNIndex,:);




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

for iInd=1:length(PSurviveCountMatrixIndex)
    [II,JJ] = ind2sub(size(PSurviveCountMatrix),PSurviveCountMatrixIndex(iInd));
    if ~(JJ>=II)
        Row={Dos142_WithName{II,4},Dos142_WithName{JJ,4},II,JJ,DMNTMatrix(II,JJ),DMNPMatrix(II,JJ)};
        Table=[Table;Row];
    end
end



NodeDegree=sum(PSurviveCountMatrix);

%Create Node File
fid = fopen(['FullDMN.node'],'w');
for i=1:size(Dos142_WithName,1)
    Temp=Dos142_WithName{i,4};
    Index=strfind(Temp,' ');
    Temp(Index)='-';
    fprintf(fid,'%g\t%g\t%g\t1\t%g\t%s\n',Dos142_WithName{i,1},Dos142_WithName{i,2},Dos142_WithName{i,3},NodeDegree(i),Temp);
end
fclose(fid);

save(['FullDMN.edge'], 'PSurviveCountMatrix', '-ASCII', '-DOUBLE','-TABS')






