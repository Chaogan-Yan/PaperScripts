%Read the multiband data's sliceorder and save it as .txt files
%Read resting-state data
clc;clear;
Workdir='/mnt/Data/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/PKUSIMENS_rest';
Sublist = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Scripts/Analysis/PKUSIMENS_Sublist.txt');

AllSliceOrder=[];
for i = 1:length(Sublist)
   Datadir=dir([Workdir,'/FunRaw/',Sublist{i},'/*.dcm']);
   DataName=Datadir.name;
   DicomData=dicominfo([Workdir,'/FunRaw/',Sublist{i},'/',DataName]);
   SliceOrder=DicomData.Private_0019_1029;
   
   save([Workdir,'/',Sublist{i},'_SliceOrder.txt'],'SliceOrder','-ascii');
   
   AllSliceOrder=[AllSliceOrder,SliceOrder];
end

%Read task-state data
clc;clear;
Workdir='/mnt/Data/RfMRILab/ChenX/Rumination_project/Data/Full_Preprocessing/PKUSIEMENS_task4Pub_oldVersion';
Sublist = importdata('/mnt/Data/RfMRILab/ChenX/Rumination_project/Analysis/Analysis_majorRevision/Analysis4Pub/PKUSIEMENS_Sublist.txt');

%Dealing with the first session, sad
for i = 1:length(Sublist)
   Datadir=dir([Workdir,'/FunRaw/',Sublist{i},'/*.dcm']);
   DataName=Datadir(1).name;
   DicomData=dicominfo([Workdir,'/FunRaw/',Sublist{i},'/',DataName]);
   SliceOrder=DicomData.Private_0019_1029;
   
   save([Workdir,'/',Sublist{i},'_S1_SliceOrder.txt'],'SliceOrder','-ascii');
   
end

%Dealing with the second session, rum
for i = 1:length(Sublist)
   Datadir=dir([Workdir,'/S2_FunRaw/',Sublist{i},'/*.dcm']);
   DataName=Datadir(1).name;
   DicomData=dicominfo([Workdir,'/S2_FunRaw/',Sublist{i},'/',DataName]);
   SliceOrder=DicomData.Private_0019_1029;
   
   save([Workdir,'/',Sublist{i},'_S2_SliceOrder.txt'],'SliceOrder','-ascii');
   
end

%Dealing with the third session, dis
for i = 1:length(Sublist)
   Datadir=dir([Workdir,'/S3_FunRaw/',Sublist{i},'/*.dcm']);
   DataName=Datadir(1).name;
   DicomData=dicominfo([Workdir,'/S3_FunRaw/',Sublist{i},'/',DataName]);
   SliceOrder=DicomData.Private_0019_1029;
   
   save([Workdir,'/',Sublist{i},'_S3_SliceOrder.txt'],'SliceOrder','-ascii');
   
end