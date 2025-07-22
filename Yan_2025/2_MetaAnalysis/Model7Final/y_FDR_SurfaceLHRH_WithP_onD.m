function [Data_CorrectedLH, Data_CorrectedRH, Header, P]=y_FDR_SurfaceLHRH_WithP_onD(StatsFileLH,StatsFileRH,MaskFileLH,MaskFileRH,qThreshold,FDRSuffix)



[BrainVolumeLH, VoxelSize, FileList, HeaderLH]=y_ReadAll(StatsFileLH);
[Path, Name, Ext]=fileparts(StatsFileLH);
[BrainVolumeLH_P]=y_ReadAll(fullfile(Path, [Name,'_P',Ext]));
[BrainVolumeLH_D]=y_ReadAll(fullfile(Path, [Name,'_D',Ext]));

[BrainVolumeRH, VoxelSize, FileList, HeaderRH]=y_ReadAll(StatsFileRH);
[Path, Name, Ext]=fileparts(StatsFileRH);
[BrainVolumeRH_P]=y_ReadAll(fullfile(Path, [Name,'_P',Ext]));
[BrainVolumeRH_D]=y_ReadAll(fullfile(Path, [Name,'_D',Ext]));


[nDimVertexLH,x]=size(BrainVolumeLH);
if ~isempty(MaskFileLH)
    [MaskDataLH]=y_ReadAll(MaskFileLH);
else
    MaskDataLH=ones(nDimVertexLH,1);
end
MaskDataOneDimLH=reshape(MaskDataLH,1,[]);
MaskIndexLH = find(MaskDataOneDimLH);


[nDimVertexRH,x]=size(BrainVolumeRH);
if ~isempty(MaskFileRH)
    [MaskDataRH]=y_ReadAll(MaskFileRH);
else
    MaskDataRH=ones(nDimVertexRH,1);
end
MaskDataOneDimRH=reshape(MaskDataRH,1,[]);
MaskIndexRH = find(MaskDataOneDimRH);




%FDR
SMapLH=BrainVolumeLH_P(MaskIndexLH);
SMapRH=BrainVolumeRH_P(MaskIndexRH);
PMap=[SMapLH;SMapRH];

% Following  FDR.m	1.3 Tom Nichols 02/01/18
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));

Thresholded=zeros(size(PMap));
if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
end


AllBrainLH = zeros(nDimVertexLH,1);
AllBrainLH(MaskIndexLH) = Thresholded(1:length(MaskIndexLH));

AllBrainRH = zeros(nDimVertexRH,1);
AllBrainRH(MaskIndexRH) = Thresholded(length(MaskIndexLH)+1:end);


Data_CorrectedLH=BrainVolumeLH_D.*AllBrainLH;
Data_CorrectedRH=BrainVolumeRH_D.*AllBrainRH;

[Path, Name, Ext]=fileparts(StatsFileLH);
y_Write(Data_CorrectedLH,HeaderLH,fullfile(Path, [Name,FDRSuffix,Ext]));

[Path, Name, Ext]=fileparts(StatsFileRH);
y_Write(Data_CorrectedRH,HeaderRH,fullfile(Path, [Name,FDRSuffix,Ext]));


fprintf('\n\tFDR correction finished.\n');
