function [ResultBrain, TimeSeries, Target_Coordinate] = TMS_Targeting_DR_V4(AllVolume, SeedMap, OutputName, MaskData, PeakMode, KeepRatio, nNeighbor, SaveImage, Header, Verbose)
% Using dual regression (DR) to calculate individualized TMS targets based on resting-state fMRI and group-level templates.
% Please refer to https://www.biorxiv.org/content/10.1101/2023.03.09.531726v3.full
%
% Input:
% 	AllVolume		-	4D data matrix (DimX*DimY*DimZ*DimTimePoints) or the directory of 3D image data file or the filename of one 4D data file.
%   SeedMap        -   A spatial seed map for dual regression (i.e., 3D mask martrix (DimX*DimY*DimZ)) or the directory of it.
%	OutputName  	-	Output filename.
% 	MaskData		-   The Brain Mask matrix (DimX*DimY*DimZ) or the Brain Mask file name
%   PeakMode       -   The method of selecting peak points from DR derived maps as the final TMS target coordinates.
%                     -1. 'Peak': The peak voxel in the search area (i.e., mask) of DR derived maps
%                     -2. 'Size': The spatial centroid of the largest cluster after thresholding.
%                     -3. 'WeightedSize': The spatial centroid of the largest cluster after thresholding, and using the values from DR-derived maps as weighting factors.
%                     -4. 'CoM': Center of mass (i.e., spatial centroid) of all voxels after thresholding.
%   KeepRatio       -   Thresholding. Keep a certain number of voxels within the search area (i.e., mask) of the DR-derived maps.
%                     -1. if KeepRatio<1 (i.e. 0.1), retain a proportion of voxels at the KeepRatio (i.e.10%) level.
%                     -2. if KeepRatio>1 (i.e. 10), retain a certain number of voxels (i.e. 10 voxels).
%   nNeighbor      -   Number of voxel neighborhoods used to define a cluster, as specified by the "Pixel Connectivity" argument in the function "bwlabeln".
%   SaveImage      -   Save a NIFTI image for the derived target.
%   Header           - A NIFTI header file. Necessary when the input for the AllVolume argument is a data matrix instead of a filename.
%   Verbose         - A boolean argument that controls whether to display detailed information during execution.
%
% Output:
%	ResultBrain     -   The dual regression map
%   TimeSeries      -   The time series after first regression
%   Header          -   The NIfTI Header
%-----------------------------------------------------------
% Written by Bin Lu 2023-01-01.
% Institute of Psychology, Chinese Academy of Sciences, 16 Lincui Road, Chaoyang District, Beijing 100101, China
% larslu@foxmail.com


if Verbose
    fprintf('\n\t Performing Dual Regression...');
end

if ~isnumeric(AllVolume)
    if Verbose
        [AllVolume,VoxelSize,theImgFileList, Header] =y_ReadAll(AllVolume);
    else
        evalc('[AllVolume,VoxelSize,theImgFileList, Header] =y_ReadAll(AllVolume);');
    end
end

[nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
BrainSize = [nDim1 nDim2 nDim3];

if ischar(MaskData)
    if Verbose
        fprintf('\nLoad mask "%s".\n', MaskData);
    end
    if ~isempty(MaskData)
        [MaskData,MaskHead]=y_Read(MaskData);
        if ~all(size(MaskData)==[nDim1 nDim2 nDim3])
            error('The size of Mask (%dx%dx%d) doesn''t match the required size (%dx%dx%d).\n',size(MaskData), [nDim1 nDim2 nDim3]);
        end
        MaskData = double(logical(MaskData));
    else
        MaskData=ones(nDim1,nDim2,nDim3);
    end
end

%Get seed
if ~isnumeric(SeedMap)
    [SeedMap,~] = y_Read(SeedMap);
end
SeedMap=reshape(SeedMap,[],1);

% Convert into 2D
AllVolume=reshape(AllVolume,[],nDimTimePoints)';
MaskDataOneDim=reshape(MaskData,1,[]);
% MaskIndex = find(MaskDataOneDim' | SeedMap);
MaskIndex = find(MaskDataOneDim);
AllVolume=AllVolume(:,MaskIndex);
SeedMap=SeedMap(MaskIndex); 

%Frist regression
Cov = SeedMap;
Cov = (Cov - mean(Cov))/std(Cov); %Mean centering and variance normalization
SpatialMean = mean(AllVolume,2);
AllVolumeDeSpatialMean = (AllVolume-repmat(SpatialMean,1,size(AllVolume,2)));
TimeSeries=AllVolumeDeSpatialMean*Cov*inv(Cov'*Cov);

%Second regression
Cov=TimeSeries;
Cov = (Cov - mean(Cov))/std(Cov); %Mean centering and variance normalization
TemporalMean = mean(AllVolume);
AllVolumeDeTemporalMean = AllVolume-repmat(TemporalMean,size(AllVolume,1),1);
b=inv(Cov'*Cov)*Cov'*AllVolume;

switch PeakMode
    case 'Peak' %% directly using peak point
        ResultBrain=zeros(size(MaskDataOneDim));
        ResultBrain(1,MaskIndex)=b;
        ResultBrain=reshape(ResultBrain,nDim1, nDim2, nDim3);
        [I,J,K] = ind2sub(size(ResultBrain),find(ResultBrain==max(b)));
        Target_Coordinate = cor2mni([I,J,K], Header.mat);
    case 'Size' %% find the largest cluster first and calculated the center of mass
        % Thresholding
        ResultBrain=zeros(size(MaskDataOneDim));
        Index = find(b);
        nVoxel = length(Index);
        [~,I] = sort(b(Index));
        if KeepRatio < 1
            b(Index(I(1:end-round(nVoxel*KeepRatio)))) = 0;
        else
            b(Index(I(1:end-KeepRatio))) = 0;
        end
        ResultBrain(1,MaskIndex)=b;
        ResultBrain=reshape(ResultBrain,nDim1, nDim2, nDim3);
        
        % Only keep the largest cluster
        [cci,num] = bwlabeln(ResultBrain,nNeighbor);
        ClusterValue = arrayfun(@(Index) length(find(cci==Index)), unique(cci));
        [~,Largest_Index] = max(ClusterValue(2:end));
        ResultBrain(cci~=Largest_Index) = 0;
        
        % Calculate center of mass
        [I,J,K] = ind2sub(size(ResultBrain),find(ResultBrain));
        Weight = ResultBrain(find(ResultBrain));   % use FC as weight
        Coordinates = cor2mni([I,J,K], Header.mat);
        Coordinates_Weighted = Coordinates.*repmat(Weight,1,3);
        Target_Coordinate = sum(Coordinates_Weighted,1)./sum(Weight);
    case 'WeightedSize' %% find the largest cluster first and calculated the center of mass
        % Thresholding
        ResultBrain=zeros(size(MaskDataOneDim));
        Index = find(b);
        nVoxel = length(Index);
        [~,I] = sort(b(Index));
        if KeepRatio < 1
            b(Index(I(1:end-round(nVoxel*KeepRatio)))) = 0;
        else
            b(Index(I(1:end-KeepRatio))) = 0;
        end
        ResultBrain(1,MaskIndex)=b;
        ResultBrain=reshape(ResultBrain,nDim1, nDim2, nDim3);
        
        % Only keep the largest cluster - cluster size weighted by voxel value
        [cci,num] = bwlabeln(ResultBrain,nNeighbor);
        ClusterSizes = zeros(num,1);
        for i = 1:num
            ClusterSizes(i) = sum(ResultBrain(cci == i), 'all');
        end
        [~,Largest_Index] = max(ClusterSizes);
        ResultBrain(cci~=Largest_Index) = 0;
        
        % Calculate center of mass
        [I,J,K] = ind2sub(size(ResultBrain),find(ResultBrain));
        Weight = ResultBrain(find(ResultBrain));   % use FC as weight
        Coordinates = cor2mni([I,J,K], Header.mat);
        Coordinates_Weighted = Coordinates.*repmat(Weight,1,3);
        Target_Coordinate = sum(Coordinates_Weighted,1)./sum(Weight);
    case 'CoM' %%  center of mass of all the survival voxels
        % Thresholding
        ResultBrain=zeros(size(MaskDataOneDim));
        Index = find(b);
        nVoxel = length(Index);
        [~,I] = sort(b(Index));
        if KeepRatio < 1
            b(Index(I(1:end-round(nVoxel*KeepRatio)))) = 0;
        else
            b(Index(I(1:end-KeepRatio))) = 0;
        end
        ResultBrain(1,MaskIndex)=b;
        ResultBrain=reshape(ResultBrain,nDim1, nDim2, nDim3);
        
        % Calculate center of mass
        [I,J,K] = ind2sub(size(ResultBrain),find(ResultBrain));
        Weight = ResultBrain(find(ResultBrain));   % use FC as weight
        Coordinates = cor2mni([I,J,K], Header.mat);
        Coordinates_Weighted = Coordinates.*repmat(Weight,1,3);
        Target_Coordinate = sum(Coordinates_Weighted,1)./sum(Weight);
end

if SaveImage
    Header.pinfo = [1;0;0];
    Header.dt    =[16,0];
    y_Write(ResultBrain,Header,OutputName);
end

if Verbose
    fprintf('\n\t Dual Regression based TMS targets optimization finished\n\t ');
end



function mni = cor2mni(cor, T)
% function mni = cor2mni(cor, T)
% convert matrix coordinate to mni coordinate
%
% cor: an Nx3 matrix
% T: (optional) rotation matrix
% mni is the returned coordinate in mni space
%
% caution: if T is not given, the default T is
% T = ...
%     [-4     0     0    84;...
%      0     4     0  -116;...
%      0     0     4   -56;...
%      0     0     0     1];
%
% xu cui
% 2004-8-18
% last revised: 2005-04-30

if nargin == 1
    T = ...
        [-4     0     0    84;...
         0     4     0  -116;...
         0     0     4   -56;...
         0     0     0     1];
end

cor = round(cor);
mni = T*[cor(:,1) cor(:,2) cor(:,3) ones(size(cor,1),1)]';
mni = mni';
mni(:,4) = [];
return;

 
