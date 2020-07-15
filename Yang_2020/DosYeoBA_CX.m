clear
load /mnt/Data/RfMRILab/Yan/YAN_Program/Atlas/Dos160_WithName.mat
BaFile = '/mnt/Data/share/Software/fMRI/DPABI_V2.3_170105/Templates/brodmann.nii';
[BaImg, Header] = y_Reslice(BaFile,'', 3, 0, '/mnt/Data/RfMRILab/Lile/DFC/CC200ROI_333.nii');
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

Dos160Yeo = Dos160_YeoNetwork_YanModified;
Dos160Yeo(Network ==6) = 8;
Dos160YeoName = cell(160,1);
DosInfo = cell(160,3);
for ii = 1:160
    coord = cell2mat(Dos160_WithName(ii, 1:3));
    ijk = Header.mat \ [coord(:); 1];
    ijk = round(ijk(1:3))';
%     disp([coord, ijk])
    Dos160BA(ii,1) = BaImg(ijk(1), ijk(2), ijk(3));
    
    switch Dos160Yeo(ii,1)
        case 1
            Dos160YeoName{ii}='VN';
        case 2
            Dos160YeoName{ii}='SMN';
        case 3
            Dos160YeoName{ii}='DAN';
        case 4
            Dos160YeoName{ii}='VAN';
        case 6
            Dos160YeoName{ii}='FPN';
        case 7
            Dos160YeoName{ii}='DMN';
        case 8
            Dos160YeoName{ii}='cerebellum';
        case 10
            Dos160YeoName{ii}='subcortical';
    end
    
    DosInfo{ii,1} = Dos160BA(ii,1);
    DosInfo{ii,2} = Dos160YeoName{ii,1};
    DosInfo{ii,3} = coord;
end

[~, NodeName] = xlsread('/mnt/Data/RfMRILab/ChenX/Small_World/Weighted_BugFree/SubgroupContrast_WeightedBugfree.xlsx', 2, 'z2:z500');
NodeInfo = cell(length(NodeName),3);
for ii = 1:length(NodeName)
    dn = str2double(NodeName{ii}(5:7));
    NodeInfo(ii,:) = DosInfo(dn, :);
end