clear;clc;
%%load info
%corr
load /mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/SubInfo/SubInfo_420.mat ;
%load /mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/SubInfo/map_402in420.mat ;
Sex(Sex==-1)=0; %female
%fcp
FCP.info = importdata('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/FCP_Organized/SubInfo/600_subinfo.mat');
FCP.info.Sex(FCP.info.Sex==1)=0 %female
FCP.info.Sex(FCP.info.Sex==2)=1 %male
Age = [Age;FCP.info.Age];
Sex = [Sex;FCP.info.Sex];
Motion = [Motion;[FCP.info.MeanFDJ,FCP.info.MeanFDJ]];
SubID = [SubID;FCP.info.Subid];
Site = [Site;FCP.info.SiteID+max(Site)];


if length(unique(Site))~=max(Site)
    s = unique(Site);
    for i = 1:length(s)
        Site(Site==s(i)) =  i ;
    end
end

clear FCP;
MaskFile ='/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/Mask/Overlap38810.nii';
Mask = y_Read(MaskFile);
mask_size = size(Mask);
% mapping index
x =  find(reshape(Mask,[],1)==1);
newmap = zeros(size(reshape(Mask,[],1)));
newmap(x,1) = 1;

IndexName = {'ReHo_FunImgARCWF','ALFF_FunImgARCW','fALFF_FunImgARCW','DegreeCentrality_FunImgARCWF','FC_D142'};
outputpath = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/FCP';

% load data
%parpool(4);
Rdir ={'Results','S2_Results'}
for ses = 1:2
    mod = [double(Age),Sex,Motion(:,ses)];
    for i =1:length(IndexName)        
        com = importdata(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/',Rdir{ses},'/',IndexName{i},'_raw.mat']);
        fcp = importdata(['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/FCP/',Rdir{ses},'/',IndexName{i},'_raw.mat']);      
        raw = [com;fcp];
        
        % clear fcp;
        % clear combat;
        
        fprintf('now is harmonizing %s \n',IndexName{i} );
        outputpath = ['/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/FCP/',Rdir{ses}];

        para_adj(raw,Site, mod, outputpath, IndexName{i});
        nonpara_adj(raw,Site, mod, outputpath, IndexName{i});
        para_unadj(raw,Site, outputpath, IndexName{i});
        nonpara_unadj(raw,Site, outputpath, IndexName{i});
        % para_unadj_combat = combat(raw',Site, [],1,0)';
        % save([outputpath,'/',IndexName{i} ,'_para_unadj_combat.mat'],'para_unadj_combat');
        %
        % nonpara_unadj_combat = combat(raw',Site, [],0,0)';
        % save([outputpath,'/',IndexName{i} ,'_nonpara_unadj_combat.mat'],'nonpara_unadj_combat');
        
    end
end
function para_adj(data,Site,mod,outputpath,featurename)
para_adj_combat = combat(data',Site, mod,1,0)';
save([outputpath,'/',featurename ,'_para_adj_combat.mat'],'para_adj_combat');
end
function nonpara_adj(data,Site,mod,outputpath,featurename)
nonpara_adj_combat =combat(data',Site, mod,0,0)';;
save([outputpath,'/',featurename ,'_nonpara_adj_combat.mat'],'nonpara_adj_combat');
end
function para_unadj(data,Site,outputpath,featurename)
para_unadj_combat = combat(data',Site,[],1,0)';
save([outputpath,'/',featurename ,'_para_unadj_combat.mat'],'para_unadj_combat');
end
function nonpara_unadj(data,Site,outputpath,featurename)
nonpara_unadj_combat =combat(data',Site, [],0,0)';;
save([outputpath,'/',featurename ,'_nonpara_unadj_combat.mat'],'nonpara_unadj_combat');
end
