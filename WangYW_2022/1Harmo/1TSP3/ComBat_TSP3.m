function ComBat_TST(rawdata,outputpath,FeatureName)
if ischar(rawdata)
    raw = importdata(rawdata);
else
    raw=rawdata;
end

% for harmonization
batchnum =3;
subnum=[41,41,41];
batch = [];
%  batchnum  = input('sites number \n');

%     subnum = [];
for i = 1:batchnum
    %         subnum(1,i) = input('subject number for sites\n ')
    batch = [batch;kron(i,ones(subnum(1,i),1))];
end
Site = batch;

Age = xlsread('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/TST/TSTafterDpabi/info/subinfo.xlsx', 'F2:F42')
Age = [Age;Age;Age];
age = Age - mean(Age);

GENDER = xlsread('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/TST/TSTafterDpabi/info/subinfo.xlsx', 'D2:D42');
gender = dummyvar(GENDER);
gender = [gender;gender;gender];
mod = [age,gender(:,2)] % female vs male

fprintf('now is harmonizing %s \n',FeatureName );

para_adj_combat = combat(raw',Site, mod,1,0)';
nonpara_adj_combat= combat(raw',Site, mod,0,0)';
para_unadj_combat = combat(raw',Site, [],1,0)';
nonpara_unadj_combat = combat(raw',Site, [],0,0)';

save([outputpath,'/',FeatureName,'_para_adj_combat.mat'],'para_adj_combat');
save([outputpath,'/',FeatureName,'_nonpara_adj_combat.mat'],'nonpara_adj_combat');
save([outputpath,'/',FeatureName,'_para_unadj_combat.mat'],'para_unadj_combat');
save([outputpath,'/',FeatureName,'_nonpara_unadj_combat.mat'],'nonpara_unadj_combat');
end