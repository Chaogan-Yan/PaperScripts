clear;clc
parpool(20)
datapath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020FrameLine/2022MVPA/convergentwoexperienceMVPA/sortdataFour';
outpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020FrameLine/2022MVPA/convergentwoexperienceMVPA/MVPAresultR5';
mkdir(outpath)

maskpath='/mnt/Data3/RfMRILab/Lihuixian/DPABI_V6.1_220101/Templates/BrainMask_05_61x73x61.img';
subid=dir([datapath,'/sub*']);

parfor isub=1:size(subid,1)
    subpath=[datapath,'/',subid(isub).name];

    ds=cosmo_fmri_dataset(subpath,'mask',maskpath,...
        'targets',repmat([1;1;2;2],4,1),...  % 2condition
        'chunks',floor(((1:16)-1)/4)'+1); 
    % remove constant features
    ds=cosmo_remove_useless_data(ds);
    
    % Use a searchlight with a more-or-less constant number of voxels,
    % both near the edges of the brain and in the center of the brain.
    radius=5;
    nbrhood=cosmo_spherical_neighborhood(ds,'radius',radius);
    % Define the measure  as cosmo_crossvalidation_measure
    measure=@cosmo_crossvalidation_measure;
    
    %classifier: a function handle to cosmo_classify_lda
    %partitions: the output of cosmo_nfold_partitioner applied to the dataset
    args=struct();
    args.classifier=@cosmo_classify_lda;
    args.partitions=cosmo_nfold_partitioner(ds);
    
    %Run searchlight
    ds_cfy=cosmo_searchlight(ds,nbrhood,measure,args);
    
    output_fn=fullfile(outpath,[subid(isub).name(1:6),'lda_accuracy.nii']);
    cosmo_map2fmri(ds_cfy, output_fn);
end









