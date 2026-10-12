% Example configuration for ana_DPABIViewerFixed.m.
configDir = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(configDir);
addpath(repositoryRoot);

cfg = struct();
cfg.dataRoot = getenv('MDD_FC_RESULTS_ROOT');
cfg.templateDir = getenv('DPABI_TEMPLATE_DIR');
cfg.seedManifest = fullfile(configDir, 'fc_seed_manifest.csv');
cfg.maskFile = fullfile(cfg.templateDir, 'BrainMask_05_97x115x97.nii');
cfg.expectedSeedCount = 4;
cfg.voxelP = 0.001;
cfg.clusterP = 0.05;
cfg.isTwoTailed = true;
cfg.connectivity = 26;
cfg.PositiveDirectionLabel = 'MDD>HC';
cfg.NegativeDirectionLabel = 'MDD<HC';
cfg.SMOOTHNESS_SOURCE = 'header_strict';
% Optional overrides when metadata are absent from the NIfTI header:
cfg.MANUAL_DF1 = [];
cfg.MANUAL_DLH = [];
cfg.MANUAL_FWHM = [];

ana_DPABIViewerFixed(cfg);
