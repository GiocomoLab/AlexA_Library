session_name = 'npI4_0414_baseline_2';
load('Y:\giocomo\attialex\NP_DATA\npI4_0424_baseline_2.mat')
ks_dir = 'I:\npI4_0424_gaincontrast_g0\npI4_0424_gaincontrast_g0_imec0';
channelMapFile = 'C:\code\KiloSort2\configFiles\neuropixPhase3B1_kilosortChanMap.mat';
ks = Neuropixel.KiloSortDataset(ks_dir,'channelMap',channelMapFile);
ks.load()
metrics = ks.computeMetrics();


%%
verify_spatialMaps;

%%
trial_pre=37:56;
trial_post = 56:75;
image_save_dir = strcat('Y:\giocomo\attialex\');
image_save_dir = fullfile(image_save_dir,'images',session_name,sprintf('rasters_waveforms_%d_%d',trial_pre(1),trial_post(end)));
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end

extract_plot_raw_data;

image_save_dir = strcat('Y:\giocomo\attialex\');
image_save_dir = fullfile(image_save_dir,'images',session_name,sprintf('rasters_waveforms_similarity_%d_%d',trial_pre(1),trial_post(end)));
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end
raw_data_similarity;
%%
figure
metrics.plotDriftmap()  

figure
cluster_ids = metrics.cluster_ids(1:5:metrics.nClusters);
metrics.plotClusterDriftmap('cluster_ids', cluster_ids);

figure
metrics.plotClusterDriftmap('cluster_ids', [486 479 469 461 457]);
xlim([5396 5396+max(sp.st)])
