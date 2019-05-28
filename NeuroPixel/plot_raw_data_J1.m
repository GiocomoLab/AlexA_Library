addpath(genpath('C:/code/neuropixel-utils'))
addpath(genpath('C:/code/AlexA_Library'))
session_name = 'npJ1_0520_baseline_2';
load('Y:\giocomo\attialex\NP_DATA\npJ1_0520_baseline_2')
ks_dir = 'F:\J1\npJ1_0520_baseline_g0\npJ1_0520_baseline_g0_imec0';
channelMapFile = 'C:\code\KiloSort2\configFiles\neuropixPhase3B1_kilosortChanMap.mat';
ks = Neuropixel.KiloSortDataset(ks_dir,'channelMap',channelMapFile);
ks.load()
metrics = ks.computeMetrics();


%%
verify_spatialMaps;


%%
trial_pre=64:75;
trial_post = 75:79;
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
metrics.plotClusterDriftmap('cluster_ids', ks.clusters_good([51]));
hold on
for ii=1:49; plot([trial_onsets_real(ii) trial_onsets_real(ii)],[0 2000],'r');end
%xlim([5396 5396+max(sp.st)])
