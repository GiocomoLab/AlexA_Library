addpath(genpath('C:/code/neuropixel-utils'))
addpath(genpath('C:/code/AlexA_Library'))
%load the processed data
session_name = 'npJ1_0520_baseline_2';
load('Y:\giocomo\attialex\NP_DATA\npJ1_0520_baseline_2')

ks_dir = 'E:\1016_contrasttrack_gainchanges_1';
imec_dir = 'E:\npF2_L_1016_contrasttrack_gainchanges_1';
channelMapFile = 'E:\1016_contrasttrack_gainchanges_1\forPRBimecP3opt4_npF2_1016_contrasttrack_gainchanges_1.mat';
ks = Neuropixel.KiloSortDataset(ks_dir,'channelMap',channelMapFile,'imecDataset',imec_dir);
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
