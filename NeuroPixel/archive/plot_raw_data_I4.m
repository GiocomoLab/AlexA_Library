addpath(genpath('C:/code/neuropixel-utils'))
addpath(genpath('C:/code/AlexA_Library'))
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
xlim([sp.vr_session_offset 5396+max(sp.st)])
l1=min(post(trial==trial_pre(1)))+sp.vr_session_offset;
xline(l1,'k')
l1=min(post(trial==trial_post(1)))+sp.vr_session_offset;
xline(l1,'k')
l1=min(post(trial==trial_post(end)))+sp.vr_session_offset;
xline(l1,'k')

figure
metrics.plotClusterDriftmap('cluster_ids', [486 479 469 461 457]);
xlim([5396 5396+max(sp.st)])

%%
figure
metrics.plotDriftmap()
ylim([0 2500])
pause
xlim([sp.vr_session_offset 5396+max(sp.st)])
pause
l1=min(post(trial==trial_pre(1)))+sp.vr_session_offset;
xline(l1,'g')
l1=min(post(trial==trial_post(1)))+sp.vr_session_offset;
xline(l1,'g')
l1=min(post(trial==trial_post(end)))+sp.vr_session_offset;
xline(l1,'g')



%%
figure
metrics.plotClusterDriftmap('cluster_ids', [448 585]);
xlim([sp.vr_session_offset 5396+max(sp.st)])
l1=min(post(trial==trial_pre(1)))+sp.vr_session_offset;
xline(l1,'k')
l1=min(post(trial==trial_post(1)))+sp.vr_session_offset;
xline(l1,'k')
l1=min(post(trial==trial_post(end)))+sp.vr_session_offset;
xline(l1,'k')

%%
for cid = [554 550 496 448 683 585]
    h = figure('Position',[100 100 380 500]); hold on;
    
    idx = find(ks.clusters_good==cid);
    spike_id=sp.clu==cid;
    spike_t = sp.st(spike_id);
    [~,~,spike_idx] = histcounts(spike_t,post);
    
   
    subplot(2,2,[1 3])
    scatter(posx(spike_idx),trial(spike_idx),2)
    title(sprintf('cluid = %d',cid))
    subplot(2,2,2)
    imagesc(squeeze(correlation_All(:,:,idx)),[-.1 .8])
    axis image
end


