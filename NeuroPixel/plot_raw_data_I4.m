session_name = 'npI4_0414_baseline_2';
load('Y:\giocomo\attialex\NP_DATA\npI4_0424_baseline_2.mat')
ks_dir = 'Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I4\neuropixels_data\npI4_0424_gaincontrast_g0\npI4_0424_gaincontrast_g0_imec0';
channelMapFile = 'C:\code\KiloSort2\configFiles\neuropixPhase3B1_kilosortChanMap.mat';
ks = Neuropixel.KiloSortDataset(ks_dir,'channelMap',channelMapFile);
ks.load()
ks.computeMetrics()

image_save_dir = strcat('Y:\giocomo\attialex\');
image_save_dir = fullfile(image_save_dir,'images',session_name,'rasters_waveforms');
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end
%%
verify_spatialMaps;


trial_pre=42:62;
trial_post = 62:82;

extract_plot_raw_data;
raw_data_similarity;
