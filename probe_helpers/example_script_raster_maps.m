addpath(genpath('C:/code/neuropixel-utils'))
addpath(genpath('C:/code/AlexA_Library'))
%load the processed data
session_name = 'npJ1_0520_baseline_2';
vr_data = load('Y:\giocomo\attialex\NP_DATA\npJ1_0520_baseline_2');

ks_dir = 'F:\J1\npJ1_0520_baseline_g0\npJ1_0520_baseline_g0_imec0';
channelMapFile = 'C:\code\KiloSort2\configFiles\neuropixPhase3B1_kilosortChanMap.mat';
ks = Neuropixel.KiloSortDataset(ks_dir,'channelMap',channelMapFile);
ks.load()
metrics = ks.computeMetrics();
image_save_dir = 'F:\images\npJ1_0520_baseline_2\rasters_Depth_Amplitude';
%%

rastersDepth_Amplitude(vr_data,metrics,image_save_dir)

final_image = combine_images('F:\images\npJ1_0520_baseline_2\rasters_Depth_Amplitude',8);
imwrite(final_image,['F:\images\combined\' session_name '_depth_amp.png'])