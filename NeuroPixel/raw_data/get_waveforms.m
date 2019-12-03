function [mean_waveforms,mean_temlate_waveforms] = get_waveforms(path)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
ks_dir = 'E:\1016_contrasttrack_gainchanges_1';
imec_dir = 'E:\npF2_L_1016_contrasttrack_gainchanges_1';
channelMapFile = 'E:\1016_contrasttrack_gainchanges_1\forPRBimecP3opt4_npF2_1016_contrasttrack_gainchanges_1.mat';
ks = Neuropixel.KiloSortDataset(ks_dir,'channelMap',channelMapFile,'imecDataset',imec_dir);
ks.load()
metrics = ks.computeMetrics();

snippetSetPre = ks.getWaveformsFromRawData('cluster_ids', clusterID,'num_waveforms', Inf, 'best_n_channels', 1, 'car', true,'spike_idx',idxPre);
meanwfPre=mean(snippetSetPre.data,3);
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

