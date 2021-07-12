ks_dir = 'Z:\giocomo\export\data\Projects\AlexA_NP\AA_190709_5\neuropixels_data\AA_190709_5_0808_gaincontrast_g0\AA_190709_5_0808_gaincontrast_g0_imec0';
[parent_dir,~,~] =fileparts(ks_dir); 
NIDAQ_file = dir(fullfile(parent_dir,'*nidq.bin'));
NIDAQ_file = fullfile(parent_dir,NIDAQ_file(1).name);

NIDAQ_config = dir(fullfile(parent_dir,'*nidq.meta'));
NIDAQ_config = fullfile(parent_dir,NIDAQ_config(1).name);

animal = 'AA_190709_5'
vr_file_full = {'Z:\giocomo\export\data\Projects\AlexA_NP\AA_190709_5\neuropixels_data\AA_190709_5_0808_gaincontrast_g0\MOV_190808_12-01-19.log'};
%%

sync_vrPanda_to_np(ks_dir,NIDAQ_file,NIDAQ_config,vr_file_full,animal,'F:\Alex\new_2','')