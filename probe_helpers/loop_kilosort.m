addpath(genpath('C:\code\KiloSort2')) % path to kilosort folder
addpath('C:\code\npy-matlab')
addpath(genpath('C:\code\AlexA_Library'))
addpath('C:\code\SiliconProbeCode\')

% directories={'E:\npI5_0418_playback_1_g0\npI5_0418_playback_1_g0_imec0',...
%             'E:\npI1_0418_mismatch_g0\npI1_0418_mismatch_g0_imec0'};
directories={'Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I4\neuropixels_data\npI4_0420_baseline_g0\npI4_0420_baseline_g0_imec0',...
    'Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I4\neuropixels_data\npI4_0421_contrast_g0\npI4_0421_contrast_g0_imec0',...
    'Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I3\neuropixels_data\npI3_0420_baseline_g0\npI3_0420_baseline_g0_imec0',...
    'Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I3\neuropixels_data\npI3_0421_contrast_g0\npI3_0421_contrast_g0_imec0'};
for dd = 1:length(directories)
    run_ks2(directories{dd});
end