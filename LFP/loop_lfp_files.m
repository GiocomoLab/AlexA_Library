addpath(genpath('/home/users/attialex/neuropixel-utils'))
channelMapFile = fullfile('/home/users/attialex/neuropixel-utils/map_files','neuropixPhase3B1_kilosortChanMap.mat');
%ksdirs = dir(fullfile('/oak/stanford/groups/giocomo','export','data','Projects','RandomForage_NPandH3','ProcessedData','/Hanover*/*ap.bin'));
%ksdirs = dir(fullfile('/oak/stanford/groups/giocomo/export/data/Projects/AlexA_NP/**','*ap.bin'))
%ksdirs=dir('Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\**\*imec0.ap.bin');
ksdirs = dir(fullfile('/oak/stanford/groups/giocomo/export/data/Projects/AlexA_NP/AA_*/neuropixels_data/AA*/AA*/*.ap.bin'));
for iD = 1%1:length(ksdirs)
    clearvars -except ksdirs channelMapFile iD 
    im_save_dir = fullfile('/oak/stanford/groups/giocomo','attialex','images','LFP');
    myKsDir=ksdirs(iD).folder;
    sprintf('now working on %s',myKsDir)
    try
    extract_plotLFP
    catch ME
        display(ME.message)
        sprintf('Failed for %d, %s',iD,ksdirs(iD).name)
    end
    
end