channelMapFile = 'C:\code\KiloSort2\configFiles\neuropixPhase3B1_kilosortChanMap.mat';
%ksdirs = dir('E:\**\*ap.bin');
%ksdirs=dir('Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\**\*imec0.ap.bin');
for iD = 1:length(ksdirs)
    clearvars -except ksdirs channelMapFile iD
    im_save_dir = 'Y:\giocomo\attialex\images\LFP';
    myKsDir=ksdirs(iD).folder;
    sprintf('now working on %s',myKsDir)
    extract_plotLFP
    
end