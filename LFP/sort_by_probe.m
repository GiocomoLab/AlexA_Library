%metadirs=dir('Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\**\*imec0.ap.meta');
for iD = 1:length(metadirs)
    im_save_dir = 'Y:\giocomo\attialex\images\LFP';

    try
    myMetaDir=metadirs(iD).folder;
    AP_config = fullfile(myMetaDir,metadirs(iD).name);
    %sprintf('now working on %s',myMetaDir)
     %fclose(AP_config)
    f=fopen(AP_config);
    dat=textscan(f,'%s %s','Delimiter','=');
    fclose(f);
names=dat{1};
vals=dat{2};
loc=contains(names,'imDatPrb_sn');
probe_sn=vals{loc};
    
    probe_folder = fullfile('Y:\giocomo\attialex\images\LFP',probe_sn);

if ~isfolder(probe_folder)
    mkdir(probe_folder)
end
session = split(myMetaDir,filesep);
session = session{end-1};
session = session(1:end-3);
%get image file
im_file = fullfile(im_save_dir,session);
type = 'LFP';
im_file = strcat(im_file,strcat(type,'.png'));
copyfile(im_file,fullfile(probe_folder,strcat(session,type,'.png')));

    catch ME
        sprintf('Failed for %d, %s',iD,metadirs(iD).name)
    end
    
end
%%


