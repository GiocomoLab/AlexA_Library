%waveform extractor
addpath(genpath('/home/users/attialex/neuropixel-utils'))
addpath(genpath('/home/users/attialex/AlexA_Library'));
%% get files
%dir('/oak/stanford/groups/giocomo/export/data/Projects/AlexA_NP/**/AA_190830_1_0927_gaincontrast_1_g0_t0.imec0.ap.bin')

%ff=dir('/oak/stanford/groups/giocomo/export/data/Projects/AlexA_NP/**/*ap.bin')
%%
%gg = dir('/oak/stanford/groups/giocomo/export/data/Projects/ContrastExperiment_neuropixels/**/*ap.bin');
%gg = dir('/oak/stanford/groups/giocomo/export/data/Projects/ContrastExperiment_neuropixels/**/*_CAR.bin');
%%
%file_location = cat(1,ff,gg);
%save('/oak/stanford/groups/giocomo/attialex/FileLocations.mat','file_location');
%%
load('/oak/stanford/groups/giocomo/attialex/FileLocations')

% %%
% imec_dir = getRawDataPath('npI5_0413_baseline_g0_t0.imec0.ap.bin',file_location);
%
% ks_dir = fileparts(imec_dir)
% pp=get_waveforms(ks_dir,imec_dir,385)
%%
filenames = dir('/oak/stanford/groups/giocomo/attialex/NP_DATA/*.mat');
%%
savedir = '/oak/stanford/groups/giocomo/attialex/mean_waveforms';
if ~isfolder(savedir)
    mkdir(savedir)
end
%%
p=gcp('nocreate')
if isempty(p)
    p=parpool(3);
end

%%
parfor iF =1:numel(filenames)
    try
        [~,session_name]=fileparts(filenames(iF).name);
        if isfile(fullfile(savedir,[session_name '.mat']))
            a=who('-file',fullfile(savedir,[session_name '.mat']));
            if ismember('mean_waveforms',a)
                sprintf('%s already exists \n',session_name)
                continue
            end
        end
        data = load(fullfile(filenames(iF).folder,filenames(iF).name));
        dat_path = data.sp.dat_path;
        n_chan = data.sp.n_channels_dat;
        %find full path
        imec_dir = getRawDataPath(dat_path,file_location);
        ks_dir = fileparts(imec_dir);
        
        [mean_waveforms,mean_template_waveforms,amplitude,aux,clusters_good]=get_waveforms(ks_dir,imec_dir,n_chan);
        
        %data.mean_waveforms = mean_waveforms;
        %data.mean_template_waveforms = mean_template_waveforms;
        m = matfile(fullfile(savedir,session_name),'Writable',true);
        m.mean_waveforms = mean_waveforms;
        m.mean_template_waveforms = mean_template_waveforms;
        m.amplitudes = amplitude;
        m.aux = aux;
        m.clusters_good = clusters_good;
        % save(fullfile(savedir,session_name),'mean_waveforms','mean_waveforms','amplitude','aux')
    catch ME
        disp(ME.message);
        a=fullfile(savedir,sprintf('%s_error.txt',session_name));
        disp(a);
        fid = fopen(a,'w');
        fprintf(fid,ME.message)
        fclose(fid);
    end
    
end
