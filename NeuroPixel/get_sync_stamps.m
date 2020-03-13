
%%
function [ts_LFP,ts_NIDAQ]=get_sync_stamps(data_dir)
[~,main_name]=fileparts(data_dir);
animal_name = strsplit(main_name,'_');
animal_name = animal_name{1};

NIDAQ_file = dir(fullfile(data_dir,'*nidq.bin'));
NIDAQ_file = fullfile(data_dir,NIDAQ_file(1).name);
NIDAQ_config = dir(fullfile(data_dir,'*nidq.meta'));
NIDAQ_config = fullfile(data_dir,NIDAQ_config(1).name);

%NIDAQ_file = fullfile(data_dir,strcat(main_name,'_t0.nidq.bin'));
%NIDAQ_config = fullfile(data_dir,strcat(main_name,'_t0.nidq.meta'));

spike_dir = fullfile(data_dir,strcat(main_name,'_imec0'));

%% timestamp nidaq

%get the nidaq sample rate & get number of recorded nidaq channels
dat=textscan(fopen(NIDAQ_config),'%s %s','Delimiter','=');
names=dat{1};
vals=dat{2};
loc=contains(names,'niSampRate');
sync_sampling_rate=str2double(vals{loc});

loc2=contains(names,'nSavedChans');
n_channels_nidaq=str2double(vals{loc2});

% get neuropixels sync pulse times
fpNIDAQ=fopen(NIDAQ_file);
datNIDAQ=fread(fpNIDAQ,[n_channels_nidaq,Inf],'*int16');
fclose(fpNIDAQ);
syncDatNIDAQ=datNIDAQ(1,:)>1000;
%% timestamp lfp
LFP_config = dir(fullfile(spike_dir,'*.lf.meta'));
LFP_config = fullfile(LFP_config.folder,LFP_config.name);

LFP_file = dir(fullfile(spike_dir,'*.lf.bin'));
LFP_file = fullfile(LFP_file.folder,LFP_file.name);

dat=textscan(fopen(LFP_config),'%s %s','Delimiter','=');
names=dat{1};
vals=dat{2};
loc=contains(names,'imSampRate');
lfp_sampling_rate=str2double(vals{loc});

% get neuropixels sync pulse times
fpLFP=fopen(LFP_file);
fseek(fpLFP,384*2,0);
ftell(fpLFP);
datLFP=fread(fpLFP,[1,inf],'*int16',384*2);
fclose(fpLFP);
%figure
%plot(datLFP)
syncDatLFP=datLFP(1,:)>10;
%%
ts_NIDAQ = strfind(syncDatNIDAQ,[0 1])/sync_sampling_rate;
ts_LFP = strfind(syncDatLFP,[0 1])/lfp_sampling_rate;

end
