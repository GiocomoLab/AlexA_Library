%function sync_vr_to_np(data_dir)
addpath(genpath('C:\code\spikes'));
addpath(genpath('C:\code\npy-matlab'));

% location of data
data_dir = 'Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I4\neuropixels_data\npI4_0424_gaincontrast_g0';
adata_dir = 'Y:\giocomo\attialex\NP_DATA';

[~,main_name]=fileparts(data_dir);
animal_name = strsplit(main_name,'_');
animal_name = animal_name{1};
NIDAQ_file = fullfile(data_dir,strcat(main_name,'_t0.nidq.bin'));
NIDAQ_config = fullfile(data_dir,strcat(main_name,'_t0.nidq.meta'));
session_name = '0424_baseline_2';
spike_dir = fullfile(data_dir,strcat(main_name,'_imec0'));

%

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
syncDat=datNIDAQ(2,:)>1000;


frame_times_np = find(abs(diff(syncDat))==1)+1;
frame_times_np = frame_times_np/sync_sampling_rate;

% read vr position data
formatSpec = '%f%f%f%f%f%[^\n\r]';
delimiter = '\t';
fid = fopen(fullfile(data_dir,strcat(session_name,'_position.txt')),'r');
dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fid);
vr_position_data = cat(2,dataArray{1:5});
%vr_position_data = vr_position_data(1:49334,:);
nu_entries = nnz(~isnan(vr_position_data(1,:)));
%vr_position_data=vr_position_data(49335:end,:);
vr_ttl=vr_position_data(:,nu_entries); %assuming TTL in last and timestamp in second last column
frame_times_vr=vr_position_data(:,nu_entries-1);

% read vr trial data
fid = fopen(fullfile(data_dir,strcat(session_name,'_trial_times.txt')),'r');
vr_trial_data = fscanf(fid, '%f', [4,inf])';
fclose(fid);
trial_contrast = [100; vr_trial_data(:,2)];
trial_gain = [1; vr_trial_data(:,3)];
num_trials = numel(trial_gain);

% read vr licking data
fid = fopen(fullfile(data_dir,strcat(session_name,'_licks.txt')),'r');
vr_lick_data = fscanf(fid, '%f', [2,inf])';
fclose(fid);
lickx = vr_lick_data(:,1);
lickt = vr_lick_data(:,2);

% set vr frame times to be the time of neuropixels pulses
% make sure the number of frames matches (can be off by one because of
% odd/even numbers of frames)
%%
tmp_diff=diff(frame_times_np);
[mm,step_idx]=find(tmp_diff>2);
sess_length=diff([0 step_idx length(frame_times_np)]);
midpoint = ([0 step_idx] + [step_idx length(frame_times_np)])/2;
%step_idx=step_idx+1;
frametimes_nlOld = frame_times_np;
[~,ml]=min(abs(sess_length-numel(frame_times_vr)));
if length(mm)>=1
    figure;
    subplot(2,1,1)
    plot(frame_times_np)
    subplot(2,1,2)
    plot(tmp_diff)
    title(sprintf('found %d blocks',numel(step_idx+2)))
    hold on
    plot(step_idx,tmp_diff(step_idx),'ro')
    
    for im=1:numel(midpoint)
    text(midpoint(im),max(tmp_diff),sprintf('%d',sess_length(im)))
    end
    
    sess=input(sprintf('Which session do you want to extract (suggesting %d)',ml));
    is_mismatch = input('Is this a MM or PB sesion [0/1]?');
    step_idx = [0 step_idx length(frame_times_np)];
    %frame_times_np=frame_times_np(ii+2:end);
    idx_start=step_idx(sess)+1;
    idx_stop = step_idx(sess+1);
    frame_times_np=frame_times_np(idx_start:idx_stop);
else
    is_mismatch=0;
end
%%
if abs(numel(frame_times_np) - numel(frame_times_vr)) <= 1
    idx=1:min(numel(frame_times_np),numel(frame_times_vr)); %use shorter index
    post = frame_times_np(idx)';
    vr_position_data=vr_position_data(idx,:);
    posx = vr_position_data(:,1);
    
    % transform lick times into neuropixels reference frame
    beta = [ones(size(post)) frame_times_vr(idx)]\post;
    lickt = beta(1) + lickt*beta(2);
else
    disp('ERROR: number of sync pulses does not match number of frames.')
end
figure
scatter(diff(post),diff(frame_times_vr(idx)),2,1:length(idx)-1)
%%
% load spike times
%sp = loadKSdir(spike_dir);
load(fullfile(adata_dir,strcat(animal_name,'_',session_name,'.mat')),'sp');

offset = post(1);

sp.vr_session_offset = offset;

save(fullfile(adata_dir,strcat(animal_name,'_',session_name,'.mat')),'sp','-append');
