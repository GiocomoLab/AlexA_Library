

addpath('Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\Analysis\MatlabImportExport_v6.0.0')
addpath('C:\code\spikes\preprocessing\phyHelpers');
addpath('C:\code\npy-matlab');
%eventsPath='F:\H1\2019-03-14_10-03-20\Events.nev';
eventsPath='F:\H1\2019-03-13_13-49-50\Events.nev';
cscPath ='Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\H1\CambridgeNeurotech\2019-03-13_13-49-50\CSC1.ncs';

%[ev_times, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV(eventsPath, [1 1 1 1 1], 1, 1, [] );
[ev_times, EventIDs, TTLs, Header] = Nlx2MatEV(eventsPath, [1 1 1 0 0], 1, 1, [] );
ev_times=ev_times-ev_times(1); %bc it starts at some weird non zero value
[sample_times,samples,header]=Nlx2MatCSC(cscPath, [1 0 0 0 1], 1, 1, [] );
sample_times=sample_times-sample_times(1);
%%
data_dir = 'F:\H1\2019-03-13_13-49-50';
session_name = '0313_contrast_2';

fid = fopen(fullfile(data_dir,strcat(session_name,'_position.txt')),'r');
vr_position_data = fscanf(fid, '%f', [3,inf])';
fclose(fid);
vr_ttl=vr_position_data(:,3);
vr_time=vr_position_data(:,2);

%%
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

%%
%frame_times_nl=ev_times(TTLs==1)/1000000;
frame_times_nl=ev_times(find(TTLs==1,1,'first'):find(TTLs==1,1,'last'))/1000000; %finds first up ttl and last up ttl and takes all timestamps in between (up and down)
frame_times_vr=vr_time;
%%
[mm,ii]=max(diff(frame_times_nl));
frametimes_nlOld = frame_times_nl;
if mm>1
    disp('Found large step in np sync data, assuming two sessions')
    %frame_times_nl=frame_times_nl(1:ii);
    frame_times_nl=frame_times_nl(ii+2:end);
end
figure
plot(diff(frame_times_nl),diff(vr_time(1:end-1)),'.')

%%
if abs(numel(frame_times_nl) - numel(frame_times_vr)) <= 1
    idx=1:min(numel(frame_times_nl),numel(frame_times_vr));
    post = frame_times_nl(idx)';
    posx = vr_position_data(idx,1);  
    
    % transform lick times into neuropixels reference frame
    beta = [ones(size(post)) frame_times_vr(idx)]\post;
    lickt = beta(1) + lickt*beta(2);
else
    disp('ERROR: number of sync pulses does not match number of frames.')
end
figure
plot(diff(post),diff(vr_time(idx)),'.')

%%
% load spike times
sp = loadKSdir(data_dir);

[aa,bb]=(max(diff(sample_times)));
start_offset=sample_times(bb)/1000000; %add this number to all sp

timegap=(sample_times(bb+1)-sample_times(bb))/1000000; %timegap between recordings, add this number to all spiketimes after offset
sp.st(sp.st>start_offset)=sp.st(sp.st>start_offset)+timegap;
%%
% shift everything to start at zero
offset = post(1);
post = post - offset;
lickt = lickt - offset;
sp.st = sp.st - offset;

%%
% resample position to have constant time bins
posx = interp1(post,posx,(0:0.02:max(post))');
post = (0:0.02:max(post))';
posx([false;diff(posx)<-2])=round(posx([false;diff(posx)<-2])/400)*400; % handle teleports
%%

% compute trial number for each time bin
trial = [1; cumsum(diff(posx)<-100)+1];

% throw out bins after the last trial
keep = trial<=num_trials;
trial = trial(keep);
posx = posx(keep);
post = post(keep);

% throw out licks before and after session
keep = lickt>min(post) & lickt<max(post);
lickx = lickx(keep);
lickt = lickt(keep);

% cut off all spikes before and after vr session
keep = sp.st >= 0 & sp.st <= post(end);
sp.st = sp.st(keep);
sp.spikeTemplates = sp.spikeTemplates(keep);
sp.clu = sp.clu(keep);
sp.tempScalingAmps = sp.tempScalingAmps(keep);

% save processed data
save(fullfile(data_dir,strcat(session_name,'.mat')),'sp','post','posx','lickt','lickx','trial','trial_contrast','trial_gain');
