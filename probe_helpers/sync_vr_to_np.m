% location of data
data_dir = 'F:\G4\1204_mismatch_1';
session_name = '1204_mismatch_1';

% get neuropixels sync pulse times
syncDat = extractSyncChannel(data_dir,277,277);
sync_sampling_rate = 2500; % from the LFP data file
frametimes_np = find(abs(diff(syncDat))==1)+1;
frametimes_np = frametimes_np/sync_sampling_rate;

% read vr position data
fid = fopen(fullfile(data_dir,strcat(session_name,'_position.txt')),'r');
vr_position_data = fscanf(fid, '%f', [3,inf])';
fclose(fid);

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

[mm,ii]=max(diff(frametimes_np));
frametimes_npOld = frametimes_np;
if mm>1
    disp('Found large step in np sync data, assuming two sessions')
    frametimes_np=frametimes_np(1:ii);
end

if abs(numel(frametimes_np) - numel(vr_position_data(:,2))) <= 1
    post = frametimes_np(1:numel(vr_position_data(:,2)))';
    posx = vr_position_data(:,1);  
    
    % transform lick times into neuropixels reference frame
    beta = [ones(size(post)) vr_position_data(:,2)]\post;
    lickt = beta(1) + lickt*beta(2);
else
    disp('ERROR: number of sync pulses does not match number of frames.')
end

% load spike times
sp = loadKSdir(data_dir);

% shift everything to start at zero
offset = post(1);
post = post - offset;
lickt = lickt - offset;
sp.st = sp.st - offset;

% resample position to have constant time bins
posx = interp1(post,posx,(0:0.02:max(post))');
post = (0:0.02:max(post))';
posx([false;diff(posx)<-2])=round(posx([false;diff(posx)<-2])/400)*400; % handle teleports

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