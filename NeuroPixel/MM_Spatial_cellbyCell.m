% depthBinSize = 80; % in units of the channel coordinates, in this case µm
% timeBinSize = 0.01; % seconds
% bslWin = [-0.2 -0.05]; % window in which to compute "baseline" rates for normalization
% psthType = 'norm'; % show the normalized version
% eventName = 'stimulus onset'; % for figure labeling
% 
% [timeBins, depthBins, allP, normVals] = psthByDepth(sp.st, spikeDepths, ...
%     depthBinSize, timeBinSize, eventTimes, window, bslWin);
% 
% figure;
% plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);
%%

addpath(genpath('C:\Users\giocomolab\Dropbox\Work\neuropixels\functions'));
addpath(genpath('C:\Users\giocomolab\Dropbox\Work\neuropixels\spikes'));

% where to find data.
% these are the sessions with histology
data_dir = 'C:\Users\giocomolab\Dropbox\Work\neuropixels\data\';
%%
session_num =1;
session_name = {'npG2_1211_gain_1',...
    'npG2_1212_gaincontrast_1',...
    'npG5_1207_gain_1',...
    'npG5_1210_gaincontrast_2'};

% some params
params = readtable('C:\Users\giocomolab\Dropbox\Work\neuropixels\UniversalParams.xlsx');
filenames = {'G2/1211_mismatch_1/1211_mismatch_1.mat',...
    'G2/1212_mismatch_1/1212_mismatch_1.mat',...
    'G5/1207_mismatch_1/1207_mismatch_1.mat',...
    'G5/1210_mismatch_1/1210_mismatch_1.mat'
    };

%%
    load(fullfile(data_dir,strcat(session_name{session_num},'.mat')));

trials=[1:max(trial)];
spC=[];
dwell_time=[];
for iT=1:length(trials)
    idxVR=trial==trials(iT);
    t_time=post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=sp.st<stop & sp.st>=start;
    edges=[0:2:402];
    edges(1)=-5;
    [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',edges);
    spC=cat(3,spC,spM);
    dwell_time=cat(1,dwell_time,dT);
end
%cellIDX=find(sp.cgs>=1);
spC=spC(sp.cids+1,:,:);
dt=dwell_time';

dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spC,1)
    spC(ii,:,:)=spC(ii,:,:)./dt;
end
spC=spC/0.02;

load(fullfile('F:',filenames{session_num}))

speed=true_speed;
speed_t=0.05;
% figure('Name',filenames{iF});; plot(speed)
% 
all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
run_periods=smooth(speed,25)>speed_t;
run_window=-30:30;
possibles=strfind(run_periods',ones(1,length(run_window)))+floor(.5*length(run_window));


mm_trigs=all_mm_trigs(ismember(all_mm_trigs,possibles));

window = [-1 4]; % look at spike times from 0.3 sec before each event to 1 sec after

% if your events come in different types, like different orientations of a
% visual stimulus, then you can provide those values as "trial groups",
% which will be used to construct a tuning curve. Here we just give a
% vector of all ones. 
trialGroups = [ones(size(mm_trigs))];
eventTimes=[post(mm_trigs)'];
%%

psthViewer_Alex(sp, eventTimes, window, trialGroups,spC);

%cid=sp.cids(idx);plot(sp.st(sp.clu==cid),sp.tempScalingAmps(sp.clu==cid),'.')