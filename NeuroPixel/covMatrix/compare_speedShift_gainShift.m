ops.chunksize=100; %in bins,so thats 200 cm
ops.stride_start = 1;
ops.binsize=2;
ops.stride = 99;
startVec = ops.stride_start:ops.stride:(400/ops.binsize-ops.chunksize+1);
chunksPerTrials = numel(startVec);
ops.region = '';
ops.gain_to_look_at = .8;
ops.contrast = 100;
ops.tt=(-6:9);
ops.template_trials = 1:6;
ops.stability_threshold=0.2;
ops.smoothSigma = 4;
%%
[filenames,triggers] = getFilesCriteria(ops.region,ops.contrast,ops.gain_to_look_at,'/oak/stanford/groups/giocomo/attialex/NP_DATA');
%%
%shift_path = 'Z:\giocomo\attialex\speed_shiftFilter';
shift_path = '/oak/stanford/groups/giocomo/attialex/speed_shiftFilter';
params = load(fullfile(shift_path,'parameters.mat'));
ALLSHIFTS = [];
REGIONS = {};
ALLFACTORS = [];
STABILITY = [];
ESTIMATED = [];
ID=[];
for iF=1:numel(filenames)
    [~,sn]=fileparts(filenames{iF});
    if ~isfile(fullfile(shift_path,[sn '.mat']))
        continue
    end
   


data = load(filenames{iF});
shift_data = load(fullfile(shift_path,[sn '.mat']));

trials = triggers{iF}(1)+ops.tt;
[PEAKS,SHIFTS,speed_mat,stability]=calculatePeakShift(data,trials,ops);
dat = shift_data.stability;
[ma,mi]=max(dat,[],2);
best_delay = params.ops.factors(mi); %in seconds
shifts = nanmean(SHIFTS,3);
shifts = shifts(:,7:10)*ops.binsize; % shifts for each cell during gain trials in cm
bl_speed = mean(mean(speed_mat(1:6,:)));
gain_speed = mean(speed_mat(7:10,:),2);
delta_speed=gain_speed-bl_speed;
ALLSHIFTS = cat(1,ALLSHIFTS,SHIFTS(:,7:10,:)*ops.binsize);
REGIONS = cat(1,REGIONS,shift_data.region');
ALLFACTORS = cat(1,ALLFACTORS,best_delay');
est=(gain_speed-bl_speed)*best_delay;
ESTIMATED = cat(1,ESTIMATED,est');
STABILITY=cat(1,STABILITY,stability');
ID = cat(1,ID,iF*ones(numel(stability),1));
end
%%
figure
idx = STABILITY>0 & ismember(REGIONS,'VISp');
for iT=1:4
    subplot(1,4,iT)
    
    scatter(nanmean(ALLSHIFTS(idx,iT,:),3),ESTIMATED(idx,iT),'.')
    axis image
    xlim([-10 10])
    ylim([-10 10])
end
