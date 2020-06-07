%%
%calculate firing rate
fr = calcFRVsTime(good_cells,data,ops,ops_here.trials);

% threshold by running speed and trials that we want to look at
speed = calcSpeed(data.posx,ops);
trial_idx = ismember(data.trial,ops_here.trials);
idx = speed>ops.SpeedCutoff & trial_idx;
X = fr(:,idx);
trial = data.trial(idx);
posx = data.posx(idx);
speed = speed(idx);

% normalize each cell to have FR between 0 and 1
X = (X-nanmin(X,[],2))./repmat(nanmax(X,[],2)-nanmin(X,[],2),1,size(X,2));

% set nans to 0
X(isnan(X)) = 0;

% mean subtract
Xtilde = X - nanmean(X,2);

%prepare y
num_components = size(Xtilde,1);
[~,~,posbin] = histcounts(posx,ops.xbinedges);
posbin(posbin==0) = 1;

%subsample for other classifiers
sub_sample_idx = false(size(trial));
sub_sample_idx(1:1:end)=true; %1s p s


train_trials = ops_here.trials(ops.trials_train);
train_trial_idx=ismember(trial,train_trials);


tc = nan(num_components,ops.nBins,1);
for i = 1:ops.nBins
    tc(:,i) = mean(Xtilde(:,posbin==i & train_trial_idx),2);
end

% decode position in gain change trials based on baseline1 trials
dot_prod = tc' * Xtilde; % predict position
[~,max_bin] = max(dot_prod);
pred_pos = ops.xbincent(max_bin);
train_trial_idx = train_trial_idx & sub_sample_idx;
%Mdl = fitlm(Xtilde(:,train_trial_idx)',posbin(train_trial_idx));
Mdl = fitcecoc(Xtilde(:,train_trial_idx)',posbin(train_trial_idx));
%yhat{iFold} = ops.xbincent(predict(Mdl,Xtilde(:,test_trial_idx)'));
yhat = ops.xbincent(mod(round(predict(Mdl,Xtilde')),ops.track_length/2)+1);


decoder.cluID = good_cells;
decoder.region = decode_region;
decoder.pred_pos = pred_pos;
decoder.true_pos = posx;
decoder.trial = trial;
decoder.speed = speed;
decoder.yhat=yhat;
