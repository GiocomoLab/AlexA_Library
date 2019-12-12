params = readtable('UniversalParams.xlsx');

good_cells = data.sp.cids(data.sp.cgs==2);
[~,speed_raw]=calcSpeed(data.posx,params);
posxhat = data.posx+-0.2*speed_raw;
posxhat = mod(posxhat,400);
figure
for cellIDX =b'%1:numel(good_cells)
subplot(2,1,1)

spike_id=data.sp.clu==good_cells(cellIDX);
spike_t = data.sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,data.post);
scatter(data.posx(spike_idx),data.trial(spike_idx),2)
xlim([0 400])
subplot(2,1,2)
scatter(posxhat(spike_idx),data.trial(spike_idx),2)
xlim([0 400])
pause
clf
end

%%
trialMap = nan(1,212);
cntr = 1;
bl_trials = find(data.trial_gain==1 & data.trial_contrast==100);
for iT =1:212
    if ismember(iT,bl_trials)
        trialMap(iT)=cntr;
        cntr=cntr+1;
    end
end
%%
trial_sorted = nan(size(data.trial));
for iT=1:numel(trial_sorted)
    trial_sorted(iT)=trialMap(data.trial(iT));
end
%%
figure

data.trial_sorted = trial_sorted;
subplot(2,1,1)
data.posWindow = [0 400];
plotRasterSpeedSort(data,params,[],good_cells(38),find(data.trial_gain==1),subplot(1,1,1))
posx=data.posx;
[~,speed_raw]=calcSpeed(data.posx,params);



%%
factors = [-.5: 0.05:0.5];

figure;
for iFactor = 1:numel(factors)%1:numel(good_cells)
ax=subplot(1,1,1);
posxhat = posx+factors(iFactor)*speed_raw;
posxhat = mod(posxhat,400);
data.posx=posxhat;
plotRasterSpeedSort(data,params,[],good_cells(38),find(data.trial_gain==1),ax)

% 
% spike_id=data.sp.clu==good_cells(cellIDX);
% spike_t = data.sp.st(spike_id);
% [~,~,spike_idx] = histcounts(spike_t,data.post);
% scatter(data.posx(spike_idx),data.trial(spike_idx),2)
% xlim([0 400])
% subplot(2,1,2)
% scatter(posxhat(spike_idx),data.trial(spike_idx),2)
% xlim([0 400])
% pause

pause
clf
end

%%
good_cells=data.sp.cids(data.sp.cgs==2);
params = readtable('UniversalParams.xlsx');
figure
edges = 0:1:400;

trialMap = nan(1,212);
cntr = 1;
bl_trials = find(data.trial_gain==1 & data.trial_contrast==100);
for iT =1:212
    if ismember(iT,bl_trials)
        trialMap(iT)=cntr;
        cntr=cntr+1;
    end
end

trial_sorted = nan(size(data.trial));
for iT=1:numel(trial_sorted)
    trial_sorted(iT)=trialMap(data.trial(iT));
end
[~,speed_raw]=calcSpeed(data.posx,params);
for cellIDX=30:numel(good_cells);

spike_id=data.sp.clu==good_cells(cellIDX);
spike_t = data.sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,data.post);
posx=mod(data.posx,400);
spike_loc = discretize(posx,edges);

idx=triu(true(nnz(data.trial_gain==1)),1);

spMat = zeros(nnz(data.trial_gain==1),400);
for ii=1:numel(spike_idx)
    if ~isnan(trial_sorted(spike_idx(ii)))
    spMat(trial_sorted(spike_idx(ii)),spike_loc (spike_idx(ii)))=spMat(trial_sorted(spike_idx(ii)),spike_loc(spike_idx(ii)))+1;
    end
end
variability = var(spMat);
variability=mean(variability);
cc=corr(spMat');
stability=mean(cc(idx));

factors = [-.25: 0.01:0.25];
VAR=zeros(size(factors));
STAB=VAR;
for iFactor = 1:numel(factors)
spMatHat = zeros(size(spMat));
posxhat = posx+factors(iFactor)*speed_raw;
posxhat = mod(posxhat,400);
spike_loc_hat = discretize(posxhat,edges);
for ii=1:numel(spike_idx)
    if ~isnan(trial_sorted(spike_idx(ii)))
    spMatHat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))=spMat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))+1;
    end
end

variability_hat = var(spMatHat);
variability_hat=mean(variability_hat);
cc=corr(spMatHat');
stability=mean(cc(idx));

VAR(iFactor)=variability_hat;
STAB(iFactor)=stability;
end
subplot(4,1,1)
plot(factors,STAB)

[ma,mi]=max(STAB);
posxhat = posx+factors(mi)*speed_raw;
posxhat = mod(posxhat,400);
posxhat2 = posx+-0.1*speed_raw;
posxhat2=mod(posxhat2,400);

subplot(4,1,2)
scatter(posx(spike_idx),trial_sorted(spike_idx),2);
xlim([0 400])
ylim([0 176])
subplot(4,1,3)
scatter(posxhat(spike_idx),trial_sorted(spike_idx),2);
xlim([0 400])
ylim([0 176])

subplot(4,1,4)
scatter(posxhat2(spike_idx),trial_sorted(spike_idx),2);
xlim([0 400])
ylim([0 176])
pause
clf
end
%%
