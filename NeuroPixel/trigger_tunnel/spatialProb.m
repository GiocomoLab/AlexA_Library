
good_cells = data.sp.cids(data.sp.cgs==2);
%good_cells = data.sp.ks_cluster.cluster_id(startsWith(data.sp.ks_cluster.KSLabel,'good'))
pos_idx = discretize(data.sp.st,data.post);
data.trial = cumsum([0;diff(data.posx)<-100]);
%stim = round(data.vr_data_resampled.Object);
stim = ones(size(data.posx));
%%
figure
cm = brewermap(3,'Set2');
for iC=1:numel(good_cells)
    idx = data.sp.clu==good_cells(iC);
    this_idx = pos_idx(idx);
    
    scatter(data.posx(this_idx),data.trial(this_idx),12,cm(stim(this_idx)+2,:))
    pause
    cla
end
   

%%
figure
scatter(data.posx,data.trial,12,data.vr_data_resampled.Stim)
%%
opt = load_mismatch_opt;
opt.TrackEnd=round(max(data.posx/10))*10;
opt.track_length=opt.TrackEnd;
opt.max_lag = 30;
cm = trialCorrMat(good_cells,0:19,data,opt);

frMat = calcTrialFRMat(good_cells,0:max(data.trial)-1,data,opt);
%%
stim_type = ones(1,max(data.trial));
for iT=1:max(data.trial)
    idx=data.trial==(iT-1);
    stim_type(iT)=max(data.vr_data_resampled.Stim(idx));
end
%%
stab = nanmean(nanmean(cm,2),3);
[~,sid] = sort(stab,'descend','MissingPlacement','last');
cmap = brewermap(3,'Set2');
cmap = parula(3);
figure
% [~,~,stim]=unique(round(data.mismatch_trigger*10));
% stim=zeros(size(data.mismatch_trigger));
% stim(data.mismatch_trigger==0)=1;
% stim(data.mismatch_trigger==0.5)=2;
% stim(data.mismatch_trigger == 0.8)=3;
stim = round(data.vr_data_resampled.Stim);
cmap = parula(numel(unique(stim)));
for iC=1:numel(good_cells)
    idx = data.sp.clu==good_cells(sid(iC));
    this_idx = pos_idx(idx);
    subplot(2,1,1)
    scatter(data.posx(this_idx),data.trial(this_idx),12,cmap(stim(this_idx)+1,:),'.')
    subplot(2,1,2)
    mean_bl=squeeze(nanmean(frMat(sid(iC),stim_type==2,:),2));
    mean_other = squeeze(nanmean(frMat(sid(iC),stim_type==1,:),2));
    plot(mean_bl)
    hold on
    plot(mean_other);
    pause
    cla
end

%%
figure
for ii=1:2
    trigs = strfind(data.vr_data_resampled.Stim==ii,[0 1]);
    [spikeTimes,~,aux,~,count_vec]=extract_triggered_spikeTimes(data.sp,data.post(trigs),'cluIDs',good_cells,'win',opt.extract_win);
    hold on
    plot(opt.time_vecs,nanmean(count_vec))
end