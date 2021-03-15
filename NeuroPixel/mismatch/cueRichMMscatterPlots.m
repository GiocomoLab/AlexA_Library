
good_cells = data.sp.cids(data.sp.cgs==2);
%good_cells = data.sp.ks_cluster.cluster_id(startsWith(data.sp.ks_cluster.KSLabel,'good'))
pos_idx = discretize(data.sp.st,data.post);
data.trial = cumsum([0;diff(data.posx)<-100]);
%stim = round(data.vr_data_resampled.Object);




opt = load_mismatch_opt;
opt.TrackEnd=round(max(data.posx/10))*10;
opt.track_length=opt.TrackEnd;
opt.max_lag = 30;
cm = trialCorrMat(good_cells,0:19,data,opt);

frMat = calcTrialFRMat(good_cells,0:max(data.trial)-1,data,opt);

%%
stab = nanmean(nanmean(cm,2),3);
[~,sid] = sort(stab,'descend','MissingPlacement','last');
cmap = brewermap(3,'Set2');
cmap = parula(3);
figure
% [~,~,stim]=unique(round(data.mismatch_trigger*10));
data.mismatch_trigger = data.vr_data_resampled.MM;
stim=zeros(size(data.mismatch_trigger));
stim(data.mismatch_trigger==0)=1;
stim(data.mismatch_trigger==0.5)=2;
% stim(data.mismatch_trigger == 0.8)=3;
cmap = parula(numel(unique(stim)));
cmap = flipud(cmap);
%cmap = brewermap(numel(unique(stim)),'Set3')
for iC=1:numel(good_cells)
    idx = data.sp.clu==good_cells(sid(iC));
    this_idx = pos_idx(idx);
    scatter(data.posx(this_idx),data.trial(this_idx),12,cmap(stim(this_idx)+1,:),'.')
   
    pause
    cla
end

