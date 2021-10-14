mec_clu = data.anatomy.cluster_id(startsWith(data.anatomy.cluster_parent,'RHP'));
good_cells = data.sp.cids(data.sp.cgs==2);
good_cells = good_cells(ismember(good_cells,mec_clu));
%good_cells = data.sp.ks_cluster.cluster_id(startsWith(data.sp.ks_cluster.KSLabel,'good'))
pos_idx = discretize(data.sp.st,data.post);
data.trial = cumsum([0;diff(data.posx)<-100]);
%stim = round(data.vr_data_resampled.Object);
stim = ones(size(data.posx));

opt = load_mismatch_opt;
opt.TrackEnd=round(max(data.posx/10))*10;
opt.track_length=opt.TrackEnd;
opt.max_lag = 30;
cm = trialCorrMat(good_cells,0:max(data.trial)-1,data,opt);

frMat = calcTrialFRMat(good_cells,0:max(data.trial)-1,data,opt);
%%
if isfield(data.vr_data_resampled,'Object')
    object_field_name = 'Object';
else
    object_field_name = 'triggerObject';
end
stim_type = ones(1,max(data.trial));
for iT=1:max(data.trial)
    idx=data.trial==(iT-1);
    stim_type(iT)=max(data.vr_data_resampled.(object_field_name)(idx));
end
%%
stab = nanmean(nanmean(cm(:,2:10,2:10),2),3);
[~,sid] = sort(stab,'descend','MissingPlacement','last');
cmap = brewermap(3,'Set2');
cmap = parula(3);
figure

stim = round(data.vr_data_resampled.(object_field_name));
stim(data.vr_data_resampled.visible>0.5)=-1;
cmap = parula(numel(unique(stim)));
for iC=1:numel(good_cells)
    idx = data.sp.clu==good_cells(sid(iC));
    this_idx = pos_idx(idx);
    subplot(2,1,1)
    scatter(data.posx(this_idx),data.trial(this_idx),12,cmap(stim(this_idx)+2,:),'.')
    xlim([0,opt.TrackEnd])
    ylim([0 max(data.trial)])
    subplot(2,1,2)
    mean_bl=squeeze(nanmean(frMat(sid(iC),stim_type==1,:),2));
    mean_other = squeeze(nanmean(frMat(sid(iC),stim_type==0,:),2));
    x_vec = 0:2:opt.TrackEnd;
    x_vec = x_vec(1:size(mean_bl,1));
    plot(x_vec,mean_bl)
    hold on
    plot(x_vec,mean_other);
    pause
    cla
    
end

%%
figure
for ii=0:1
    trigs = strfind(data.vr_data_resampled.(object_field_name)==ii,[0 1]);
    [spikeTimes,~,aux,~,count_vec]=extract_triggered_spikeTimes(data.sp,data.post(trigs),'cluIDs',good_cells,'win',opt.extract_win);
    hold on
    plot(opt.time_vecs,nanmean(count_vec))
end
%%
figure
[~,sid]=sort(stim_type);
imagesc(squeeze(nanmean(cm(stab>0.5,sid,sid))),[0 0.7])
