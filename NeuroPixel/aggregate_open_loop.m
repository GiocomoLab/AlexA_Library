files = dir('Z:\giocomo\attialex\NP_DATA\CL_OL_Towers\*_mismatch*.mat');


mismatch=load(fullfile(files(iF).folder,files(iF).name));
playback=load(fullfile(files(iF).folder,strrep(files(iF).name,'mismatch','playback')));
good_cells = mismatch.sp.cids(mismatch.sp.cgs==2);

corrs_baseline = zeros(1,max(mismatch.trial));

for ii = 1:max(mismatch.trial)
    idx = mismatch.trial==ii & mismatch.posx>2 & mismatch.posx<398;
    corrs_baseline(ii)=corr(mismatch.true_speed(idx),[0;diff([mismatch.posx(idx)])]);
end

corrs_playback = zeros(1,max(playback.trial));

for ii = 1:max(playback.trial)
    idx = playback.trial==ii & playback.posx>2 & playback.posx<398;
    corrs_playback(ii)=corr(playback.true_speed(idx),[0;diff([playback.posx(idx)])]);
end

difference_by_location_playback = zeros(max(playback.trial),100);
difference_by_location_mismatch = zeros(max(mismatch.trial),100);


tmp_diff = playback.true_speed-[0;diff(playback.posx)];
for iT=1:max(playback.trial)
    for iP = 1:numel(pos_bins)
        idx = playback.trial == iT & discrete_pos == iP;
        difference_by_location_playback(iT,iP)=mean(abs(tmp_diff(idx)));
    end
end

tmp_diff = mismatch.true_speed-[0;diff(mismatch.posx)];
  discrete_pos = discretize(mismatch.posx,pos_bins);

for iT=1:max(mismatch.trial)
    for iP = 1:numel(pos_bins)
        idx = mismatch.trial == iT & discrete_pos == iP;
        difference_by_location_mismatch(iT,iP)=mean(abs(tmp_diff(idx)));
    end
end


[a,b]=sort(corrs_playback,'descend');
ranking=1:numel(b);
ranking(b) = ranking;


[~,spatialMap_mm] = get_spatial_map(mismatch);
[~,spatialMap_pb] = get_spatial_map(playback);
%% stability in baseline
tc_1=mean(spatialMap_mm(:,2:98,1:40),3);
tc_2=mean(spatialMap_mm(:,2:98,41:end),3);
tmp = corr(tc_1',tc_2');
stability=diag(tmp);
%[~,stability_sort]=sort(stability);



%%
[a,b]=sort(corrs_playback,'descend');
t_c_mm=mean(spatialMap_mm,3);
similar_idx = b(1:10);
diff_idx = b(end-9:end);
t_c_similar = mean(spatialMap_pb(:,:,similar_idx),3);
t_c_different = mean(spatialMap_pb(:,:,diff_idx),3);

%% stability in baseline
tc_1=mean(spatialMap_mm(:,2:98,1:40),3);
tc_2=mean(spatialMap_mm(:,2:98,41:end),3);
tmp = corr(tc_1',tc_2');
stability=diag(tmp);
[~,stability_sort]=sort(stability);
%%

score_similar = zeros(1,numel(good_cells));
score_different = zeros(1,numel(good_cells));
for ii = 1:numel(good_cells)
    score_similar(ii)=corr(t_c_mm(ii,:)',t_c_similar(ii,:)');
    score_different(ii)=corr(t_c_mm(ii,:)',t_c_different(ii,:)');
end