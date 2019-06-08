%mismatch=load('Z:\giocomo\attialex\NP_DATA\mismatch\npI1_0418_mismatch_4.mat');
%playback=load('Z:\giocomo\attialex\NP_DATA\mismatch\npI1_0418_playback_4.mat');
%mismatch=load('Z:\giocomo\attialex\NP_DATA\CL_OL_Towers\npAA1R_0419_mismatchTowers_1.mat');
%playback=load('Z:\giocomo\attialex\NP_DATA\CL_OL_Towers\npAA1R_0419_PlaybackTowers_1.mat');
%gain=load('Z:\giocomo\attialex\NP_DATA\mismatch\npI1_0418_gain_3.mat');
%%
[~,maps_mm] = compute_spatial_similarity(mismatch);
%% combine rasters
good_cells = mismatch.sp.cids(mismatch.sp.cgs==2);
figure;
for iC=1:length(good_cells)
spike_id = mismatch.sp.clu==good_cells(iC);
    spike_t = mismatch.sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,mismatch.post);
     
    scatter(mismatch.posx(spike_idx),mismatch.trial(spike_idx),2,'b')
    
    spike_id = playback.sp.clu==good_cells(iC);
    spike_t = playback.sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,playback.post);
     hold on
    scatter(playback.posx(spike_idx),playback.trial(spike_idx)+max(mismatch.trial),2,'r')
    
    pause
    cla;
end
% %%
% idx = mismatch.trial==5;
% figure
% plot(mismatch.true_speed(idx))
% hold on
% plot([0;diff(mismatch.posx(idx))])
%%
delays=zeros(1,max(mismatch.trial));
corrs_baseline = zeros(1,max(mismatch.trial));
run_speed_mismatch=[0;diff(mismatch.posx)];
invalid=find(run_speed_mismatch<-50);
run_speed_mismatch(invalid) = 0.5*run_speed_mismatch(invalid-1)+0.5*run_speed_mismatch(invalid+1);
%run_speed(invalid)=interp1(1:length(run_speed),run_speed,invalid);
run_speed_playback=[0;diff(playback.posx)];
invalid=find(run_speed_playback<-5 | run_speed_playback>5);
run_speed_playback(invalid) = 0.5*run_speed_playback(invalid-1)+0.5*run_speed_playback(invalid+1);
for ii = 1:max(mismatch.trial)
    idx = mismatch.trial==ii;
    delays(ii)=finddelay(mismatch.true_speed(idx),[0;diff([mismatch.posx(idx)])],5);
    corrs_baseline(ii)=corr(mismatch.true_speed(idx),[0;diff([mismatch.posx(idx)])]);
    
end
figure
subplot(2,1,1)    
scatter(mismatch.true_speed,run_speed_mismatch)
    xlim([-.5 3])
    ylim([-.5 3])
    title('correlation speed, diff pos baseline')
    subplot(2,1,2)
    scatter(playback.true_speed,run_speed_playback)
    xlim([-.5 3])
    ylim([-.5 3])
        title('playback')

  %%
  corrs_playback = zeros(1,max(playback.trial));
  difference_by_location_playback = zeros(max(playback.trial),100);
  difference_by_location_mismatch = zeros(max(mismatch.trial),100);

  pos_bins =0:2:400;
  pos_bins(1)=-Inf;
  pos_bins(end)=Inf;
  discrete_pos = discretize(playback.posx,pos_bins);
for ii = 1:max(playback.trial)
    idx = playback.trial==ii & playback.posx>2 & playback.posx<398;
    corrs_playback(ii)=corr(playback.true_speed(idx),[0;diff([playback.posx(idx)])]);
    tmp_d=abs(playback.true_speed(idx)-[0;diff(playback.posx(idx))]);
    
    
end
tmp_diff = playback.true_speed-[0;diff(playback.posx)];
for iT=1:max(playback.trial)
    for iP = 1:numel(pos_bins)
        idx = playback.trial == iT & discrete_pos == iP;
        difference_by_location_playback(iT,iP)=mean(abs(tmp_diff(idx)));
    end
end
figure
imagesc(difference_by_location_playback,[-2 2])
title('Difference run speed, visual speed')
tmp_diff = mismatch.true_speed-[0;diff(mismatch.posx)];
  discrete_pos = discretize(mismatch.posx,pos_bins);

for iT=1:max(mismatch.trial)
    for iP = 1:numel(pos_bins)
        idx = mismatch.trial == iT & discrete_pos == iP;
        difference_by_location_mismatch(iT,iP)=mean(abs(tmp_diff(idx)));
    end
end
figure
imagesc(difference_by_location_mismatch,[-2 2])
title('Difference run speed, visual speed baseline')
%% sort and rank trials from most similar to least similar
[a,b]=sort(corrs_playback,'descend');
ranking=1:numel(b);
ranking(b) = ranking;
figure
imagesc(difference_by_location_playback(b,:),[-2 2])
title('Difference run speed, visual speed')

%% plot for sanity check
figure
for iR=1:3
    idx = playback.trial == b(end-iR);
    %idxCL = mismatch.trial == b(iR);
    t_pb = playback.post(idx)-playback.post(find(idx,1));
    %t_fb = mismatch.post(idxCL)-mismatch.post(find(idxCL,1));
    subplot(2,3,iR)
    plot(t_pb,playback.true_speed(idx))
hold on
plot(t_pb,[0;diff(playback.posx(idx))])
legend({'run speed','vis speed'})
idx = playback.trial == b(iR);
    %idxCL = mismatch.trial == b(iR);
    t_pb = playback.post(idx)-playback.post(find(idx,1));
    ylim([0 3])

    subplot(2,3,iR+3)
    plot(t_pb,playback.true_speed(idx))
hold on
plot(t_pb,[0;diff(playback.posx(idx))])
ylim([0 3])
%plot(t_fb,[0;diff(mismatch.posx(idxCL))]);
end

%%
trials_ranked = ranking(playback.trial);
figure
scatter(playback.post,playback.posx,2,trials_ranked)
title('colored by similarity rank')
%%
%% combine rasters
% good_cells = mismatch.sp.cids(mismatch.sp.cgs==2);
figure;
for iC=1:length(good_cells)
spike_id = mismatch.sp.clu==good_cells(iC);
    spike_t = mismatch.sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,mismatch.post);
     
    scatter(mismatch.posx(spike_idx),mismatch.trial(spike_idx),2,'b')
    
    spike_id = playback.sp.clu==good_cells(iC);
    spike_t = playback.sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,playback.post);
     hold on
    scatter(playback.posx(spike_idx),trials_ranked(spike_idx)+max(mismatch.trial),2,'r')
    
    pause
    cla;
end
%%
pb_ranked = playback;
pb_ranked.trial = trials_ranked;
[sim,sim_s]=compute_spatial_similarity(pb_ranked);
%% top bottom 20%

%%
[~,spatialMap_smoothed] = get_spatial_map(playback);
[~,spatialMap_mm] = get_spatial_map(mismatch);
%% stability in baseline
tc_1=mean(spatialMap_mm(:,2:98,1:40),3);
tc_2=mean(spatialMap_mm(:,2:98,41:end),3);
tmp = corr(tc_1',tc_2');
stability=diag(tmp);
[~,stability_sort]=sort(stability);

figure
plot(tc_1(stability_sort(end-3),:))
hold on
plot(tc_2(stability_sort(end-3),:))

%%
[a,b]=sort(corrs_playback,'descend');
t_c_mm=mean(spatialMap_mm,3);
similar_idx = b(1:10);
diff_idx = b(end-9:end);
t_c_similar = mean(spatialMap_smoothed(:,:,similar_idx),3);
t_c_different = mean(spatialMap_smoothed(:,:,diff_idx),3);

% diff_similar = mean(t_c_similar-t_c_mm);
% diff_different = mean(t_c_different-t_c_mm);
% figure
% plot(diff_similar)
% hold ont
% plot(diff_different)
for ii=1:10
    cellIDX = stability_sort(end-ii+1);
    figure
    plot(t_c_mm(cellIDX,:))
    hold on
    plot(t_c_similar(cellIDX,:))
    plot(t_c_different(cellIDX,:))
    legend({'baseline','similar','different'})
end
%%
score_similar = zeros(1,numel(good_cells));
score_different = zeros(1,numel(good_cells));
for ii = 1:numel(good_cells)
    score_similar(ii)=corr(t_c_mm(ii,:)',t_c_similar(ii,:)');
    score_different(ii)=corr(t_c_mm(ii,:)',t_c_different(ii,:)');
end
figure
speed_difference_similar = nanmean(difference_by_location_playback(similar_idx,:));
speed_difference_different = nanmean(difference_by_location_playback(diff_idx,:));
subplot(1,2,1)
scatter(score_similar,score_different,4,stability)
xlabel('similar trials')
ylabel('different trials')
axis image
grid on
subplot(1,2,2)
plot(speed_difference_similar(3:end))
hold on
plot(speed_difference_different(3:end))
plot(nanmean(difference_by_location_mismatch(:,3:end)))
title('absolute difference in speed')
legend({'similar','different','baseline'})
    