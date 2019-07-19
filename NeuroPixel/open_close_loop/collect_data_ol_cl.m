good_cells = mismatch.sp.cids(mismatch.sp.cgs==2);
pos_bins = [0:4:400];
pos_bins(1)=-0.01;
pos_bins(end)=404;
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


tmp_diff = playback.true_speed-[0;diff(playback.posx)];
  discrete_pos = discretize(playback.posx,pos_bins);

for iT=1:max(playback.trial)
    for iP = 1:numel(pos_bins)
        idx = playback.trial == iT & discrete_pos == iP;
        difference_by_location_playback(iT,iP)=mean(abs(tmp_diff(idx)));
    end
end



[~,spatialMap_mm] = get_spatial_map(mismatch);
[~,spatialMap_pb] = get_spatial_map(playback);

difference_by_position = zeros(max(mismatch.trial),400);
difference_by_position_time = zeros(max(mismatch.trial),400);

for iT=1:max(playback.trial)
    idx = playback.trial==iT;
    loc_diff = playback.posx(idx)-cumsum(playback.true_speed(idx));
    loc_time = playback.post(idx);
    loc_time = loc_time-loc_time(1);
    interp = interp1(loc_time,loc_diff,linspace(0,loc_time(end),400));
    difference_by_position(iT,:)=interp;
    difference_by_position_time(iT,:)=linspace(0,loc_time(end),400);
end




