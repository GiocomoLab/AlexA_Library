function [PEAKS,SHIFTS,corrMat,shiftMat_all,SPEED,SPEED2]=calculatePeakShiftSession_new(data,trials,ops)

good_cells = data.sp.cids(data.sp.cgs==2);
[corrMat,frMat,shiftMat_all]=trialCorrMat(good_cells,trials,data,ops);
speed = calcSpeed(data.posx,ops);
fast_idx = speed >= ops.SpeedCutoff;



template1_trials = 1:ops.num_tr_bl;
template1={};
cntr = 0;



for ii=template1_trials
    cntr = cntr+1;
    tmp_idx =setdiff(template1_trials,ii);

    tmp = squeeze(nanmean(frMat(:,tmp_idx,:),2));
    template1{cntr}=tmp;
end
template1{numel(template1_trials)+1}=squeeze(nanmean(frMat(:,template1_trials,:),2));


%%
nBins = ops.nBins;
startVec = ops.stride_start:ops.stride:(nBins-ops.chunksize+1);
nReps = numel(startVec);
maxlag = ops.max_lag/ops.SpatialBin;
nTrials = ops.num_tr_tot;
nCells = numel(good_cells);
SHIFTS = nan(nCells,nTrials,nReps);
PEAKS = SHIFTS;
SPEED = nan(nTrials,nReps);
repidx = 0;

trial2templateMap=zeros(1,nTrials);
cntr=1;
aa=numel(template1_trials)+1;
for iT = 1:numel(trial2templateMap)
    if ~ismember(iT,template1_trials)
        trial2templateMap(iT)=aa;
    else
        trial2templateMap(iT)=cntr;
        cntr=cntr+1;
    end
end
    shift_all = -ops.max_lag:ops.SpatialBin:ops.max_lag;
for iStart = startVec
    repidx = repidx+1;
    startbin = iStart;
    stopbin = iStart+ops.chunksize-1;
    %xcorr relative to template trials: 4 before gain onset for all except
    %these 4
    [corrMat_snippet,shiftMat] = calc_xcorr_snippet(frMat,template1,startbin,stopbin,maxlag,trial2templateMap);
    %[xcorrs2,~] = calc_xcorr_snippet(spatialMap,template2,startbin,stopbin,maxlag);
    sm=shift_all(shiftMat);
    sm(isnan(corrMat_snippet))=nan;
    
    PEAKS(:,:,repidx)=corrMat_snippet;
    SHIFTS(:,:,repidx)=sm;
    for iT=1:numel(trials)
        t_idx = data.trial ==trials(iT);
        p_idx = data.posx>=ops.SpatialBin*startbin & data.posx<=stopbin*ops.SpatialBin;
        SPEED(iT,repidx)=mean(speed(t_idx & p_idx));
        SPEED2(iT,repidx)=mean(speed(t_idx & p_idx & fast_idx));
    end
end
