function [corrMat,frMat,shiftMat] = trialCorrMat(cell_id,trials,dat,opt)
% takes max over lags from -opt.max_lag to +opt.max_lag (in cm)
% MGC 12/15/2019

corrMat = nan(numel(cell_id),numel(trials),numel(trials));
shiftMat = nan(numel(cell_id),numel(trials),numel(trials));
max_lag = getOr(opt,'max_lag',0);
shift_all = -max_lag:opt.SpatialBin:max_lag;

frMat = calcTrialFRMat(cell_id,trials,dat,opt); % single trial fr mat
for i = 1:numel(cell_id)
    if numel(cell_id)==1
        fr_this = frMat';
    else
        fr_this = squeeze(frMat(i,:,:))';
    end
    fr_this = fr_this-nanmean(fr_this); % subtract mean on each trial
    xcorr_this = xcorr(fr_this,max_lag/opt.SpatialBin,'coeff');
    [xcorr_this,max_idx] = max(xcorr_this,[],1); % take max over lags
    shift_this = shift_all(max_idx);
    shift_this(isnan(xcorr_this)) = nan;
    xcorr_this = reshape(xcorr_this,numel(trials),numel(trials));
    shiftMat(i,:,:) = reshape(shift_this,numel(trials),numel(trials)); 
    corrMat(i,:,:) = xcorr_this+diag(nan(numel(trials),1)); % make diagonals nan
end
corrMat = squeeze(corrMat);

end