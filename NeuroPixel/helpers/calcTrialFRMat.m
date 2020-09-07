function frMat = calcTrialFRMat(cell_id,trials,dat,opt)
% computes single trial firing rate matrix for single cell
% MGC 12/15/2019
% edit aa: removed calculation for contiguous trials as this interfered
% with shifted posx

num_xbin = opt.track_length/opt.SpatialBin;

% contiguous trials
try
    % contiguous trials
if max(diff(trials))==1
    frContig = calcFRVsDist(cell_id,trials,dat,opt);
    frMat = reshape(frContig,numel(cell_id),num_xbin,numel(trials));
    frMat = permute(frMat,[1 3 2]);
else % go trial by trial
    frMat = nan(numel(cell_id),numel(trials),num_xbin);
    for tr_idx = 1:numel(trials)
        frMat(tr_idx,:) = calcFR(cell_id,trials(tr_idx),dat,opt);
    end
end
    
catch

    frMat = nan(numel(cell_id),numel(trials),num_xbin);
    for tr_idx = 1:numel(trials)
        frMat(:,tr_idx,:) = calcFR(cell_id,trials(tr_idx),dat,opt);
    end
end

frMat = squeeze(frMat);

end