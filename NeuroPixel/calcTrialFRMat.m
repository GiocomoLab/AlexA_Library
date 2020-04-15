function frMat = calcTrialFRMat(cell_id,trials,dat,opt)
% computes single trial firing rate matrix for single cell
% MGC 12/15/2019
% edit aa: removed calculation for contiguous trials as this interfered
% with shifted posx

num_xbin = opt.track_length/opt.SpatialBin;

% contiguous trials

    frMat = nan(numel(cell_id),numel(trials),num_xbin);
    for tr_idx = 1:numel(trials)
        frMat(:,tr_idx,:) = calcFR(cell_id,trials(tr_idx),dat,opt);
    end


frMat = squeeze(frMat);

end