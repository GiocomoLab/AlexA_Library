function opt = load_mismatch_opt()

opt = struct();
opt.speed_t=0.05;
opt.extract_win = [-2 3];
opt.aux_win = [-50 50];
opt.TimeBin = 0.02;
opt.time_bins =-2:0.02:3;
opt.time_vecs = opt.time_bins(1:end-1)*0.5+opt.time_bins(2:end)*0.5;
opt.extract_win = [-2 3];
opt.aux_win = [-50 50];
opt.TimeBin = 0.02;
opt.speed_filt_win = 61
opt.speedSigma = 10;
opt.smoothSigma_time = 0.2; % in sec; for smoothing fr vs time
opt.smoothSigma_dist = 5; % in cm; for smoothing fr vs distance
opt.SpatialBin = 2;
opt.run_window = -30:30;

% opt = struct();
% 
% opt.TimeBin = 0.02;
% opt.SpatialBin = 2; % cm
% opt.SmoothSigmaFR = 4; % cm
opt.TrackStart = 0; % cm
opt.TrackEnd = 400; % cm
opt.SpeedCutoff = 2; % cm/s
% opt.cell_info_file = 'cell_info';
% opt.brain_region = 'MEC';
% opt.save_figs = true;
% opt.save_results = true;
% opt.stab_thresh = 0.5;
% opt.max_num_cells = 537; % max num cells in a session
% opt.xbinedges = opt.TrackStart:opt.SpatialBin:opt.TrackEnd;
% opt.xbincent = opt.xbinedges(1:end-1)+opt.SpatialBin/2;
opt.track_length = opt.TrackEnd-opt.TrackStart;
% opt.gains_all = [0.8 0.7 0.6 0.5];
% opt.contr_all = [100 50 20 10 5 2 0];
% opt.max_lag = 30; % cm; for xcorr peak
% opt.num_tr_gc = 4;
% opt.num_tr_bl = 6;
% opt.num_tr_tot = opt.num_tr_gc + 2*opt.num_tr_bl;
% opt.xcorr_chunk = 200; % cm
% opt.num_pos_chunks = opt.track_length/opt.xcorr_chunk;
% opt.pos_chunk_idx = 1:numel(opt.xbincent);
% opt.pos_chunk_idx = reshape(opt.pos_chunk_idx,numel(opt.xbincent)/opt.num_pos_chunks,opt.num_pos_chunks)';
% % opt.num_lags = 2*opt.max_lag/opt.SpatialBin+1;
% opt.num_chunk_tot = opt.num_tr_tot*opt.num_pos_chunks;
% opt.tol = 0.0001; % for comparing floats
% opt.trials_per_contrast = 10; % trials per contrast in range of contrasts sessions
% opt.smoothSigma_dist = 5; % in cm; for smoothing fr vs distance
% opt.smoothSigma_time = 0.2; % in sec; for smoothing fr vs time
% 
% % for rep clustering
% opt.min_num_stab_cells = 5;
% opt.min_frac_stab_cells = 0; % had 0.2 before, decided to remove this requirement for simplicity
% opt.num_clust = 3; % for k means clustering
% opt.num_pcs_clustering = 3; % num PCs to keep for k means clustering
% 
% [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
% opt.colors.MEC=cb(1,:);
% opt.colors.VISp = cb(4,:);
% opt.colors.RS = cb(5,:);
