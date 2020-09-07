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