files = dir('Z:\giocomo\attialex\matfiles_new\*sine*.mat')
opt = load_mismatch_opt;
%%
MM=[];
MMSquare=[];
SID=[];
SPIKE_TIMES=cell(numel(files),1);
RUN_TRACES = cell(numel(files),1);
CLUIDS = RUN_TRACES;
for iF=1:numel(files)
    sine_dat = load(fullfile(files(iF).folder,files(iF).name));
    square_name = strrep(files(iF).name,'Sine','Square');
    if ~isfile(fullfile(files(iF).folder,square_name))
        continue
    end
    square_dat = load(fullfile(files(iF).folder,square_name));
    
    
    
    mismatch_trigger = sine_dat.vr_data_resampled.Move<0.5;
    mm_trig_square = square_dat.vr_data_resampled.Move<0.5;
    good_cells = sine_dat.sp.cids(sine_dat.sp.cgs==2);
    
    all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
    true_speed = sine_dat.vr_data_resampled.Speed;
    if iscolumn(true_speed)
        speed=true_speed';
    else
        speed=true_speed;
    end
    filt = gausswin(61); %61 pretty close to what we use in other
    filt = filt/sum(filt);
    smooth_speed = conv(speed,filt,'same');
    run_periods=smooth_speed>opt.speed_t;
    run_window=-30:30;
    possibles=strfind(run_periods,ones(1,length(run_window)))+floor(.5*length(run_window));
    
    mm_trigs=all_mm_trigs(ismember(all_mm_trigs,possibles));
    
    
    all_mm_trigs_square=strfind(mm_trig_square>0.9,[0 0 1 1])+2;
    true_speed_square = square_dat.vr_data_resampled.Speed;
    
    smooth_speed_square = conv(true_speed_square,filt,'same');
    run_periods_square=smooth_speed_square>opt.speed_t;
    run_window=-30:30;
    possibles_square=strfind(run_periods_square,ones(1,length(run_window)))+floor(.5*length(run_window));
    
    mm_trigs_square=all_mm_trigs_square(ismember(all_mm_trigs_square,possibles_square));
    
    [spikeTimes,~,aux]=extract_triggered_spikeTimes(sine_dat.sp,sine_dat.post(mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[sine_dat.post' ;smooth_speed],'aux_win',opt.aux_win);
    
    [spikeTimesSquare]=extract_triggered_spikeTimes(square_dat.sp,square_dat.post(mm_trigs_square),'cluIDs',good_cells,'win',opt.extract_win);
    
    trial_vec =cat(1,spikeTimes{:});
    trial_vec_square = cat(1,spikeTimesSquare{:});
    CLUIDS{iF}=good_cells;
    count_vec = zeros(numel(good_cells),numel(opt.time_bins)-1);
    count_vec_square = count_vec;
    for iC=1:numel(good_cells)
        idx = trial_vec(:,2)==good_cells(iC);
        [spike_count]=histcounts(trial_vec(idx,1),opt.time_bins);
        count_vec(iC,:)=spike_count;
        idx = trial_vec_square(:,2)==good_cells(iC);
        [spike_count]=histcounts(trial_vec_square(idx,1),opt.time_bins);
        count_vec_square(iC,:)=spike_count;
        
    end
    n_trigs_included = numel(unique(trial_vec(:,3)));
    n_trigs_included_square = numel(unique(trial_vec_square(:,3)));
    count_vec = count_vec/n_trigs_included/opt.TimeBin;
    count_vec_square = count_vec_square/n_trigs_included_square/opt.TimeBin;
    
    MM=cat(1,MM,count_vec);
    SID = cat(1,SID,ones(numel(good_cells),1)*iF);
    MMSquare=cat(1,MMSquare,count_vec_square);
end
