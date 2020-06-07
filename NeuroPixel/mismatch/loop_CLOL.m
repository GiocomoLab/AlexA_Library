opt = load_mismatch_opt;
opt.time_bins =-2:0.02:3;
opt.time_vecs = opt.time_bins(1:end-1)*0.5+opt.time_bins(2:end)*0.5;
opt.extract_win = [-2 3];
opt.TimeBin = 0.02;
opt.smoothSigma_time = 0.0; % in sec; for smoothing fr vs time
pb_files = dir('/Users/attialex/NP_DATA_2/*_playback_*.mat');
MM=[];
PB=[];
PB_slow=[];
SID=[];
for iF=6%1:numel(pb_files)
        mm_name = strrep(pb_files(iF).name,'playback','mismatch');
    if ~isfile(fullfile(pb_files(iF).folder,mm_name))
        continue
    end
    data_pb = load(fullfile(pb_files(iF).folder,pb_files(iF).name));
    data_mm = load(fullfile(pb_files(iF).folder,mm_name));
        mismatch_trigger = data_mm.mismatch_trigger;

    if size(mismatch_trigger,1) ~=1
        mismatch_trigger=mismatch_trigger';
    end
    if nnz(mismatch_trigger==1)>nnz(mismatch_trigger==0)
        %in some old files mismatch_trigger was actually the move variable,
        %i.e. move ==0 is mismatch
        mismatch_trigger = mismatch_trigger<0.1;
    end
    
    
    all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
    true_speed = data_mm.true_speed;
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
    [~,sp] = calcSpeed(data_mm.posx,load_default_opt);
    [spikeTimesMM,~,aux]=extract_triggered_spikeTimes(data_mm.sp,data_mm.post(mm_trigs),'win',opt.extract_win,'aux',[data_mm.post' ;sp'],'aux_win',[-50 50]);
    mm_trigs_pb = zeros(size(mm_trigs));
    for iT=1:numel(mm_trigs)
    [~,mm_trigs_pb(iT)] = min(abs(data_pb.post-data_mm.post(mm_trigs(iT))));
    end
    
    true_speed = data_pb.true_speed;
    if iscolumn(true_speed)
        speed=true_speed';
    else
        speed=true_speed;
    end
    speed = conv(speed,filt,'same');
    
    [~,sp] = calcSpeed(data_pb.posx,load_default_opt);
    [spikeTimesPB,~,aux]=extract_triggered_spikeTimes(data_pb.sp,data_pb.post(mm_trigs_pb),'win',opt.extract_win,'aux',[data_pb.post';sp';speed],'aux_win',[-50 50]);
    run_speed=squeeze(mean(aux(2,:,25:75),3));
    [~,sidx_pb]=sort(run_speed);
    slow_trials=sidx_pb(1:floor(numel(mm_trigs)/3));
    
    trial_vec =cat(1,spikeTimesMM{:});
    trial_vec_pb = cat(1,spikeTimesPB{:});
    idx_slow = ismember(trial_vec_pb(:,3),slow_trials);
    
    trial_rank = 1:max(trial_vec_pb(:,3));
    trial_rank(sidx_pb)=trial_rank;
    trial_vec_pb_sorted = trial_rank(trial_vec_pb(:,3));

    
    
    good_cells = data_mm.sp.cids(data_mm.sp.cgs==2);
    count_vec = zeros(numel(good_cells),numel(opt.time_bins)-1);
    count_vec_pb = count_vec;
    count_vec_pb_slow = count_vec;
    
    for iC=1:numel(good_cells)
        idx = trial_vec(:,2)==good_cells(iC);
        [spike_count]=histcounts(trial_vec(idx,1),opt.time_bins);
        count_vec(iC,:)=spike_count;
        idx = trial_vec_pb(:,2)==good_cells(iC);
        [spike_count]=histcounts(trial_vec_pb(idx,1),opt.time_bins);
        count_vec_pb(iC,:)=spike_count;
        spike_count_slow = histcounts(trial_vec_pb(idx & idx_slow,1),opt.time_bins);
        count_vec_pb_slow(iC,:)=spike_count_slow;
    end
    n_trigs_included = numel(unique(trial_vec(:,3)));
    n_trigs_included_pb = numel(unique(trial_vec_pb(:,3)));
    
    count_vec = count_vec/n_trigs_included/opt.TimeBin;
    count_vec_pb = count_vec_pb/n_trigs_included_pb/opt.TimeBin;

    MM=cat(1,MM,count_vec);
    PB = cat(1,PB,count_vec_pb);
    PB_slow = cat(1,PB_slow,count_vec_pb_slow);
    SID = cat(1,SID,ones(numel(good_cells),1)*iF);
    figure
    mm_resp = mean(count_vec(:,105:130),2)-mean(count_vec(:,75:100),2);
    [~,cell_sort]=sort(mm_resp,'descend');
    for iC=1:10
    subplot(1,2,1)
            idx = trial_vec(:,2)==good_cells(cell_sort(iC));
            scatter(trial_vec(idx,1),trial_vec(idx,3),15,'.')
            subplot(1,2,2)
                        idx = trial_vec_pb(:,2)==good_cells(cell_sort(iC));

            scatter(trial_vec_pb(idx,1),trial_vec_pb_sorted(idx),15,'.')
        pause
        clf
    end
end
%%
mm_resp = mean(count_vec(:,105:130),2)-mean(count_vec(:,75:100),2);
[~,sidx]=sort(mm_resp,'desc');
figure
for iC=1:numel(sidx)
    plot(count_vec(sidx(iC),:))
    hold on
    plot(count_vec_pb(sidx(iC),:));
    plot(count_vec_pb_slow(sidx(iC),:));
    pause
    cla
end