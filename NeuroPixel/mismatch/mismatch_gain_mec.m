%get mismatch files
mm_files = dir('/Volumes/T7/attialex/NP_DATA_corrected/*_mismatch_*.mat');
gain_path = '/Volumes/T7/attialex/NP_DATA_corrected';

opt = load_mismatch_opt;
gain_opt = load_default_opt;
%for each mismatch file, find gain file
gain_files = cell(numel(mm_files),1);
for iF=1:numel(mm_files)
    locs = strfind(mm_files(iF).name,'_');
    animal_date = mm_files(iF).name(1:locs(2));
    
    gain_file = dir(fullfile(gain_path,[animal_date 'g*.mat']));
    if isempty(gain_file)
        gain_file = dir(fullfile(gain_path,[animal_date 'c*.mat']));
    end
    
    if isempty(gain_file)
        disp(mm_files(iF).name)
    else
        gain_files{iF}=gain_file(1).name;
    end
    %disp(numel(gain_file))
end

%loop through pairs,
for iF=1:numel(mm_files)
    if isempty(gain_files{iF})
        continue
    end
    mm_data = load(fullfile(mm_files(iF).folder,mm_files(iF).name));
    gain_data = load(fullfile(gain_path,gain_files{iF}));
    
    
    if ~isfield(gain_data,'anatomy')
        disp(mm_files(iF).name)
        disp('no anatomy')
        nC= numel(gain_data.sp.cgs);
        reg = repmat({'UNK'},1,nC);
        depth = nan(nC,1); % go via channel number
        %continue
        valid_region = true(size(reg));
    else
        
        if ~isfield(gain_data.anatomy,'depth')
            dd =gain_data.anatomy.z2-gain_data.anatomy.tip_distance;
            gain_data.anatomy.depth = dd';
        end
        if isfield(gain_data.anatomy,'parent_shifted')
            reg = gain_data.anatomy.parent_shifted;
            depth = gain_data.anatomy.depth_shifted';
        else
            reg = gain_data.anatomy.cluster_parent;
            depth = gain_data.anatomy.depth';
        end
        if iscolumn(reg)
            reg = reg';
        end
        %valid_region = startsWith(reg,{'MEC','ECT'});
        valid_region= true(size(reg));
    end
    if nnz(valid_region)==0
        continue
    end
    
    
    mismatch_trigger = mm_data.mismatch_trigger;
    if size(mismatch_trigger,1) ~=1
        mismatch_trigger=mismatch_trigger';
    end
    if nnz(mismatch_trigger==1)>nnz(mismatch_trigger==0)
        %in some old files mismatch_trigger was actually the move variable,
        %i.e. move ==0 is mismatch
        mismatch_trigger = mismatch_trigger<0.1;
    end
    
    good_cells = gain_data.sp.cids(gain_data.sp.cgs==2 & valid_region);
    depth_this = depth(gain_data.sp.cgs==2 & valid_region);
    region_this = reg(gain_data.sp.cgs==2 & valid_region);
    all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
    true_speed = mm_data.true_speed;
    %out = fitGLM_mismatch(data,good_cells);
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
    tmp = run_periods & mismatch_trigger<0.9;
    possibles_random = strfind(tmp,ones(1,length(run_window)))+floor(.5*length(run_window));
    possibles_random = randsample(possibles_random,numel(mm_trigs));
    if all(size(gain_data.post)==size(speed))
        speed = speed';
    end
    %[~,sp] = calcSpeed(gain_data.posx,load_default_opt);
    [spikeTimes,~,aux,~,count_vec]=extract_triggered_spikeTimes(mm_data.sp,mm_data.post(mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[mm_data.post' ;smooth_speed],'aux_win',opt.aux_win);
    
    firing_rate_mm = zeros(numel(good_cells),1);
    for iC=1:numel(good_cells)
        idx = mm_data.sp.clu==good_cells(iC);
        firing_rate_mm(iC)=nnz(idx);
    end
    firing_rate_mm = firing_rate_mm/(mm_data.post(end)-mm_data.post(1));
    
    firing_rate_gain = zeros(numel(good_cells),1);
    for iC=1:numel(good_cells)
        idx = gain_data.sp.clu==good_cells(iC);
        firing_rate_gain(iC)=nnz(idx);
    end
    firing_rate_gain = firing_rate_gain/(gain_data.post(end)-gain_data.post(1));
    
    
    mm_resp = mean(count_vec(:,opt.time_vecs>0.1 & opt.time_vecs<0.6),2)-mean(count_vec(:,opt.time_vecs<-.1 & opt.time_vecs>-.6),2);
    mm_resp = mm_resp./firing_rate_mm;
    count_vecN = smoothdata(count_vec-mean(count_vec(:,opt.time_vecs<-.1 & opt.time_vecs>-.6),2),2,'gaussian',5)./firing_rate_mm;
    [~,sid]=sort(mm_resp,'descend','MissingPlacement','last');
    
    [corr_mat,fr_mat,shift_mat]=trialCorrMat(good_cells,1:max(gain_data.trial),gain_data,gain_opt);
    
    g_ons = strfind(gain_data.trial_gain' < 1 & gain_data.trial_contrast' == 100 ,[0 1]);
    
    pre_mat = fr_mat(:,g_ons,:);
    post_mat = fr_mat(:,g_ons+1,:);
    pm=squeeze(mean(pre_mat,2));
    pom = squeeze(mean(post_mat,2));
    figure('Name',mm_files(iF).name);
    subplot(3,1,1)
    plot(mean(pm(sid(1:20),:)),'b')
    hold on
    plot(mean(pom(sid(1:20),:)),'r')
    plot(nanmean(pm(sid(21:end),:)),'b--')
    plot(nanmean(pom(sid(21:end),:)),'r--')
    title('spatial firing rate')
    
    gain_speed =calcSpeed(gain_data.posx,gain_opt);
    gain_speed = gain_speed./gain_data.trial_gain(gain_data.trial);
    
    run_ons = strfind(gain_speed'>1,[zeros(1,30),ones(1,50)])+30;
    all_trials = 1:numel(gain_data.trial_gain);
    bl_trials = all_trials(gain_data.trial_contrast==100 & gain_data.trial_gain==1);
    bl_onset = run_ons(ismember(gain_data.trial(run_ons),bl_trials));
    gain_trials = all_trials(gain_data.trial_contrast ==100 & gain_data.trial_gain<1);
    gain_onset = run_ons(ismember(gain_data.trial(run_ons),gain_trials));
    
    [spikeTimes,~,aux,~,count_vec_run_bl]=extract_triggered_spikeTimes(gain_data.sp,gain_data.post(bl_onset),'cluIDs',good_cells,'win',opt.extract_win);
    if ~isempty(gain_onset)
        [spikeTimes,~,aux,~,count_vec_run_gain]=extract_triggered_spikeTimes(gain_data.sp,gain_data.post(gain_onset),'cluIDs',good_cells,'win',opt.extract_win);
    else
        count_vec_run_gain = nan(size(count_vec_run_bl));
    end
    
    subplot(3,1,2)
    plot(opt.time_vecs,mean(smoothdata(count_vec_run_bl(sid(1:20),:),2,'gaussian',10)))
    hold on
    plot(opt.time_vecs,mean(smoothdata(count_vec_run_gain(sid(1:20),:),2,'gaussian',10)))
    plot(opt.time_vecs,mean(smoothdata(count_vecN(sid(1:20),:),2,'gaussian',10)))
    title('running onsets')
    legend({'during baseline','during gain'})
    last_bl_trials = strfind(gain_data.trial_contrast' == 100 & gain_data.trial_gain' < 1,[0 1]);
    gain_onsets = strfind(gain_data.trial_contrast(gain_data.trial)'==100 & gain_data.trial_gain(gain_data.trial)'<1,[0 1]);
    
    last_bl_onsets = strfind(ismember(gain_data.trial,last_bl_trials)',[0 1]);
    
    [spikeTimes,~,aux,~,count_vec_bl_onsets]=extract_triggered_spikeTimes(gain_data.sp,gain_data.post(last_bl_onsets),'cluIDs',good_cells,'win',opt.extract_win);
    [spikeTimes,~,aux,~,count_vec_gain_onsets]=extract_triggered_spikeTimes(gain_data.sp,gain_data.post(gain_onsets),'cluIDs',good_cells,'win',opt.extract_win);
    subplot(3,1,3)
    count_vec_bl_onsets=count_vec_bl_onsets-mean(count_vec_bl_onsets(:,opt.time_vecs>-.6 & opt.time_vecs<-.1),2);
        count_vec_gain_onsets=count_vec_gain_onsets-mean(count_vec_gain_onsets(:,opt.time_vecs>-.6 & opt.time_vecs<-.1),2);

    plot(opt.time_vecs,mean(smoothdata(count_vec_bl_onsets(sid(1:20),:),2,'gaussian',10)))
    hold on
    plot(opt.time_vecs,mean(smoothdata(count_vec_gain_onsets(sid(1:20),:),2,'gaussian',10)))
    plot(opt.time_vecs,nanmean(smoothdata(count_vecN(sid(1:20),:),2,'gaussian',10)))
    title('firing rate beginning of trial')
    legend({'baseline pre','first gain trial'})
    drawnow
end
%get spatial firing rates

%compare