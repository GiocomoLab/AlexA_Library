matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA_2/*mismatch*.mat');

matfiles = matfiles(~cellfun(@(x) contains(x,'tower'), {matfiles.name}));
opt = load_mismatch_opt;
opt.time_bins =-2:0.02:3;
opt.extract_win = [-2 3];
MM=[];
avgMM = [];
avgMMR = [];
SID = [];
MM_R=[];
for iF=1:numel(matfiles)
    disp(iF)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    mismatch_trigger = data_out.mismatch_trigger;
    if size(mismatch_trigger,1) ~=1
        mismatch_trigger=mismatch_trigger';
    end
    if nnz(mismatch_trigger==1)>nnz(mismatch_trigger==0)
        %in some old files mismatch_trigger was actually the move variable,
        %i.e. move ==0 is mismatch
        mismatch_trigger = mismatch_trigger<0.1;
    end
    
    
    all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
    true_speed = data_out.true_speed;
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
    if all(size(data_out.post)==size(speed))
        speed = speed';
    end
    %[~,sp] = calcSpeed(data_out.posx,load_default_opt);
    %% MM
    %[spike_mat,win,adata]=extract_triggered_spikes(data_out.sp,data_out.post(mm_trigs),'win',[-4 4],'aux',[data_out.post'; [speed];sp'],'aux_win',[-200 200]);
    %[spike_mat_random,~,adata_random]=extract_triggered_spikes(data_out.sp,data_out.post(possibles_random),'win',[-4 4],'aux',[data_out.post'; [speed]],'aux_win',[-200 200]);
    [spikeTimes]=extract_triggered_spikeTimes(data_out.sp,data_out.post(mm_trigs),'win',opt.extract_win);
    [spikeTimesRandom]=extract_triggered_spikeTimes(data_out.sp,data_out.post(possibles_random),'win',opt.extract_win);
    
    trial_vec =cat(1,spikeTimes{:});
    trial_vec_random = cat(1,spikeTimesRandom{:});
    good_cells = data_out.sp.cids(data_out.sp.cgs==2);
    count_vec = zeros(numel(good_cells),numel(opt.time_bins)-1);
    count_vec_random = count_vec;
    for iC=1:numel(good_cells)
        idx = trial_vec(:,2)==good_cells(iC);
        [spike_count]=histcounts(trial_vec(idx,1),opt.time_bins);
        count_vec(iC,:)=spike_count;
         idx = trial_vec_random(:,2)==good_cells(iC);
        [spike_count]=histcounts(trial_vec_random(idx,1),opt.time_bins);
        count_vec_random(iC,:)=spike_count;      
        
    end
    n_trigs_included = numel(unique(trial_vec(:,3)));
        n_trigs_included_random = numel(unique(trial_vec_random(:,3)));

    count_vec = count_vec/n_trigs_included;
    count_vec_random = count_vec_random/n_trigs_included_random;
    
    MM=cat(1,MM,count_vec);
    MM_R = cat(1,MM_R,count_vec_random);
    SID = cat(1,SID,ones(numel(good_cells),1)*iF);
    avgMM = cat(1,avgMM,mean(count_vec));
    avgMMR = cat(1,avgMMR,mean(count_vec_random));
    %MMR=cat(1,MMR,count_vec_random);
end
%%
figure
MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
imagesc((MM_ms),[-1 1])
cmap = cbrewer('div','RdBu',20);
cmap=flipud(cmap);
colormap(cmap);
%%
params=struct();
params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
params.xLim=[-2 3];
figure
plotAVGSEM(MM',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0)
plotAVGSEM(MM_R',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[.5 .5 .5])

%%
