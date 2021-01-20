%get mismatch files
mm_files = dir('/Volumes/T7/attialex/NP_DATA_corrected/*_mismatch_*.mat');
delay_path = '/Volumes/T7/attialex/speed_filtered_correctedData_shortidx2';

opt = load_mismatch_opt;
gain_opt = load_default_opt;
%for each mismatch file, find gain file
delay_files = cell(numel(mm_files),1);
for iF=1:numel(mm_files)
    locs = strfind(mm_files(iF).name,'_');
    animal_date = mm_files(iF).name(1:locs(2));
    
    gain_file = dir(fullfile(delay_path,[animal_date 'g*.mat']));
    if isempty(gain_file)
        gain_file = dir(fullfile(delay_path,[animal_date 'c*.mat']));
    end
    
    if isempty(gain_file)
        disp(mm_files(iF).name)
    else
        delay_files{iF}=gain_file(1).name;
    end
    %disp(numel(gain_file))
end
%%
DELAY = [];
MM = [];
RUN = [];
FIRING_RATE = [];
SID = [];
REGION={};
%loop through pairs,
for iF=1:numel(mm_files)
    if isempty(delay_files{iF})
        continue
    end
    mm_data = load(fullfile(mm_files(iF).folder,mm_files(iF).name));
    delay_data = load(fullfile(delay_path,delay_files{iF}));
    
    
    if ~isfield(delay_data,'region')
        disp(mm_files(iF).name)
        disp('no anatomy')
        nC= numel(delay_data.sp.cgs);
        reg = repmat({'UNK'},1,nC);
        depth = nan(nC,1); % go via channel number
        %continue
        valid_region = true(size(reg));
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
    
    good_cells = mm_data.sp.cids(mm_data.sp.cgs==2);
    if numel(good_cells) ~= numel(delay_data.CID)
        error('number of clusters does not match')
    end
    
    all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
    true_speed = mm_data.true_speed;
    
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
    
    %[~,sp] = calcSpeed(gain_data.posx,load_default_opt);
    [spikeTimes,~,aux,~,count_vec]=extract_triggered_spikeTimes(mm_data.sp,mm_data.post(mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[mm_data.post' ;smooth_speed],'aux_win',opt.aux_win);
    
    firing_rate_mm = zeros(numel(good_cells),1);
    for iC=1:numel(good_cells)
        idx = mm_data.sp.clu==good_cells(iC);
        firing_rate_mm(iC)=nnz(idx);
    end
    firing_rate_mm = firing_rate_mm/(mm_data.post(end)-mm_data.post(1));
    
    
    
    
    run_ons = strfind(smooth_speed>opt.speed_t,[zeros(1,30),ones(1,50)])+30;
    
    [spikeTimes,~,aux,~,count_vec_run]=extract_triggered_spikeTimes(mm_data.sp,mm_data.post(run_ons),'cluIDs',good_cells,'win',opt.extract_win,'aux',[mm_data.post' ;smooth_speed],'aux_win',opt.aux_win);
    
    stab = delay_data.all_stability;
    stable_blocks = stab>=.5;
    tmp_factors = delay_data.all_factors;
    tmp_factors(~stable_blocks)=nan; %set shifts to nan where stability <.5
    
    enough_data_idx = sum(~isnan(tmp_factors),1)>2; %only use data where there are more than 2 stable blocks
    
    dat = nanmean(tmp_factors,1)';
    dat(~enough_data_idx)=nan;
    
    DELAY = cat(1,DELAY,dat);
    MM = cat(1,MM,count_vec);
    RUN = cat(1,RUN,count_vec_run);
    FIRING_RATE = cat(1,FIRING_RATE,firing_rate_mm);
    SID = cat(1,SID,iF*ones(size(firing_rate_mm)));
    REGION = cat(1,REGION,delay_data.region');
end
%get spatial firing rates

%%
rr={'VISp','RS','MEC'};
for iR=1:numel(rr)
    valid = ~isnan(DELAY) & startsWith(REGION,rr{iR}) & FIRING_RATE>1;
    [~,dd]=sort(DELAY(valid));
    count_vec_n = MM(valid,:)-mean(MM(valid,opt.time_vecs>-.6 & opt.time_vecs<-.1),2);
    count_vec_n = smoothdata(count_vec_n,2,'gaussian',3)./FIRING_RATE(valid);
    figure
    subplot(2,3,[1 2 4 5])
    imagesc(opt.time_vecs,1:nnz(valid),count_vec_n(dd,:),[-2 2])
    xlim([-.5 2.5])
    xlabel('time from onset')
        title(rr{iR})

    subplot(2,3,[3 6])
    DELAY_valid = DELAY(valid);
    plot(DELAY_valid(dd),1:nnz(valid),'.')
    set(gca,'YDir','reverse')
    ylim([0 nnz(valid)])
    colormap(brewermap(14,'*RdYlGn'))
    xlabel('delay')
end
%%
rr={'VISp','RS','MEC'};
figure
for iR=1:numel(rr)
    subplot(1,numel(rr),iR)
    hold on
    quantiles = 3;
    valid = ~isnan(DELAY) & startsWith(REGION,rr{iR}) & FIRING_RATE>1;
    
    nSamples = nnz(valid);
    %cmap = cbrewer('div','RdBu');
    cmap = brewermap(quantiles,'Set1');
    xl =[-1 2.5];
    [~,sid]=sort(DELAY(valid));
    MM_smooth = MM-mean(MM(:,opt.time_vecs>-.6 & opt.time_vecs<-.1),2);
    MM_smooth = smoothdata(MM_smooth,2,'gaussian',7)./FIRING_RATE(:);
    MM_valid = MM_smooth(valid,:);
    MM_other = MM_smooth(isnan(DELAY) & startsWith(REGION,rr{iR}) & FIRING_RATE>1,:);
    DELAY_valid = DELAY(valid);
    for iC=1:quantiles%nChunks
        sub_idx = round((iC-1)/quantiles*nSamples + 1):round(iC/quantiles*nSamples);
        IDX = sid(sub_idx);
        boundedline(opt.time_vecs,nanmean(MM_valid(IDX,:))*100,nanstd(MM_valid(IDX,:)*100/sqrt(nnz(IDX))),'alpha','cmap',cmap(iC,:))
        %plot(opt.time_vecs,nanmean(MM_valid(IDX,:)))
        xlim(xl)
        disp(mean(DELAY_valid(IDX)))
        %     subplot(2,1,2)
        %     boundedline(1:82,mean(Waveform_valid(IDX,:)),std(Waveform_valid(IDX,:)))
        %     figure
        %     imagesc(opt.time_vecs,1:numel(IDX),smoothdata(MM_valid(IDX,:),2,'gaussian',11),[-1 3])
        
    end
    boundedline(opt.time_vecs,nanmean(MM_other)*100,nanstd(MM_other)*100/sqrt(size(MM_other,1)),'alpha','cmap',[.2 .2 .2])
    legend({'smallest','middle','highest','other'})
    title(rr{iR})
    xlim([-.5 2.5])
    xlabel('time from onset')
    ylabel('% change FR')
end
%%
valid = ~isnan(DELAY) & startsWith(REGION,{'VISp','RS','MEC'}) & FIRING_RATE>1;
mm_resp = mean(MM(valid,opt.time_vecs>.05 & opt.time_vecs<0.6,:),2)-mean(MM(valid,opt.time_vecs>-.6 & opt.time_vecs<-.1),2);
mm_resp = mm_resp./FIRING_RATE(valid);
r_valid = REGION(valid);
r_valid(startsWith(r_valid,'VISp'))={'VIS'};
r_valid(startsWith(r_valid,'RS'))={'RSC'};

[u,s,d]=unique(r_valid);
figure
scatter(mm_resp,DELAY(valid),15,d)
colormap(brewermap(numel(u),'Set1'))