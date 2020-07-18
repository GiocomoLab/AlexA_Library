matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA/mismatch/*mismatch*');

dark_folder = '/Volumes/Samsung_T5/attialex/distance_tuning_xcorr_only';

%matfiles = dir('Z:\giocomo\attialex\NP_DATA\mismatch\*mismatch*.mat');

%matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA_corrected/*mismatch*.mat');
%matfiles = dir('/Users/attialex/NP_DATA_2/*mismatch*.mat');
%matfiles = dir('/Users/attialex/mismatch/*mismatch*.mat');
matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA/mismatch/*mismatch*');
%matfiles = matfiles(~cellfun(@(x) contains(x,'tower'), {matfiles.name}));
%%
plotfig= false;
if plotfig
    imsavedir  = '/Users/attialex/images';
    if ~isfolder(imsavedir)
        mkdir(imsavedir)
    end
end
opt = load_mismatch_opt;
opt.time_bins =-2:0.02:3;
opt.time_vecs = opt.time_bins(1:end-1)*0.5+opt.time_bins(2:end)*0.5;
opt.extract_win = [-2 3];
opt.aux_win = [-50 50];
opt.TimeBin = 0.02;
opt.smoothSigma_time = 0.0; % in sec; for smoothing fr vs time

MM=[];

SID = [];
MM_R=[];
THETA_POWER = [];
cmap_fr = cbrewer('seq','BuPu',20);
SPIKE_TIMES=cell(numel(matfiles),1);
RUN_TRACES = cell(numel(matfiles),1);
CLUIDS = RUN_TRACES;
DEPTH = [];
DARK_TUNED = [];
REGION = {};
for iF=1:numel(matfiles)

    parts = strsplit(matfiles(iF).name,'_');
    session_part = strcat(parts{1},'_',parts{2});
    
    poss=dir(fullfile(dark_folder,strcat(session_part,'*')));
    
    if isempty(poss)
        disp('no dark data')
        continue
    end
    dark_data = load(fullfile(dark_folder,poss(1).name));
    
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    
    
    if ~isfield(data_out,'anatomy')
        disp('no anatomy')
        
        %continue
        valid_region = data_out.sp.cgs==2;
        depth = nan(size(data_out.sp.cgs));
        reg = zeros(size(data_out.sp.cgs));
        %depth_this = depth(data_out.sp.cgs==2 & valid_region);
        %region_this = reg(data_out.sp.cgs==2 & valid_region);
    else
        
        if ~isfield(data_out.anatomy,'depth')
            dd =data_out.anatomy.z2-data_out.anatomy.tip_distance;
            data_out.anatomy.depth = dd';
        end
        if isfield(data_out.anatomy,'parent_shifted')
            reg = data_out.anatomy.parent_shifted;
            depth = data_out.anatomy.depth_shifted';
        else
            reg = data_out.anatomy.cluster_parent;
            depth = data_out.anatomy.depth';
        end
        if iscolumn(reg)
            reg = reg';
        end
        valid_region = startsWith(reg,{'MEC'});
    end
    if nnz(valid_region)==0
        continue
    end
    
    
    mismatch_trigger = data_out.mismatch_trigger;
    if size(mismatch_trigger,1) ~=1
        mismatch_trigger=mismatch_trigger';
    end
    if nnz(mismatch_trigger==1)>nnz(mismatch_trigger==0)
        %in some old files mismatch_trigger was actually the move variable,
        %i.e. move ==0 is mismatch
        mismatch_trigger = mismatch_trigger<0.1;
    end
    
    good_cells = data_out.sp.cids(data_out.sp.cgs==2 & valid_region);
    depth_this = depth(data_out.sp.cgs==2 & valid_region);
    if isrow(depth_this)
        depth_this = depth_this';
    end
    region_this = reg(data_out.sp.cgs==2 & valid_region);
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
    [spikeTimes,~,aux]=extract_triggered_spikeTimes(data_out.sp,data_out.post(mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data_out.post' ;smooth_speed],'aux_win',opt.aux_win);
    [spikeTimesAll,~,auxAll]=extract_triggered_spikeTimes(data_out.sp,data_out.post(all_mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data_out.post' ;smooth_speed],'aux_win',opt.aux_win);
    
    SPIKE_TIMES{iF}=spikeTimesAll;
    RUN_TRACES{iF}=squeeze(auxAll);
    [spikeTimesRandom]=extract_triggered_spikeTimes(data_out.sp,data_out.post(possibles_random),'cluIDs',good_cells,'win',opt.extract_win);
    
    trial_vec =cat(1,spikeTimes{:});
    trial_vec_random = cat(1,spikeTimesRandom{:});
    CLUIDS{iF}=good_cells;
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
    
    count_vec = count_vec/n_trigs_included/opt.TimeBin;
    count_vec_random = count_vec_random/n_trigs_included_random/opt.TimeBin;
    %theta
    tB=0:opt.TimeBin:max(data_out.sp.st);
    
    frMat= zeros(numel(good_cells),numel(tB)-1);
    
    for iC=1:numel(good_cells)
        idx = data_out.sp.clu==good_cells(iC);
        frMat(iC,:)=histcounts(data_out.sp.st(idx),tB);
    end
    
    %[PxxSpikes,FSpikes]=pwelch(frMat',size(frMat,2),[],[],1/opt.TimeBin);
    
    
    %     theta_range=[4 12];
    %     theta_idx = FSpikes>theta_range(1) & FSpikes<=theta_range(2);
    %     rest_idx = ~theta_idx & FSpikes>1;
    %     thetaPower = mean(PxxSpikes(theta_idx,:));
    %     restPower = mean(PxxSpikes(rest_idx,:));
    
    %[~,theta_sort] = sort(thetaPower./(thetaPower+restPower),'descend');
    nC=size(frMat,1);
    maxLag=50;
    xcorrs=zeros(nC,2*maxLag+1);
    for iC=1:nC
        xcorrs(iC,:)=xcorr(frMat(iC,:),maxLag);
    end
    thetaPower = xcorrs(:,58);
    thetaPower_low = xcorrs(:,54);
    [~,theta_sort] = sort(thetaPower,'descend');
    theta = thetaPower-thetaPower_low;
    %     r = 1:length(thetaPower);
    %     r(theta_sort) = r;
    %     r=r/length(thetaPowerN);
    %collect data
    MM=cat(1,MM,count_vec);
    REGION = cat(2,REGION,region_this);
    MM_R = cat(1,MM_R,count_vec_random);
    THETA_POWER = cat(1,THETA_POWER,[thetaPower,thetaPower_low]);
    SID = cat(1,SID,ones(numel(good_cells),1)*iF);
    DEPTH = cat(1,DEPTH,depth_this);
    %MMR=cat(1,MMR,count_vec_random);
    [~,sn] = fileparts(matfiles(iF).name);
    mm_resp = mean(count_vec(:,105:125),2)-mean(count_vec(:,75:100),2);
    %save(['/Users/attialex/temp/' sn '.mat'],'count_vec','trial_vec','trial_vec_random','good_cells','theta','mm_resp');
    is_tuned = false(nC,1);
    peak_loc_pf = zeros(nC,1);
    MAX_RANGE_BIN = 300;
    upper_limits = squeeze(dark_data.acg_upper_limits(:,:,end));
    for iN = 1:nC
        
        
        a=strfind(diff(dark_data.ACG(iN,:))>0,[0 1]); %find elbow point/location of first minimum
        if numel(a)<1 %no elbow point detected, skip
            continue
        end
        
        if a(1)>MAX_RANGE_BIN % elbow point beyond 500cm, probably noise, skip
            continue
        end
        [pks,loc]=findpeaks(dark_data.ACG(iN,1:MAX_RANGE_BIN),'SortStr','descend','NPeaks',1,'MinPeakProminence',.2);
        if upper_limits(iN,loc)<pks
            peak_loc_pf(iN)=loc;
        end
        
        sub = dark_data.ACG(iN,a(1):MAX_RANGE_BIN); %only search within between elbow point and MAX_RANGE_BIN
        [ma,mi]=max(sub);
        mi = mi+a(1)-1;
        is_tuned(iN) = upper_limits(iN,mi)<ma;
%         mod_depth(iN)=ma-ACG(iN,a(1));
%         peak_loc(iN)=mi;
        
        
    end
    DARK_TUNED = cat(1,DARK_TUNED,is_tuned);
    
    if plotfig
        imsavedir_this = fullfile(imsavedir,sn);
        if ~isfolder(imsavedir_this)
            mkdir(imsavedir_this)
        end
        
        idx_pre = opt.time_vecs<0 & opt.time_vecs>-.5;
        idx_post = opt.time_vecs>0.1 & opt.time_vecs<0.6;
        mm_response = mean(count_vec(:,idx_post),2)-mean(count_vec(:,idx_pre),2);
        %frMat = calcTrialFRMat(good_cells,1:max(data_out.trial)-1,data_out,load_default_opt);
        [cm,frMat,~]=trialCorrMat(good_cells,1:max(data_out.trial)-1,data_out,load_default_opt);
        figure
        imagesc(squeeze(nanmean(cm)),[0 .7])
        axis image
        saveas(gcf,fullfile(imsavedir_this,'sim_map.png'));
        close(gcf)
        [~,sidx]=sort(mm_response,'descend');
        figure
        for ii=1:10
            subplot(2,2,[1 3])
            imagesc(1:2:400,1:max(data_out.trial)-1,squeeze(frMat(sidx(ii),:,:)))
            hold on
            plot(data_out.posx(mm_trigs),data_out.trial(mm_trigs),'ro')
            colormap(cmap_fr)
            hold on
            subplot(2,2,2)
            idx = trial_vec(:,2)==good_cells(sidx(ii));
            scatter(trial_vec(idx,1),trial_vec(idx,3),15,'.')
            subplot(2,2,4)
            plot(opt.time_vecs,count_vec(sidx(ii),:));
            saveas(gcf,fullfile(imsavedir_this,sprintf('%d.png',good_cells(sidx(ii)))));
            clf
        end
        close(gcf)
    end
end
%%
figure
MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
imagesc((MM_ms),[-20 20])
cmap = cbrewer('div','RdBu',20);
cmap=flipud(cmap);
colormap(cmap);
%%
params=struct();
params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
params.xLim=[-2 3];

figure
plotAVGSEM(MM(logical(DARK_TUNED),:)',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[0 0 1])
plotAVGSEM(MM(~logical(DARK_TUNED),:)',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[1 0 0])

legend({'Tuned','no tuned'})

%%
params=struct();
params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
params.xLim=[-2 3];
uS= unique(SID);
for iS=1:numel(uS)
figure
plotAVGSEM(MM(logical(DARK_TUNED) & SID == uS(iS),:)',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[0 0 1])
plotAVGSEM(MM(~logical(DARK_TUNED) & SID == uS(iS),:)',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[1 0 0])

legend({'Tuned','no tuned'})
end

