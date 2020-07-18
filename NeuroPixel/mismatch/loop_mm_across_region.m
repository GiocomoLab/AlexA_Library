
%matfiles = dir('Z:\giocomo\attialex\NP_DATA\mismatch\*mismatch*.mat');

%matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA_corrected/*mismatch*.mat');
matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA/*mismatch*.mat');
%matfiles = dir('/Users/attialex/NP_DATA_2/*mismatch*.mat');
%matfiles = dir('/Users/attialex/mismatch/*mismatch*.mat');
%matfiles = matfiles(~cellfun(@(x) contains(x,'tower'), {matfiles.name}));

plotfig= false;
if plotfig
    imsavedir  = '/Users/attialex/images';
    if ~isfolder(imsavedir)
        mkdir(imsavedir)
    end
end
opt = load_mismatch_opt;
opt.TimeBin = 0.02;

opt.time_bins =-1:opt.TimeBin:1;
opt.time_vecs = opt.time_bins(1:end-1)*0.5+opt.time_bins(2:end)*0.5;
opt.extract_win = [-2 3];
opt.aux_win = [-50 50];
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
REGION = {};
for iF=1:numel(matfiles)
    disp(iF)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    if ~isfield(data_out,'anatomy')
        disp('no anatomy')
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
    valid_region = startsWith(reg,{'VISp','MEC','RS'});
    if nnz(valid_region)==0
        continue
    end
    good_cells = data_out.sp.cids(data_out.sp.cgs==2 & valid_region);
    depth_this = depth(data_out.sp.cgs==2 & valid_region);
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
    %save(['/Users/attialex/temp/' sn '.mat'],'count_vec','trial_vec','trial_vec_random','good_cells','theta','mm_resp');
    if plotfig
            mm_resp = mean(count_vec(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)-mean(count_vec(:,opt.time_vecs>-0.5 & opt.time_vecs<-0.1),2);

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
[uS]=unique(SID);
frac = .2;
mm_resp = mean(MM(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)-mean(MM(:,opt.time_vecs<-.1),2);
MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
MM_All=[];
Reg_All=[];

figure('Color','white')
subplot(1,2,1)
hold on
regions = {'VISp','MEC','RS'};
for iS=1:numel(uS)
    idx_this = SID==uS(iS);
    resp_this = mm_resp(idx_this);
    mm_this = MM_ms(idx_this,:);
    [~,sid]=sort(resp_this,'descend');
    N=round(frac*nnz(idx_this));
    if startsWith(REGION{find(idx_this,1)},'MEC')
        col = [0 0 1];
        r = 2;
    elseif strcmp(REGION{find(idx_this,1)},'VISp')
        col = [1 0 0];
        r =1;
    elseif startsWith(REGION{find(idx_this,1)},'RS')
        col = [0 1 0];
        r =3;
    else
        continue
    end
    Reg_All = cat(1,Reg_All,r);
    tmp = median(mm_this(sid(1:N),:));
    MM_All = cat(1,MM_All,tmp);
    plot(opt.time_vecs,tmp,'Color',col)
end
subplot(1,2,2)
plot(opt.time_vecs,nanmedian(MM_All(Reg_All==2,:)),'b')
hold on
plot(opt.time_vecs,median(MM_All(Reg_All==1,:)),'r')
plot(opt.time_vecs,nanmedian(MM_All(Reg_All==3,:)),'k')

legend({'MEC','VISp','RS'})
%% RESp time V1
figure
frac = .5;
idx_v1= startsWith(REGION,'VISp');
mm_resp = mean(MM(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)-mean(MM(:,opt.time_vecs<-.1),2);

mm_this = MM_ms(idx_v1,opt.time_vecs>0);

depth_this = DEPTH(idx_v1);

idx_up = depth_this<median(depth_this);
idx_down = ~idx_up;
    [a,b]=max(mm_this,[],2);
    %histogram(b,'Normalization','probability')

    scatter(opt.time_vecs(b(idx_up))-opt.time_vecs(1),a(idx_up))
    hold on
        scatter(opt.time_vecs(b(~idx_up))-opt.time_vecs(1),a(~idx_up))

    %%
regions = {'VISp','MEC'};
mm_resp = mean(MM(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)-mean(MM(:,opt.time_vecs<-.1),2);
MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);

figure('Color','white')
hold on
for iR = 1:numel(regions)
    idx_this = startsWith(REGION,regions{iR});
    resp_this = mm_resp(idx_this);
    [~,sid]=sort(resp_this,'descend');
    N=round(numel(sid)*.2);
    mm_this = MM_ms(idx_this,opt.time_vecs>0);
    [a,b]=max(mm_this(sid(1:N),:),[],2);
    %histogram(b,'Normalization','probability')
    subplot(2,1,1)
    hold on
    scatter(opt.time_vecs(b)-opt.time_vecs(1),a)
    subplot(2,1,2)
    hold on
    histogram(opt.time_vecs(b)-opt.time_vecs(1),'Normalization','probability','BinEdges',[0:0.05:1])
end
subplot(2,1,1)
xlabel('time to max response')
ylabel('max response')
    legend(regions)
subplot(2,1,2)
xlabel('time to max response')
legend(regions)