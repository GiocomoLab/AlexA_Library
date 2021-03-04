%addpath(genpath('F:\code\cortexlab_spikes\analysis\'))

matfiles = dir('Z:\giocomo\attialex\NP_DATA\*mismatch*.mat');
wave_path = 'I:\mean_waveforms\mean_waveforms';
matfiles = dir('I:\mismatch_mec\np*mismatch*.mat');
%matfiles = dir('/Users/attialex/NP_DATA_2/*mismatch*.mat');
matfiles = dir('/Volumes/T7/attialex/mismatch_mec/np*mismatch*.mat');
%matfiles = matfiles(~cellfun(@(x) contains(x,'tower'), {matfiles.name}));
histo_table = process_histology();
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
WAVEFORM = [];
SID = [];
MM_R=[];
THETA_POWER = [];
FIRING_RATE = [];
%cmap_fr = cbrewer('seq','BuPu',20);
cmap_fr = summer(20);
SPIKE_TIMES=cell(numel(matfiles),1);
RUN_TRACES = cell(numel(matfiles),1);
CLUIDS = RUN_TRACES;
DEPTH = [];
REGION = {};
FR=[];
POS3D=[];
SPIKE_CHANNEL = [];
templateDuration = [];
templateWaveform = [];
MM_THETA={};
for iF=1:numel(matfiles)
    disp(iF)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
        templatePositionsAmplitudes(data_out.sp.temps, data_out.sp.winv, data_out.sp.ycoords, data_out.sp.spikeTemplates, data_out.sp.tempScalingAmps);
    nC=numel(data_out.sp.cgs);
    
    spike_channel = nan(nC,1);
    for ii=1:nC
        spike_channel(ii)=nanmedian(spikeDepths(data_out.sp.clu==data_out.sp.cids(ii)));
    end
    if ~isfield(data_out,'anatomy')
        disp(matfiles(iF).name)
        disp('no anatomy')
        reg = repmat({'UNK'},1,nC);
        depth = nan(nC,1); % go via channel number
        %continue
        valid_region = true(size(reg));
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
        valid_region = startsWith(reg,{'MEC','ECT'});
    end
    if nnz(valid_region)==0
        continue
    end
    
    parts = split(matfiles(iF).name,'_');
    animal = parts{1};
    date = parts{2};
    
    session_num=find(startsWith(histo_table.animal,animal) & startsWith(histo_table.date,date));
    if numel(session_num) ==1
    pos3D = histo_table.origin(session_num,:)+spike_channel.*histo_table.unit_vector(session_num,:);
    else
        pos3D = nan(nC,3);
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
    region_this = reg(data_out.sp.cgs==2 & valid_region);
    pos3d_this = pos3D(data_out.sp.cgs==2 & valid_region,:);
    all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
    true_speed = data_out.true_speed;
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
    if all(size(data_out.post)==size(speed))
        speed = speed';
    end
    %[~,sp] = calcSpeed(data_out.posx,load_default_opt);
    [spikeTimes,~,aux,~,count_vec]=extract_triggered_spikeTimes(data_out.sp,data_out.post(mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data_out.post' ;smooth_speed],'aux_win',opt.aux_win);
    [spikeTimesAll,~,auxAll,~,count_vec_all]=extract_triggered_spikeTimes(data_out.sp,data_out.post(all_mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data_out.post' ;smooth_speed],'aux_win',opt.aux_win);
    
    SPIKE_TIMES{iF}=spikeTimesAll;
    RUN_TRACES{iF}=squeeze(auxAll);
    [spikeTimesRandom,~,~,~,count_vec_random]=extract_triggered_spikeTimes(data_out.sp,data_out.post(possibles_random),'cluIDs',good_cells,'win',opt.extract_win);
    
    trial_vec =cat(1,spikeTimes{:});
    trial_vec_random = cat(1,spikeTimesRandom{:});
    CLUIDS{iF}=good_cells;
    
    %theta
    tB=0:opt.TimeBin:max(data_out.sp.st);
    
    frMat= zeros(numel(good_cells),numel(tB)-1);
    firing_rate = zeros(numel(good_cells),1);
    for iC=1:numel(good_cells)
        idx = data_out.sp.clu==good_cells(iC);
        firing_rate(iC)=nnz(idx);
        frMat(iC,:)=histcounts(data_out.sp.st(idx),tB);
    end
    firing_rate = firing_rate/(data_out.post(end)-data_out.post(1));
    %[PxxSpikes,FSpikes]=pwelch(frMat',size(frMat,2),[],[],1/opt.TimeBin);
    
    
    %     theta_range=[4 12];
    %     theta_idx = FSpikes>theta_range(1) & FSpikes<=theta_range(2);
    %     rest_idx = ~theta_idx & FSpikes>1;
    %     thetaPower = mean(PxxSpikes(theta_idx,:));
    %     restPower = mean(PxxSpikes(rest_idx,:));
    trigs = [];
    [p,s,t]=pspectrum(frMat(iC,:),1/opt.TimeBin,'spectrogram','FrequencyLimits',[0 10],'FrequencyResolution',[1]);
    for iT=1:numel(mm_trigs)
        [~,idx]=min(abs(data_out.post(mm_trigs(iT))-t));
        trigs(end+1)=idx;
    end
    power_estimate = zeros(numel(good_cells),81);
    power_estimate_low = power_estimate;
    mean_specgram = [];
    for iC=1:numel(good_cells)
        [p,s,t]=pspectrum(frMat(iC,:),1/opt.TimeBin,'spectrogram','FrequencyLimits',[0 10],'FrequencyResolution',[1]);
        snps = extract_snps(p,trigs,'win',[-40 40]);
        est = mean(mean(snps(s>6 & s<8,:,:),1),3);
        est_l= mean(mean(snps(s< 6 & s>=1,:,:),1),3);
        power_estimate(iC,:)=est;
        power_estimate_low(iC,:)=est_l;
        mean_specgram = cat(3,mean_specgram,mean(snps,3));
    end
    dt = mean(diff(t));
    time_vec = linspace(-40*dt,40*dt,81);
    P.power_estimate = power_estimate;
    P.power_estimate_low = power_estimate_low;
    P.time_vec = time_vec;
    P.mean_specgram = mean_specgram;
    MM_THETA=cat(1,MM_THETA,{P});
    %[~,theta_sort] = sort(thetaPower./(thetaPower+restPower),'descend');
    nC=size(frMat,1);
    maxLag=50;
    xcorrs=zeros(nC,2*maxLag+1);
    for iC=1:nC
        xcorrs(iC,:)=xcorr(frMat(iC,:),maxLag,'coeff');
    end
    thetaPower = xcorrs(:,58);
    thetaPower_low = xcorrs(:,54);
    [~,theta_sort] = sort(thetaPower,'descend');
    theta = thetaPower-thetaPower_low;
    %     r = 1:length(thetaPower);
    %     r(theta_sort) = r;
    %     r=r/length(thetaPowerN);
    %collect data
    if isfile(fullfile(wave_path,matfiles(iF).name))
        wave_data = load(fullfile(wave_path,matfiles(iF).name));
        %p = ismember(wave_data.clusters_good,good_cells);
        %temp = data_out.sp.cgs==2;
        tmp_reg = reg(data_out.sp.cgs==2);
        %tmp_idx = startsWith(tmp_reg,'VISp');
        tmp_idx = true(size(tmp_reg));
        waveform = wave_data.mean_waveforms(tmp_idx,:);
        
    else
        waveform = nan(numel(good_cells),82);
    end
    %tempPeakWF
    %cluID_this = unique(data_out.sp.clu);
    tmpWF = nan(numel(good_cells),82);
    tmpDur = nan(numel(good_cells),1);
    for iC=1:numel(good_cells)
        temp_id=unique(data_out.sp.spikeTemplates(data_out.sp.clu==good_cells(iC)));
        if ~isempty(temp_id)
            temp_id = temp_id(1)+1;
            tmpWF(iC,:)=tempPeakWF(temp_id,:);
            tmpDur(iC)=tempDur(temp_id);
        end
    end
    
    templateWaveform = cat(1,templateWaveform,tmpWF);
    templateDuration = cat(1,templateDuration,tmpDur);
    WAVEFORM = cat(1,WAVEFORM,waveform);
    MM=cat(1,MM,count_vec);
    REGION = cat(2,REGION,region_this);
    MM_R = cat(1,MM_R,count_vec_random);
    FR= cat(1,FR,firing_rate);
    THETA_POWER = cat(1,THETA_POWER,[thetaPower,thetaPower_low]);
    SID = cat(1,SID,ones(numel(good_cells),1)*iF);
    DEPTH = cat(1,DEPTH,depth_this);
    POS3D = cat(1,POS3D,pos3d_this);
    %MMR=cat(1,MMR,count_vec_random);
    [~,sn] = fileparts(matfiles(iF).name);
    mm_resp = mean(count_vec(:,105:125),2)-mean(count_vec(:,75:100),2);
    freq_idx = s>6 & s<8;
    theta_power = mean(mean(P.mean_specgram(freq_idx,:,:),1),2);
    theta_powerN = squeeze(theta_power./mean(mean(P.mean_specgram(~freq_idx,:,:),1),2));
    theta = theta_power;
    theta_norm = theta_powerN;
    save(['/Users/attialex/temp/' sn '.mat'],'count_vec','trial_vec','trial_vec_random','good_cells','theta','mm_resp','theta_norm','firing_rate');
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
[~,sid]=sort(abs(templateDuration));
MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
imagesc((MM_ms(sid,:)),[-20 20])
%cmap = cbrewer('div','RdBu',20);
%cmap=flipud(cmap);
%colormap(cmap);
%%
params=struct();
params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
params.xLim=[-2 3];
figure
plotAVGSEM(MM',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0)
plotAVGSEM(MM_R',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[.5 .5 .5])

%%
RANK=[];
[uS]=unique(SID);
for iS=1:numel(uS)
    idx = SID==uS(iS);
    resp = diff(THETA_POWER(idx,:),[],2);
    ranking=1:numel(resp);
    [~,sid]=sort(resp);
    ranking(sid)=ranking/numel(resp);
    RANK = cat(1,RANK,ranking');
    
end

%%
[~,sidx]=sort(diff(THETA_POWER,[],2));
MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);

chunksize=203;
nChunks = floor(size(MM,1)/chunksize);
%cmap = cbrewer('div','RdBu');
cmap = brewermap(21,'RdBu');
cmap=flipud(cmap);
for iC=[1 nChunks]%nChunks
    figure
    sub_idx=(iC-1)*chunksize+(1:chunksize);
    IDX = sidx(sub_idx);
    %IDX=RANK<.1
    params=struct();
    params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
    params.xLim=[-2 3];
    subplot(2,1,1)
    plotAVGSEM(MM(IDX,:)',gca,'parameters',params,'ms',false,'baseline',opt.time_bins>=-.5 & opt.time_bins<0)
    plotAVGSEM(MM_R(IDX,:)',gca,'parameters',params,'ms',false,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[.5 .5 .5])
    
    
    
    subplot(2,1,2)
    imagesc(opt.time_bins,1:nnz(IDX),MM_ms(IDX,:),opt.aux_win)
    
    colormap(cmap);
    
end
% subplot(1,2,2)
% MM_ms = MM_R-mean(MM_R(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
% imagesc(MM_ms(IDX,:),[-1 1])
% cmap = cbrewer('div','RdBu',20);
% cmap=flipud(cmap);
% colormap(cmap);
%%
MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);

chunksize=203;
nChunks = floor(size(MM,1)/chunksize);
cmap = brewermap(21,'RdBu');
cmap=flipud(cmap);

for iC=[0 .9;.1 1]%nChunks
    figure
    
    IDX = RANK>= iC(1) & RANK<=iC(2);
    %IDX=RANK<.1
    params=struct();
    params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
    params.xLim=[-2 3];
    subplot(2,1,1)
    plotAVGSEM(MM(IDX,:)',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0)
    plotAVGSEM(MM_R(IDX,:)',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[.5 .5 .5])
    
    
    
    subplot(2,1,2)
    imagesc(opt.time_bins,1:nnz(IDX),MM_ms(IDX,:),[-50 50])
    
    colormap(cmap);
    
end


%% double peak V1
%cmap = cbrewer('div','RdBu',20);
%cmap=flipud(cmap);
cmap = summer(20);
t2 = opt.time_vecs>.5 & opt.time_vecs<1;
bl = opt.time_vecs>-.25 & opt.time_vecs<0;
t1 = opt.time_vecs>.1 & opt.time_vecs<.5;
resp = mean(MM(:,t2),2)-mean(MM(:,t1),2);
[~,sidx]=sort(resp);
figure
imagesc(opt.time_vecs,1:numel(sidx),MM(sidx,:)-mean(MM(sidx,bl),2),[-10 10])
colormap(cmap)

MM_ms = MM-mean(MM(:,bl),2);

chunksize=250;
nChunks = floor(size(MM,1)/chunksize);

for iC=[1:nChunks]%nChunks
    figure
    sub_idx=(iC-1)*chunksize+(1:chunksize);
    IDX = sidx(sub_idx);
    %IDX=RANK<.1
    params=struct();
    params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
    params.xLim=[-2 3];
    subplot(3,1,1)
    plotAVGSEM(MM(IDX,:)',gca,'parameters',params,'ms',false,'baseline',opt.time_bins>=-.5 & opt.time_bins<0)
    plotAVGSEM(MM_R(IDX,:)',gca,'parameters',params,'ms',false,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[.5 .5 .5])
    
    
    
    subplot(3,1,2)
    imagesc(opt.time_bins,1:nnz(IDX),MM_ms(IDX,:),opt.aux_win)
    
    colormap(cmap);
    subplot(3,1,3)
    plot(nanmean(WAVEFORM(IDX,:)))
end
%%
chunksize=100;
nChunks = floor(size(MM,1)/chunksize);
figure
hold on
for iC=1:nChunks
    sub_idx=(iC-1)*chunksize+(1:chunksize);
    IDX = sidx(sub_idx);
    plot(nanmean(WAVEFORM(IDX,:)))
end
%% todo even/odd, waveforms, depth vs sorting second peak
