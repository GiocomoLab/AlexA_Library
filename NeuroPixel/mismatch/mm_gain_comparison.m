
%matfiles = dir('Z:\giocomo\attialex\NP_DATA\mismatch\*mismatch*.mat');

%matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA_2/*mismatch*.mat');
matfiles = dir('/Users/attialex/NP_DATA_2/*mismatch*.mat');
matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA_corrected/*mismatch*.mat');
matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA/mismatch/*mismatch*.mat');

%matfiles = dir('/Users/attialex/mismatch/*mismatch*.mat');
matfiles = matfiles(cellfun(@(x) contains(x,'tower'), {matfiles.name}));

plotfig= false;
if plotfig
    imsavedir  = '/Users/attialex/images_stab';
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
MM_half1=[];
MM_half2 = [];
SID = [];
MM_R=[];
THETA_POWER = [];
CorrMat={};
cmap_fr = cbrewer('seq','BuPu',20);
SPIKE_TIMES=cell(numel(matfiles),1);
RUN_TRACES = cell(numel(matfiles),1);
CLUIDS = RUN_TRACES;
DEPTH = [];
GAIN = {};
MM_Trials = {};
for iF=1:numel(matfiles)
    disp(iF)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    if ~isfield(data_out,'anatomy')
        fprintf('no anatomy for %s',matfiles(iF).name);
        continue
        %valid_region = data_out.sp.cgs==2;
    else
        if isfield(data_out.anatomy,'parent_shifted')
            reg = data_out.anatomy.parent_shifted;
            %depth = data_out.anatomy.depth_shifted';
        else
            reg = data_out.anatomy.cluster_parent;
            %depth = data_out.anatomy.depth';
        end
        if iscolumn(reg)
            reg = reg';
        end
        valid_region = startsWith(reg,'MEC');
    end
    if nnz(valid_region)==0
        continue
    end
    
    good_cells = data_out.sp.cids(data_out.sp.cgs==2 & valid_region);
    mismatch_trigger = data_out.mismatch_trigger;
    if size(mismatch_trigger,1) ~=1
        mismatch_trigger=mismatch_trigger';
    end
    if nnz(mismatch_trigger==1)>nnz(mismatch_trigger==0)
        %in most files, mm trigger is actually move variable,
        %i.e. move ==0 is mismatch
        mismatch_trigger = mismatch_trigger<0.1;
    end
    
    all_mm_trigs=strfind(mismatch_trigger,[0 0 1 1])+2;
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
    p = [all_mm_trigs;ismember(all_mm_trigs,possibles);data_out.posx(all_mm_trigs)';data_out.trial(all_mm_trigs)'];
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
    [spikeTimes,~,aux,mm_trigs_idx]=extract_triggered_spikeTimes(data_out.sp,data_out.post(mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data_out.post' ;smooth_speed],'aux_win',opt.aux_win);
    mm_trigs = mm_trigs(mm_trigs_idx);
    [spikeTimesAll,~,auxAll]=extract_triggered_spikeTimes(data_out.sp,data_out.post(all_mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data_out.post' ;smooth_speed],'aux_win',opt.aux_win);
    first_half = data_out.posx(mm_trigs)<=200;
    second_half = data_out.posx(mm_trigs)>200;
    SPIKE_TIMES{iF}=spikeTimesAll;
    RUN_TRACES{iF}=squeeze(auxAll);
    
    trial_vec =cat(1,spikeTimes{:});
    trial_vec_1 =cat(1,spikeTimes{first_half});
    
    trial_vec_2 =cat(1,spikeTimes{second_half});
    
    CLUIDS{iF}=good_cells;
    count_vec = zeros(numel(good_cells),numel(opt.time_bins)-1);
    count_vec_1 = count_vec;
    count_vec_2 = count_vec;
    for iC=1:numel(good_cells)
        idx = trial_vec(:,2)==good_cells(iC);
        [spike_count]=histcounts(trial_vec(idx,1),opt.time_bins);
        count_vec(iC,:)=spike_count;
        
        idx = trial_vec_1(:,2)==good_cells(iC);
        [spike_count]=histcounts(trial_vec_1(idx,1),opt.time_bins);
        count_vec_1(iC,:)=spike_count;
        
        idx = trial_vec_2(:,2)==good_cells(iC);
        [spike_count]=histcounts(trial_vec_2(idx,1),opt.time_bins);
        count_vec_2(iC,:)=spike_count;
        
        
    end
    n_trigs_included = numel(unique(trial_vec(:,3)));
    n_trigs_included_1 = numel(unique(trial_vec_1(:,3)));
    n_trigs_included_2 = numel(unique(trial_vec_2(:,3)));
    
    count_vec = count_vec/n_trigs_included/opt.TimeBin;
    count_vec_1 = count_vec_1/n_trigs_included_1/opt.TimeBin;
    count_vec_2 = count_vec_2/n_trigs_included_2/opt.TimeBin;
    
    %theta
    tB=0:opt.TimeBin:max(data_out.sp.st);
    
    %% trial corr mat
    [corrMat_this,frMat_this,shiftMat_this] = trialCorrMat(good_cells,1:(max(data_out.trial)-1),data_out,load_default_opt);
    
    trial_ends = find(diff(data_out.posx)<-100);
    trial_gain = zeros(size(trial_ends));
    stop = 0;
    for iT = 1:numel(trial_ends)
        start = stop+1;
        stop = trial_ends(iT);
        trial_gain(iT)=mode(data_out.mismatch_trigger(start:stop));
        
    end
    
    
    %%
    GAIN{end+1}=trial_gain;
    MM_Trials{end+1}=p;
    MM=cat(1,MM,count_vec);
    MM_half1 = cat(1,MM_half1,count_vec_1);
    MM_half2 = cat(1,MM_half2,count_vec_2);
    CorrMat{end+1} = corrMat_this;
    SID = cat(1,SID,ones(numel(good_cells),1)*iF);
    %DEPTH = cat(1,DEPTH,depth_this);
    %MMR=cat(1,MMR,count_vec_random);
    [~,sn] = fileparts(matfiles(iF).name);
    mm_resp = mean(count_vec(:,105:125),2)-mean(count_vec(:,75:100),2);
    %save(['/Users/attialex/temp/' sn '.mat'],'count_vec','trial_vec','trial_vec_random','good_cells','theta','mm_resp');
    if plotfig
        imsavedir_this = fullfile(imsavedir,sn);
        if ~isfolder(imsavedir_this)
            mkdir(imsavedir_this)
        end
        
        idx_pre = opt.time_vecs<0 & opt.time_vecs>-.5;
        idx_post = opt.time_vecs>0.1 & opt.time_vecs<0.6;
        mm_response = mean(count_vec(:,idx_post),2)-mean(count_vec(:,idx_pre),2);
        %frMat = calcTrialFRMat(good_cells,1:max(data_out.trial)-1,data_out,load_default_opt);
        %[cm,frMat,~]=trialCorrMat(good_cells,1:max(data_out.trial)-1,data_out,load_default_opt);
        figure
        subplot(1,10,1)
        imagesc(1,1:numel(trial_gain),trial_gain)
        hold on
        plot(ones(size(all_mm_trigs)),data_out.trial(all_mm_trigs),'ro')
        plot(ones(size(mm_trigs)),data_out.trial(mm_trigs),'g.')
        subplot(1,10,2:10)
        imagesc(squeeze(nanmean(corrMat_this)))
        axis image
        saveas(gcf,fullfile(imsavedir_this,'sim_map.png'));
        close(gcf)
        %[~,sidx]=sort(mm_response,'descend');
        tmp = nanmean(nanmean(corrMat_this,2),3);
        [~,sidx]=sort(tmp,'descend');
        figure
        for ii=1:10
            subplot(2,2,[1 3])
            imagesc(1:2:400,1:max(data_out.trial)-1,squeeze(frMat_this(sidx(ii),:,:)))
            hold on
            plot(data_out.posx(mm_trigs),data_out.trial(mm_trigs),'ro')
            colormap(cmap_fr)
            hold on
            subplot(2,2,2)
            idx = trial_vec(:,2)==good_cells(sidx(ii));
            scatter(trial_vec(idx,1),trial_vec(idx,3),15,'.')
            subplot(2,2,4)
            plot(opt.time_vecs,count_vec(sidx(ii),:));
            hold on
            plot(opt.time_vecs,count_vec_1(sidx(ii),:),'r--')
            plot(opt.time_vecs,count_vec_2(sidx(ii),:),'g--')
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
%% average all
params=struct();
params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
params.xLim=[-2 3];
figure
plotAVGSEM(MM',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0)
%% first/second
params=struct();
params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
params.xLim=[-2 3];
figure('Color','white')
subplot(1,3,1)
plotAVGSEM(MM_half1',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0)
plotAVGSEM(MM_half2',gca,'parameters',params,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[ 1 0 0])

bl = opt.time_vecs>-.25 & opt.time_vecs<0;
t1 = opt.time_vecs>.1 & opt.time_vecs<.5;
resp_1 = mean(MM_half1(:,t1),2)-mean(MM_half1(:,bl),2);
resp_2 = mean(MM_half2(:,t1),2)-mean(MM_half2(:,bl),2);
subplot(1,3,2)
STAB = [];
for iS=1:numel(CorrMat)
    
    ons = strfind(GAIN{iS}'==0.5,[0 1])+2;
    if isempty(ons)
        ons=size(CorrMat{iS},2);
        
    end
    stab = squeeze(nanmean(nanmean(CorrMat{iS}(:,1:ons,1:ons),2),3));
    idx = SID==iS;
    STAB = cat(1,STAB,stab);
end
scatter(resp_1,resp_2,15,STAB,'filled')
set(gca,'CLim',[0 .7])
axis image

xlabel('response first half')
ylabel('response second half')

tidx = 80:150;
n_samples_per_cond = numel(tidx);

mm_resp = MM_half1(:,tidx)-mean(MM_half1(:,75:100),2);
pb_resp = MM_half2(:,tidx)-mean(MM_half2(:,75:100),2);

smooth_method = 'gaussian';
smooth_window = 9;
mm_resp = smoothdata(mm_resp,2,smooth_method,smooth_window);
pb_resp = smoothdata(pb_resp,2,smooth_method,smooth_window);
mm_resp = zscore(mm_resp,[],2);
pb_resp = zscore(pb_resp,[],2);
r_idx = randperm(size(pb_resp,1));
%pb_resp=pb_resp(r_idx,:);


data_matrix = cat(1,mm_resp',pb_resp');
color = cat(1,hsv(n_samples_per_cond),hsv(n_samples_per_cond));

[coeff,score,latent,tsquared,explained,mu] = pca(data_matrix);
score(n_samples_per_cond,:)=nan;
subplot(1,3,3)
scatter3(score(:,1),score(:,2),score(:,3),45,color,'.')
hold on
patch([score(:,1); nan],[score(:,2); nan],[score(:,3); nan],[opt.time_vecs(tidx) opt.time_vecs(tidx) nan],'FaceColor','none','EdgeColor','interp','LineWidth',2)
colormap(parula)
%% average corrmap around trial 10 with/without stability
AllC=[];
for iS=1:numel(CorrMat)
    this = CorrMat{iS}(:,1:20,1:20);
    stab_pre = nanmean(nanmean(this(:,6:10,6:10),2),3);
    stab_post = nanmean(nanmean(this(:,11:15,11:15),2),3);
    idx = stab_pre>.2 & stab_post>.2;
    AllC = cat(3,AllC,squeeze(nanmean(this(idx,:,:))));
end
figure
imagesc(squeeze(nanmean(AllC,3)))
%%
[uS]=unique(SID);
AllC=[];
bl = opt.time_vecs>-.25 & opt.time_vecs<0;
t1 = opt.time_vecs>.1 & opt.time_vecs<.5;
resp = mean(MM(:,t1),2)-mean(MM(:,bl),2);
MMResp = [];
STAB=[];
for iS=1:numel(uS)
    if ~ismember(0.5,GAIN{iS})
        continue
    end
    ons = strfind(GAIN{iS}'==0.5,[0 1])+2;
    
    pre=ons-6:(ons-1);
    post = ons:(ons+3);
    snip = CorrMat{iS}(:,ons-6:ons+3,ons-6:ons+3);
    stab = [squeeze(nanmean(nanmean(snip(:,1:5,1:5),2),3)) squeeze(nanmean(nanmean(snip(:,1:5,6:7),2),3))];
    idx = SID==iS;
    STAB = cat(1,STAB,stab);
    MMResp = cat(1,MMResp,resp(idx));
    %AllC=cat(3,AllC,squeeze(nanmean(CorrMat{iS}(:,ons-6:ons+3,ons-6:ons+3))));
end
%figure
%imagesc(nanmean(AllC,3))
figure('Color','White')
subplot(1,2,1)
scatter(STAB(:,1),STAB(:,2),15,MMResp,'filled')
cmap = cbrewer('div','RdBu',20);
cmap=flipud(cmap);
colormap(cmap)
set(gca,'CLim',[-10 10])
axis image
lims = [-.1 1];
xlim(lims)
ylim(lims)
grid on
xlabel('stab pre')
ylabel('stab post')
subplot(1,2,2)
scatter(STAB(:,1),MMResp,15,'filled')
xlabel('stab pre')
ylabel('MM resp')

%%
%%
[uS]=unique(SID);
AllC=[];
bl = opt.time_vecs>-.25 & opt.time_vecs<0;
t1 = opt.time_vecs>.1 & opt.time_vecs<.5;
resp = mean(MM(:,t1),2)-mean(MM(:,bl),2);
%tmp = MM./mean(MM(:,opt.time_vecs<0),2);
%resp = tmp;
resp(isinf(tmp))=nan;
%resp = nanmean(resp(:,t1),2)-nanmean(resp(:,bl),2);

MMResp = [];
STAB=[];
for iS=1:numel(uS)
    
    ons = strfind(GAIN{iS}'==0.5,[0 1])+2;
    if isempty(ons)
        ons=size(CorrMat{iS},2);
        
    end
    stab = squeeze(nanmean(nanmean(CorrMat{iS}(:,1:ons,1:ons),2),3));
    idx = SID==uS(iS);
    STAB = cat(1,STAB,stab);
    MMResp = cat(1,MMResp,resp(idx));
    %AllC=cat(3,AllC,squeeze(nanmean(CorrMat{iS}(:,ons-6:ons+3,ons-6:ons+3))));
end
%figure
%imagesc(nanmean(AllC,3))
figure('Color','White')
subplot(1,2,1)
scatter(STAB,MMResp,15,'filled')
xlabel('stab')
ylabel('MM resp')

bin_width = 0.1;
bin_starts = 0:0.05:1;
mm=nan(numel(bin_starts),2);
bin_centers = bin_starts+bin_width/2;
for iB=1:numel(bin_starts)
    idx = STAB>bin_starts(iB) & STAB<=bin_starts(iB)+bin_width;
    mm(iB,:)=[nanmean(MMResp(idx)),nanstd(MMResp(idx))/sqrt(nnz(idx))];
end
subplot(1,2,2)
hold on
errorbar(bin_centers,mm(:,1),mm(:,2))
xlabel('stab')
ylabel('MM resp')