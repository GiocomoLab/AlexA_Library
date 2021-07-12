
%matfiles = dir('Z:\giocomo\attialex\NP_DATA\mismatch\*mismatch*.mat');

%matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA_corrected/*mismatch*.mat');
%matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA/*mismatch*.mat');
%matfiles = dir('/Users/attialex/NP_DATA_2/*mismatch*.mat');
%matfiles = dir('/Users/attialex/mismatch/*mismatch*.mat');
%matfiles = matfiles(~cellfun(@(x) contains(x,'tower'), {matfiles.name}));
matfiles = dir('F:\Alex\new_2\*_MM_*');
matfiles = cat(1,matfiles,dir('F:\Alex\new_2\*_mismatchSquare_*'));
%matfiles = cat(1,matfiles,dir('Z:\giocomo\attialex\NP_DATA_corrected\*mismatch*.mat'));
savepath = 'F:\Alex\mm_out';
if ~isfolder(savepath)
    mkdir(savepath)
end
plotfig= false;
if plotfig
    imsavedir  = '/Users/attialex/images';
    if ~isfolder(imsavedir)
        mkdir(imsavedir)
    end
end
opt = load_mismatch_opt;
opt.TimeBin = 0.02;

opt.extract_win = [-2 3];
opt.aux_win = [-50 50];
opt.smoothSigma_time = 0.5; % in sec; for smoothing fr vs time

MM=[];

SID = [];
MM_R=[];
THETA_POWER = [];
CORR_M=[];
CORR_PB=[];
DELAYS=[];

cmap = brewermap(20,'BuPu');
SPIKE_TIMES=cell(numel(matfiles),1);
RUN_TRACES = cell(numel(matfiles),1);
CLUIDS = RUN_TRACES;
DEPTH = [];
REGION = {};
FR=[];
pb_file_list={};
%%
% for iF=1:numel(matfiles)
%     data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
%     if isfield(data_out,'anatomy')
%         unique(data_out.anatomy.cluster_parent)'
%          if nnz(startsWith(data_out.anatomy.cluster_parent,'ENTm'))>0
%             keyboard
%         end
%     end
% end
%%
for iF=1:numel(matfiles)
    disp(iF)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    if ~isfield(data_out,'anatomy')
        disp('no anatomy')
        continue
    end
    
    if istable(data_out.anatomy)
        good_cells = data_out.sp.cids(data_out.sp.cgs==2 & ismember(data_out.sp.cids,data_out.anatomy.cluster_id));
        clu_histo_idx = ismember(data_out.anatomy.cluster_id,good_cells);
        CLUIDS{iF}=[good_cells,data_out.anatomy.cluster_id(clu_histo_idx)];
        region_this = data_out.anatomy.cluster_parent(clu_histo_idx);
        subregion_this = data_out.anatomy.cluster_region(clu_histo_idx);
%         if nnz(startsWith(data_out.anatomy.cluster_parent,'ENTm'))>0
%             keyboard
%         end
    else
        good_cells = data_out.sp.cids(data_out.sp.cgs==2);
        CLUIDS{iF}=[good_cells,good_cells];
        if isfield(data_out.anatomy,'parent_shifted')
            reg = data_out.anatomy.parent_shifted;
            subreg = data_out.anatomy.cluster_shifted;
            %depth = data_out.anatomy.depth_shifted';
        else
            reg = data_out.anatomy.cluster_parent;
            if isfield(data_out.anatomy,'cluster_region')
                subreg = data_out.anatomy.cluster_region;
            else
                subreg = reg; %for old MEC data
            end
            %depth = data_out.anatomy.depth';
        end
        reg=reg(data_out.sp.cgs==2);
        subregion_this = subreg(data_out.sp.cgs==2);
        if ~iscolumn(reg)
            reg = reg';
            subregion_this=subregion_this';
        end
        region_this = reg;
    end
    
    
    firing_rate = nan(numel(good_cells),1);
    for iC=1:numel(good_cells)
        firing_rate(iC)=nnz(data_out.sp.clu==good_cells(iC));
    end
    firing_rate = firing_rate/(data_out.post(end)-data_out.post(1));
    
    
    [spike_times,count_vec,aux_mat]=extractMM(data_out,good_cells,opt);
    
    
        
    
    
    %     tB=0:opt.TimeBin:max(data_out.sp.st);
    %
    %     frMat= zeros(numel(good_cells),numel(tB)-1);
    %
    %     for iC=1:numel(good_cells)
    %         idx = data_out.sp.clu==good_cells(iC);
    %         frMat(iC,:)=histcounts(data_out.sp.st(idx),tB);
    %     end
    %
    %[PxxSpikes,FSpikes]=pwelch(frMat',size(frMat,2),[],[],1/opt.TimeBin);
    if isfield(data_out,'true_speed')
        true_speed = data_out.true_speed;
    else
        true_speed = data_out.vr_data_resampled.velM;
    end
    velM_s = gauss_smoothing(true_speed,opt.smoothSigma_time/opt.TimeBin)/opt.TimeBin;
    if isrow(velM_s)
        velM_s=velM_s';
    end
    
    [fr,spikeCount] = calcFRVsTime(good_cells,data_out,opt);
    if isrow(velM_s)
        velM_s=velM_s';
    end
    corr_m=corr(fr',velM_s);
    
    pb_file=find(strcmp(cat(1,files{:,2}),matfiles(iF).name));
    if ~isempty(strfind(matfiles(iF).name,'mismatch'))
        pb_name = strrep(matfiles(iF).name,'mismatch','playback');
    else
        pb_name = '';
    end
    corr_mpb = nan(size(corr_m)); %set variables here, get replaced below if pb file exists
    corr_ppb = corr_mpb;
    if ~isempty(pb_file) % if pb file in new file format
        pb_file_list{end+1}=files{pb_file,1}{1};
        data_ol = load(fullfile(matfiles(iF).folder,files{pb_file,1}{1}));
        true_speed_ol = data_ol.vr_data_resampled.velM/opt.TimeBin;
        [vis_flow_pb,sp] = calcSpeed(data_ol.posx,opt);
        velM_pb = gauss_smoothing(true_speed_ol,opt.smoothSigma_time/opt.TimeBin);
        velP_pb = gauss_smoothing(sp,opt.smoothSigma_time/opt.TimeBin);
        fr_pb = calcFRVsTime(good_cells,data_ol,opt);
        corr_mpb=corr(fr_pb',velM_pb');
        corr_ppb=corr(fr_pb',velP_pb);
    end
    
    if isfile(fullfile(matfiles(iF).folder,pb_name)) %pb file old data
        pb_file_list{end+1}=pb_name;
        data_ol = load(fullfile(matfiles(iF).folder,pb_name));
        true_speed_ol = data_ol.true_speed/opt.TimeBin;
        [vis_flow_pb,sp] = calcSpeed(data_ol.posx,opt);
        velM_pb = gauss_smoothing(true_speed_ol,opt.smoothSigma_time/opt.TimeBin);
        velP_pb = gauss_smoothing(sp,opt.smoothSigma_time/opt.TimeBin);
        fr_pb = calcFRVsTime(good_cells,data_ol,opt);
        if isrow(velM_pb)
            velM_pb=velM_pb';
        end
        corr_mpb=corr(fr_pb',velM_pb);
        corr_ppb=corr(fr_pb',velP_pb);
    end
    
    delays = zeros(numel(good_cells),1);
    for iC=1:numel(good_cells)
        %idx = data_out.sp.clu==good_cells(iC);
        %frMat(iC,:)=histcounts(data_out.sp.st(idx),tB);
        [tmp,max_c] = finddelay(fr(iC,:),velM_s,50);
        if max_c<1e-8
            tmp = nan;
        end
        delays(iC)=tmp;
    end
    mf = matfile(fullfile(savepath,matfiles(iF).name),'Writable',true);
    mf.corr_m=corr_m;
    mf.corr_pb = [corr_mpb,corr_ppb];
    mf.delays = delays;
    mf.MM = count_vec;
    mf.spike_times = spike_times;
    mf.region = region_this;
    mf.subregion=subregion_this;
    mf.fr = firing_rate;
    mf.clu_id = CLUIDS{iF};
    %CORR_M = cat(1,CORR_M,corr_m);
    %CORR_PB = cat(1,CORR_PB,[corr_mpb,corr_ppb]);
    %DELAYS = cat(1,DELAYS,delays);
    %MM=cat(1,MM,count_vec);
    %REGION = cat(1,REGION,region_this);
    %MM_R = cat(1,MM_R,count_vec_random);
    %SID = cat(1,SID,ones(numel(good_cells),1)*iF);
    %FR = cat(1,FR,firing_rate);
    
end
%%
%%
MM=[];

SID = [];
MM_R=[];
THETA_POWER = [];
CORR_M=[];
CORR_PB=[];
DELAYS=[];
REGION={};
SUBREGION = {};
FR=[];
savefiles=dir(fullfile(savepath,'*.mat'));

for iF=1:numel(savefiles)
    mf = matfile(fullfile(savefiles(iF).folder,savefiles(iF).name));
    CORR_M = cat(1,CORR_M,mf.corr_m);
    CORR_PB = cat(1,CORR_PB,mf.corr_pb);
    DELAYS = cat(1,DELAYS,mf.delays);
    MM=cat(1,MM,mf.MM);
    REGION = cat(1,REGION,mf.region);
    SUBREGION = cat(1,SUBREGION,mf.subregion);
    %MM_R = cat(1,MM_R,count_vec_random);
    SID = cat(1,SID,ones(numel(mf.corr_m),1)*iF);
    FR = cat(1,FR,mf.fr);
end

%%
REGION=strrep(REGION,'MEC','ENTm');

%% heatmap and average

[uR,~,idxR]=unique(REGION);

MM_smooth = smoothdata(MM,2,'gaussian',5);
MM_ms = MM_smooth-mean(MM_smooth(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
%MM_ms=MM_ms./FR;
for iR=1:numel(uR)
    
    idx = idxR==iR & FR>=1 & isfinite(MM(:,1));
    idx = idx & isfinite(CORR_PB(:,1));
    if nnz(idx)<200
        continue
    end
    figure('Position',[ 680   141   366   837],'Color','w')
    subplot(5,1,[1:4])
    tmp_r = MM_ms(idx,:);
    tmp_d = CORR_M(idx);
    tmp_d = CORR_PB(idx,2);
    %tmp_d = CORR_PB(idx,1);
    %tmp_d = FR(idx);
    [~,sid]=sort(tmp_d);
    imagesc(opt.time_vecs,1:numel(sid),tmp_r(sid,:),[-5 5]);
    xlim([-.5 1.5])
    colormap(brewermap(40,'*RdBu'))
    title(uR{iR})
    subplot(5,1,5)
    mm = nanmean(tmp_r);
    boundedline(opt.time_vecs,nanmean(MM_ms(idx,:)),nanstd(MM(idx,:))/sqrt(nnz(idx)))
    xlim([-.5 1.5])
    title(sprintf('%s, n=%d',uR{iR},nnz(idx) ))
    ylabel('change in FR')
    xlabel('time')
end
%% heatmap of subset
figure
hold on
%regs = {'MEC','VISp','ECT','RSPd'};
regs = {'VISp','RSPv','ENTm'};

MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
MM_ms = smoothdata(MM_ms,2,'gaussian',5);
mm=MM_ms./FR;
for iR=1:numel(regs)
    
    subplot(1,numel(regs),iR)
    idx = strcmp(REGION,regs{iR}) & FR>1 & isfinite(MM_ms(:,1));
    idx = idx & isfinite(CORR_PB(:,1));
    tmp_r = MM_ms(idx,:);
    tmp_d = CORR_M(idx);
    tmp_d = CORR_PB(idx,2);
    
    %tmp_d = FR(idx);
    [~,sid]=sort(tmp_d);
    imagesc(opt.time_vecs,1:numel(sid),tmp_r(sid,:),[-5 5]);
    xlim([-.5 1.5])
    colormap(brewermap(40,'*RdBu'))
    title(regs{iR})
    
end
%% scatter mimsatch and corrM
figure
hold on
regs = {'ENTm','VISp'};
mm_resp = mean(MM(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)-mean(MM(:,opt.time_vecs<-.1),2);
mm=mm_resp./FR;
%mm=mm_resp;
cmap = brewermap(numel(regs),'Set1');
for iR=1:numel(regs)
    idx = strcmp(REGION,regs{iR}) & FR>1;
    
    scatter(CORR_M(idx),mm(idx),35,cmap(iR,:),'.')
end
legend(regs)
grid on
%% scatter corrm ,mm, colored by delay
hold on
regs = {'ENTm','VISp','SCs'};
mm_resp = mean(MM(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)-mean(MM(:,opt.time_vecs<-.1),2);
mm=mm_resp./FR;
for iR=1:numel(regs)
    subplot(2,numel(regs),iR)
    idx = strcmp(REGION,regs{iR}) & FR>1;
    
    scatter(CORR_M(idx),mm(idx),35,DELAYS(idx),'.')
    xlabel('corr_m')
    ylabel('mmresp')
    grid on
    subplot(2,numel(regs),iR+numel(regs))
    scatter(DELAYS(idx),mm(idx),35,CORR_M(idx),'.')
    set(gca,'CLim',[-.5 .5])
end
%% responses per site

regs ={'VISp','MBmot','ECT','ENTm'};
uS = unique(SID);
responses = struct();
mm_resp = mean(MM(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)-mean(MM(:,opt.time_vecs<-.1),2);
MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
%MM_ms = MM_ms./FR;
%mm_resp = mm_resp./FR;
for iR=1:numel(regs)
    responses.(regs{iR}) = [];
    for iS=1:numel(uS)
        idx_this = strcmp(REGION,regs{iR}) & SID==uS(iS);%; & FR>1;
        resp_this = mm_resp(idx_this);
        t=prctile(resp_this,80);
        idx_this = idx_this & mm_resp>t;
        if nnz(idx_this)>5
            mm=mean(MM_ms(idx_this,:));
            responses.(regs{iR})=cat(1,responses.(regs{iR}),mm);
        end
        
    end
    
end


figure
hold on
cmap = brewermap(numel(regs),'Set1');
for iR=1:numel(regs)
    plot(opt.time_vecs,responses.(regs{iR}),'Color',cmap(iR,:))
end
figure('Color',[1 1 1])
hold on
cmap = brewermap(numel(regs),'Set1');
for iR=1:numel(regs)
    plot(opt.time_vecs,mean(responses.(regs{iR})),'Color',cmap(iR,:))
end
xlabel('Time')
ylabel('Change in FR [Hz]')
xlim([-0.1,1])
legend(regs)
%% delays


%set(gca,'CLim',[-1 1])
regs = {'ENTm','VISp'};
figure
hold on
cmap = brewermap(2,'*Set1');
for iR=1:numel(regs)
    subplot(1,3,1)
    hold on
    idx = strcmp(REGION,regs{iR}) & FR>1;
    histogram(DELAYS(idx)/50,15,'normalization','pdf')
    xlabel('Delay fr and run')
    subplot(1,3,2)
    hold on
    scatter(CORR_M(idx),DELAYS(idx)/50,15,cmap(iR,:))
    xlabel('corr run')
    ylabel('delay')
    subplot(1,3,3)
    hold on
    idx = strcmp(REGION,regs{iR}) & CORR_M>0.2;
    histogram(DELAYS(idx)/50,10,'normalization','pdf')
end
for ii=1:3
subplot(1,3,ii)
    legend(regs)
end
%%
figure
hold on
%regs = {'MEC','VISp','ECT','RSPd'};
regs = {'VISp','RSPv','ENTm'};

MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
MM_ms = smoothdata(MM_ms,2,'gaussian',5);
mm_resp = mean(MM(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)-mean(MM(:,opt.time_vecs<-.1),2);
%mm_resp = mm_resp./FR;
%mm=MM_ms./FR;
for iR=1:numel(regs)
    subplot(1,numel(regs),iR)
    idx = strcmp(REGION,regs{iR})  & isfinite(MM_ms(:,1)) & isfinite(CORR_PB(:,1));
    tmp_r = MM_ms(idx,:);
    tmp_d = CORR_M(idx);
    %scatter(CORR_M(idx),CORR_PB(idx,1),35,mm_resp(idx),'.')
    scatter(CORR_PB(idx,2),CORR_PB(idx,1),35,mm_resp(idx),'.')
    colormap(brewermap(20,'*RdYlBu'))
    %colormap jet
    set(gca,'CLIM',[-5 5])
    axis image
    grid on
    title(nnz(idx))
end
%%
figure
IDX = startsWith(REGION,'VISp') & FR>=.5;
imagesc(MM_ms(IDX,:))
%% V1 and MEC
uS=unique(SID);
for iS=1:numel(uS)
    if nnz(startsWith(REGION(SID==uS(iS)),'ECT'))>0 && nnz(startsWith(REGION(SID==uS(iS)),'ENTm'))>0
        figure
        idx_mec = startsWith(REGION,'ENTm') & SID==uS(iS)& FR>=1;
        mm_mec = nanmean(MM_ms(idx_mec,:));
        idx_vis = startsWith(REGION,'ECT') & SID==uS(iS)& FR>=1;
        mm_vis = nanmean(MM_ms(idx_vis,:));
        plot(opt.time_vecs,mm_mec)
        hold on
        plot(opt.time_vecs,mm_vis)
    end
end
%%
for iS=1:numel(uS)
    if nnz(startsWith(REGION(SID==uS(iS)),'ENTm'))>0
        unique(REGION(SID==uS(iS)))
    end
end
%% RSPd and RSPv
uS=unique(SID);
for iS=1:numel(uS)
    if nnz(startsWith(REGION(SID==uS(iS)),'RSPd'))>0 && nnz(startsWith(REGION(SID==uS(iS)),'RSPv'))>0
        figure
        idx_mec = startsWith(REGION,'RSPd') & SID==uS(iS)& FR>=1;
        mm_mec = nanmean(MM_ms(idx_mec,:));
        idx_vis = startsWith(REGION,'RSPv') & SID==uS(iS)& FR>=1;
        mm_vis = nanmean(MM_ms(idx_vis,:));
        plot(opt.time_vecs,mm_mec)
        hold on
        plot(opt.time_vecs,mm_vis)
    end
end
%%
regs = unique(REGION);
%regs = {'MEC','ECT','VISp','RSPd','RSPv','RHP','RSPagl'}
MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
MM_ms = smoothdata(MM_ms,2,'gaussian',5);
MM_ms./FR;
for iR=1:numel(regs)
    idx = strcmp(REGION,regs{iR}) & FR >1;
    if nnz(idx)>50
        figure
        %plot(opt.time_vecs,nanmean(MM_ms(idx,:)))
        boundedline(opt.time_vecs,nanmean(MM_ms(idx,:)),nanstd(MM(idx,:))/sqrt(nnz(idx)))
        title(sprintf('%s, n=%d',regs{iR},nnz(idx)))
        xlim([-.5 2])
    end
end
%%
figure
hold on
%regs = {'MEC','VISp','ECT','RSPd'};
regs = {'ENTm','VISp','RSPd'};

regs = unique(REGION);
MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
MM_ms = smoothdata(MM_ms,2,'gaussian',5);

%mm=MM_ms./FR;
for iR=1:numel(regs)
    
    idx = strcmp(REGION,regs{iR}) & FR>1;
    if nnz(idx)>50
        figure
        tmp_r = MM_ms(idx,:);
        tmp_d = CORR_M(idx);
        %tmp_d = FR(idx);
        [~,sid]=sort(tmp_d);
        imagesc(opt.time_vecs,1:numel(sid),tmp_r(sid,:),[-10 10]);
        xlim([-.5 1.5])
        colormap(brewermap(40,'*RdBu'))
        title(regs{iR})
    end
end