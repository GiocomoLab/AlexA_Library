
%matfiles = dir('Z:\giocomo\attialex\NP_DATA\mismatch\*mismatch*.mat');

%matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA_corrected/*mismatch*.mat');
%matfiles = dir('/Volumes/Samsung_T5/attialex/NP_DATA/*mismatch*.mat');
%matfiles = dir('/Users/attialex/NP_DATA_2/*mismatch*.mat');
%matfiles = dir('/Users/attialex/mismatch/*mismatch*.mat');
%matfiles = matfiles(~cellfun(@(x) contains(x,'tower'), {matfiles.name}));
matfiles = dir('F:\Alex\new_2\*_MMCueRich_*');
savepath = 'F:\Alex\mmCueRich_out';
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
            subreg = data_out.anatomy.region_shifted;
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
    data_out.trial = cumsum([0;diff(data_out.posx)<-100]);

opt.TrackEnd=round(max(data_out.posx/10))*10;
opt.track_length=opt.TrackEnd;
opt.max_lag = 30;
[corrMat,frMat] = trialCorrMat(good_cells,0:max(data_out.trial)-1,data_out,opt);
    
        
    
    
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
    
    
    corr_mpb = nan(size(corr_m)); %set variables here, get replaced below if pb file exists
    corr_ppb = corr_mpb;
   
    
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
    mf.corrMat = corrMat;
    mf.frMat = frMat;
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
STAB = [];
for iF=1:numel(savefiles)
    mf = matfile(fullfile(savefiles(iF).folder,savefiles(iF).name));
    CORR_M = cat(1,CORR_M,mf.corr_m);
    
    DELAYS = cat(1,DELAYS,mf.delays);
    MM=cat(1,MM,mf.MM);
    REGION = cat(1,REGION,mf.region);
    SUBREGION = cat(1,SUBREGION,mf.subregion);
    %MM_R = cat(1,MM_R,count_vec_random);
    SID = cat(1,SID,ones(numel(mf.corr_m),1)*iF);
    FR = cat(1,FR,mf.fr);
    corrMat = mf.corrMat;
    STAB = cat(1,STAB,nanmean(nanmean(corrMat,2),3));
end
%%
[uR,~,idxR]=unique(REGION);

MM_smooth = smoothdata(MM,2,'gaussian',5);
MM_ms = MM_smooth-mean(MM_smooth(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
%MM_ms=MM_ms./FR;
for iR=1:numel(uR)
    
    idx = idxR==iR & FR>=0 & isfinite(MM(:,1));
    %idx = idx & isfinite(CORR_PB(:,1));
    if nnz(idx)<50
        continue
    end
    figure('Position',[ 680   141   366   837],'Color','w')
    subplot(6,1,[1:4])
    tmp_r = MM_ms(idx,:);
    tmp_d = CORR_M(idx);
    %tmp_d = CORR_PB(idx,2);
    %tmp_d = CORR_PB(idx,1);
    tmp_d = STAB(idx);
    [~,sid]=sort(tmp_d);
    imagesc(opt.time_vecs,1:numel(sid),tmp_r(sid,:),[-5 5]);
    xlim([-.5 1.5])
    colormap(brewermap(40,'*RdBu'))
    title(uR{iR})
    subplot(6,1,5)
    mm = nanmean(tmp_r);
    boundedline(opt.time_vecs,nanmean(MM_ms(idx,:)),nanstd(MM(idx,:))/sqrt(nnz(idx)))
    xlim([-.5 1.5])
    title(sprintf('%s, n=%d',uR{iR},nnz(idx) ))
    ylabel('change in FR')
    xlabel('time')
    subplot(6,1,6)
    MM_tmp = MM(idx,:);
    MM_tmp = smoothdata(MM_tmp,2,'gaussian',5);
    tmp_ms = MM_tmp-mean(MM_tmp(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
    idx = sid(1:round(0.2*numel(sid)));
    boundedline(opt.time_vecs,nanmean(tmp_ms(idx,:)),nanstd(MM_tmp(idx,:))/sqrt(nnz(idx)),'alpha')
    idx = sid(numel(sid)-round(0.2*numel(sid)):numel(sid));
    boundedline(opt.time_vecs,nanmean(tmp_ms(idx,:)),nanstd(MM_tmp(idx,:))/sqrt(nnz(idx)),'cmap',[1 0 0],'alpha')
    xlim([-.5 1.5])
    ylabel('change in FR')
    xlabel('time')
    
    
end
%%

mm_resp = mean(MM(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)-mean(MM(:,opt.time_vecs<-.1),2);

[uR,~,idxR]=unique(REGION);

MM_smooth = smoothdata(MM,2,'gaussian',5);
MM_ms = MM_smooth-mean(MM_smooth(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
%MM_ms=MM_ms./FR;
for iR=1:numel(uR)
    idx = idxR==iR & FR>=1 & isfinite(MM(:,1));
    figure
    scatter(STAB(idx),mm_resp(idx),'k')
    title(uR{iR})
end
