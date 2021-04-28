%new version that uses malcolms approach for shift correction
ops = load_default_opt;
ops.trial_range = [-6:9];
ops.BinWidth = 2;
ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.smoothSigma=ops.smoothSigma_dist;
smoothSigma = ops.smoothSigma/ops.BinWidth;
ops.filter = gausswin(floor(smoothSigma*5/2)*2+1);
ops.filter = ops.filter/sum(ops.filter);
ops.max_lag = 30;
ops.maxLag = ops.max_lag;
%OAK='/oak/stanford/groups/giocomo/';
OAK = '/Volumes/T7';
%OAK = '/Volumes/Crucial X8/';
savepath = fullfile(OAK,'attialex','tbtxcorr_shift_validation');
shiftDir = fullfile(OAK,'attialex','speed_filtered_correctedData_shortidx2');
if ~isfolder(savepath)
    mkdir(savepath)
end
shift_ops = load(fullfile(shiftDir,'parameters.mat'));
shift_ops = shift_ops.ops;

%%
gain = 0.8;
contrast = 100;
regions = {'VISp','RS','MEC'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA_corrected'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end
%%

for iF=1:numel(filenames)
    
    
    [~,sn]=fileparts(filenames{iF});
    data = load(filenames{iF});
    
    shift_data = load(fullfile(shiftDir,[sn '.mat']));
    tmp = shift_data.all_stability;
    stab_idx = tmp>.5;
    tmp_f = shift_data.all_factors;
    tmp_f(~stab_idx)=nan;
    score = nan([size(shift_data.all_factors),2]);
    for iC=1:numel(shift_data.CID)
        n_stable_blocks = sum(stab_idx(:,iC));
        if n_stable_blocks<2
            continue
        end
        stable_blocks = find(stab_idx(:,iC));
        
        for iBlock = 1:numel(stable_blocks)
            dat2 = data;
            this_block = stable_blocks(iBlock);
            mask = true(1,n_stable_blocks);
            mask(iBlock) = false; %look at all other stable blocks
            shift = nanmean(tmp_f(stable_blocks(mask),iC));
            
            trials_this_block = shift_data.start_idx(this_block)+(0:9);
            idx = dat2.sp.clu==shift_data.CID(iC);
            dat2.sp.st(idx)=dat2.sp.st(idx)+shift; %%verify
            frMat = calcTrialFRMat(shift_data.CID(iC),trials_this_block,dat2,shift_ops);
            fr_bl = calcTrialFRMat(shift_data.CID(iC),trials_this_block,data,shift_ops);
            %dc = nanmean(1-pdist(fr_bl(:,shift_ops.idx),'correlation'))-nanmean(1-pdist(frMat(:,shift_ops.idx),'correlation'));
            score(this_block,iC,1)=nanmean(1-pdist(fr_bl(:,shift_ops.idx),'correlation'));
            score(this_block,iC,2)=nanmean(1-pdist(frMat(:,shift_ops.idx),'correlation'));
        end
        
    end
    if isfield(data.anatomy,'parent_shifted')
        reg = data.anatomy.parent_shifted;
    else
        reg = data.anatomy.cluster_parent;
    end
    if iscolumn(reg)
        reg = reg';
    end
    
    reg=reg(data.sp.cgs==2);
    data_out = matfile(fullfile(savepath,sprintf('%s',sn)),'Writable',true);
    
    data_out.region = reg;
    
    data_out.factors = shift_data.all_factors;
    data_out.CID = shift_data.CID;
    data_out.score = score;
end
%%
OAK = '/Volumes/T7';
savepath = fullfile(OAK,'attialex','tbtxcorr_shift_validation');
files = dir(fullfile(savepath,'*.mat'));
CORR = [];
REG = {};
DT = [];
opt.pval_cutoff = 0.01;
opt.min_prom = 0.1;
dist_cells = dist_tuning.pval<opt.pval_cutoff & dist_tuning.prom>opt.min_prom;
for iF=1:numel(files)
    data = load(fullfile(files(iF).folder,files(iF).name));
    els = strsplit(files(iF).name,'_');
    mouse = els{1};
    date = els{2};
    dist_tuned_this = nan(size(data.CID))';
    for iC=1:numel(data.CID)
        identifier = sprintf('%s_%s_c%d',mouse,date,data.CID(iC));
        tmp = dist_cells(strcmp(dist_tuning.UniqueID,identifier));
        if ~isempty(tmp)
            dist_tuned_this(iC) = tmp;
        end
    end
    cv_corr = squeeze(nanmean(data.score));
    CORR = cat(1,CORR,cv_corr);
    REG = cat(1,REG,data.region');
    DT = cat(1,DT,dist_tuned_this);
end

REG(startsWith(REG,'VISpm'))={'VIS'};
REG(startsWith(REG,'RSPagl'))={'RSA'};

%%
regions = {'MEC','VISp','RSP'};
figure
for iR = 1:3
    tmp = CORR(:,1)-CORR(:,2);
    if iR==1
    idx = startsWith(REG,regions{iR}) & isfinite(tmp) & DT==1;
    else
        idx = startsWith(REG,regions{iR}) & isfinite(tmp);
    end
    subplot(1,3,iR)
    
    histogram(-1*tmp(idx),-.31:0.02:.31)
    title(regions{iR})
    [h,p]=ttest(tmp(idx));
    xlabel(sprintf('%.4f, p=%.2e',mean(-1*tmp(idx)),p))
end
%% fractional change

regions = {'MEC','VISp','RSP'};
figure
for iR = 1:3
    tmp = (CORR(:,2)-CORR(:,1))./(CORR(:,1));
if iR==1
    idx = startsWith(REG,regions{iR}) & isfinite(tmp) & DT==0;
    else
        idx = startsWith(REG,regions{iR}) & isfinite(tmp);
end
subplot(1,3,iR)
    
    histogram(tmp(idx),[-.55:0.05:0.55])
    title(regions{iR})
    [h,p]=ttest(tmp(idx));
    %tmp = CORR(:,2)./(CORR(:,1));
    xlabel(sprintf('%.4f, p=%.2e',mean(1*tmp(idx)),p))
end