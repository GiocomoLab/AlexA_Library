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
savepath = fullfile(OAK,'attialex','tbtxcorr_constant_shift');
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
shifts = [-0.12,-0.06,0.01];
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA_corrected'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end
%%

for iF=38:numel(filenames)
    
    
    [~,sn]=fileparts(filenames{iF});
    data = load(filenames{iF});
    if isfield(data.anatomy,'parent_shifted')
        reg = data.anatomy.parent_shifted;
    else
        reg = data.anatomy.cluster_parent;
    end
    if iscolumn(reg)
        reg = reg';
    end
    
    reg=reg(data.sp.cgs==2);
    
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
        
        if startsWith(reg{iC},'VISp')
            shift = -.12;
        elseif startsWith(reg{iC},'RS')
            shift = -.06;
        elseif startsWith(reg{iC},'MEC')
            shift = 0.01;
        else
            shift = 0;
        end
        dat2 = data;
        idx = dat2.sp.clu==shift_data.CID(iC);

        dat2.sp.st(idx)=dat2.sp.st(idx)+shift; %%verify
        
        for iBlock = 1:numel(stable_blocks)
            this_block = stable_blocks(iBlock);            
            trials_this_block = shift_data.start_idx(this_block)+(0:9);
            frMat = calcTrialFRMat(shift_data.CID(iC),trials_this_block,dat2,shift_ops);
            fr_bl = calcTrialFRMat(shift_data.CID(iC),trials_this_block,data,shift_ops);
            %dc = nanmean(1-pdist(fr_bl(:,shift_ops.idx),'correlation'))-nanmean(1-pdist(frMat(:,shift_ops.idx),'correlation'));
            score(this_block,iC,1)=nanmean(1-pdist(fr_bl(:,shift_ops.idx),'correlation'));
            score(this_block,iC,2)=nanmean(1-pdist(frMat(:,shift_ops.idx),'correlation'));
        end
        
    end
    
    data_out = matfile(fullfile(savepath,sprintf('%s',sn)),'Writable',true);
    
    data_out.region = reg;
    
    data_out.factors = shift_data.all_factors;
    data_out.CID = shift_data.CID;
    data_out.score = score;
end
%%
OAK = '/Volumes/T7';
savepath = fullfile(OAK,'attialex','tbtxcorr_shift_validation');
savepath_constant = fullfile(OAK,'attialex','tbtxcorr_constant_shift');

files = dir(fullfile(savepath,'*.mat'));
files_constant = dir(fullfile(savepath_constant,'*.mat'));
CORR = [];
CORR_CONSTANT = [];
REG = {};
for iF=1:numel(files)
    data = load(fullfile(files(iF).folder,files(iF).name));
    data_constant = load(fullfile(files_constant(iF).folder,files_constant(iF).name));
    cv_corr = squeeze(nanmean(data.score));
    corr_const = squeeze(nanmean(data_constant.score));
    CORR = cat(1,CORR,cv_corr);
    CORR_CONSTANT = cat(1,CORR_CONSTANT,corr_const);
    REG = cat(1,REG,data.region');
end

REG(startsWith(REG,'VISpm'))={'VIS'};
REG(startsWith(REG,'RSPagl'))={'RSA'};

%%
regions = {'MEC','VISp','RSP'};
figure
for iR = 1:3
    tmp = CORR(:,2);
    tmp_c = CORR_CONSTANT(:,2);
    idx = startsWith(REG,regions{iR}) & isfinite(tmp);
    subplot(2,3,iR)
    
    scatter(tmp(idx),tmp_c(idx))
    title(regions{iR})
    [h,p]=ttest(tmp(idx)-tmp_c(idx));
    xlabel(sprintf('%.4f, p=%.2e',mean(tmp(idx)./tmp_c(idx)),p))
    subplot(2,3,iR+3)
    pref = (tmp-tmp_c)./(tmp+tmp_c);
    histogram(pref(idx),[-.25:0.025:0.25])
end
%%
regions = {'MEC','VISp','RSP'};
figure
for iR = 1:3
    tmp = CORR(:,2);
    tmp_c = CORR_CONSTANT(:,2);
    idx = startsWith(REG,regions{iR}) & isfinite(tmp);
    subplot(2,3,iR)
    
    scatter(tmp(idx),tmp_c(idx))
    title(regions{iR})
    [h,p]=ttest(tmp(idx)-tmp_c(idx));
    xlabel(sprintf('%.4f, p=%.2e',mean(tmp(idx)./tmp_c(idx)),p))
    subplot(2,3,iR+3)
    pref = (tmp)./(tmp_c);
    histogram(pref(idx),[.75:0.025:1.25])
end