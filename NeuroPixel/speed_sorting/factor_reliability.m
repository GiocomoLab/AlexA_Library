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
savepath = fullfile(OAK,'attialex','tbtxcorr_factor_reliability');
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
PV = [];
REG = {};

for iF=1:numel(filenames)
    
    [~,sn]=fileparts(filenames{iF});
    %data = load(filenames{iF});
    
    shift_data = load(fullfile(shiftDir,[sn '.mat']));
    tmp = shift_data.all_stability;
    stab_idx = tmp>.5;
    tmp_f = shift_data.all_factors;
    %tmp_f(~stab_idx)=nan;
    p_vals = nan(size(shift_data.CID,2),4);
    for iC=1:numel(shift_data.CID)
        n_stable_blocks = sum(stab_idx(:,iC));
        if n_stable_blocks<2
            continue
        end
        f = mean(tmp_f(stab_idx(:,iC),iC));
        if f==0
            continue
        end
        [~,p_vals(iC,1)] = ttest(tmp_f(stab_idx(:,iC)),iC);
        [~,p_vals(iC,2)]=ttest(tmp_f(:,iC));
        
        t = sign(f);
        consistency = nnz(sign(tmp_f(stab_idx(:,iC),iC))==t)/n_stable_blocks;
        p_vals(iC,3)=consistency;
        if consistency == 0
            keyboard
        end
        p_vals(iC,4)=n_stable_blocks;
    end
    
    PV = cat(1,PV,p_vals);
    REG = cat(2,REG,shift_data.region);
    %data_out = matfile(fullfile(savepath,sprintf('%s',sn)),'Writable',true);
    
    %data_out.region = shift_data.region;
    
    %data_out.factors = p_vals;
    %data_out.CID = shift_data.CID;
end
%%
REG(startsWith(REG,'VISpm'))={'VIS'};
REG(startsWith(REG,'RSPagl'))={'RSA'};

%%
regions = {'MEC','VISp','RSP'};
figure
for iR = 1:3
    tmp = PV(:,3);
    if iR==1
    idx = startsWith(REG,regions{iR}) & isfinite(tmp)';%; & DT==1;
    else
        idx = startsWith(REG,regions{iR}) & isfinite(tmp)';
    end
    subplot(2,3,iR)
    
    histogram(tmp(idx),[0:0.1:1])
    title(regions{iR})
    xlabel(sprintf('%.4f',mean(tmp(idx))))
    subplot(2,3,iR+3)
    scatter(PV(idx,4),tmp(idx))
    xlabel('nstable')
    ylabel('consistency')
end

%%
counts = PV(isfinite(PV(:,4)),4);
max_rep=10000;
C = zeros(max_rep,2);
for i_rep = 1:max_rep
    %n_stable = randi([2,20]);
    n_stable = randsample(counts,1);
    dist = -.3+0.6*rand(n_stable,1);
    f=mean(dist);
    t = sign(f);
    consistency= nnz(sign(dist)==t)/n_stable;
    
    C(i_rep,1)=n_stable;
    C(i_rep,2)=consistency;
end
figure
scatter(C(:,1),C(:,2))

%%
regions = {'MEC','VISp','RSP'};
figure
for iR = 1:3
    tmp = PV(:,3);
    if iR==1
    idx = startsWith(REG,regions{iR}) & isfinite(tmp)';%; & DT==1;
    else
        idx = startsWith(REG,regions{iR}) & isfinite(tmp)';
    end
    
    counts = PV(idx,4);
max_rep=10000;
C = zeros(max_rep,2);
for i_rep = 1:max_rep
    %n_stable = randi([2,20]);
    n_stable = randsample(counts,1);
    dist = -.3+0.6*rand(n_stable,1);
    f=mean(dist);
    t = sign(f);
    consistency= nnz(sign(dist)==t)/n_stable;
    
    C(i_rep,1)=n_stable;
    C(i_rep,2)=consistency;
end
    
    subplot(2,3,iR)
    hold on
    histogram(C(:,2),[0:0.1:1],'Normalization','probability')
    histogram(tmp(idx),[0:0.1:1],'Normalization','probability')
    title(regions{iR})
    xlabel(sprintf('%.4f, %.3e',mean(tmp(idx)),ranksum(C(:,2),tmp(idx))))
    legend({'random','data'},'Location','NW')
    subplot(2,3,iR+3)
    scatter(PV(idx,4),tmp(idx))
    xlabel('nstable')
    ylabel('consistency')
end