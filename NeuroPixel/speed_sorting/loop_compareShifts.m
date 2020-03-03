ops.factors = -.25:0.01:.25;
ops.BinWidth =2;
ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.n_preceeding = 18;
ops.trials = 3:20;
ops.TimeBin = 0.02;
ops.idx = [10:ops.BinWidth:390]/ops.BinWidth;% in bins
fi = gausswin(5);
fi=fi'/sum(fi);
ops.filter = fi;
ops.bl_pre = 1:ops.n_preceeding;
ops.gain_trials = ops.n_preceeding+[1:4];
ops.similarity_trials = ops.n_preceeding+[-5:4];
ops.plotfig = false;
%%
gain = 0.8;
contrast = 100;
region = 'VISp';
[filenames,triggers] = getFilesCriteria(region,contrast,gain,'F:/NP_DATA');
%%
p=gcp('nocreate');
if isempty(p)
    parpool();
end

%%
output = cell(numel(filenames),1);
parfor iF=1:numel(filenames)
    data = load(filenames{iF});
    ops_local = ops;
    data_out=cell(2,1);
    for iRep = 1:min(2,numel(triggers{iF}))
        preceeding_trials = triggers{iF}(iRep)+[-ops_local.n_preceeding:-1];
        if ~all(data.trial_gain(preceeding_trials)==1 & data.trial_contrast(preceeding_trials)==contrast)
            continue
        end
        ops_local.trials = triggers{iF}(iRep)+[-ops_local.n_preceeding:9];
        data_out{iRep}=compare_speedShift_gainShift(data,ops_local);
        test_trials = triggers{iF}(iRep)+[4:13];
        ops_local.trials=test_trials;
        [correlation_shifted,correlation_noshift]=test_alignement(data,ops_local,data_out{iRep}.factors,data_out{iRep}.CID);
        data_out{iRep}.correlation_shifted = correlation_shifted;
        data_out{iRep}.correlation_noshift = correlation_noshift;
        
    end
    output{iF}=data_out;
end
%% across sites
allM=[];
allS=[];
allMHat = [];
allSHat = [];
region = 'VISp';
for iF=1:numel(output)
    if isempty(output{iF})
        continue
    end
    
    for iRep =1:2
        if isempty(output{iF}{iRep})
            continue
        end
        
        data_out = output{iF}{iRep};
        dcorr = data_out.correlation_shifted - data_out.correlation_noshift;
        idx_region  = startsWith(data_out.region,region);
        idx = data_out.similarity>.4 & startsWith(data_out.region,region)' & dcorr'>0;
        if nnz(idx)>2 && nnz(idx)/nnz(idx_region)>.2

            tmp = squeeze(nanmean(data_out.corrMat(idx,:,:)));  
            tmpS = squeeze(nanmean(data_out.shiftMat(idx,:,:)));
            allM=cat(3,allM,tmp);
            allS = cat(3,allS,tmpS);
            allMHat = cat(3,allMHat,squeeze(nanmean(data_out.corrMatHat(idx,:,:))));
            allSHat = cat(3,allSHat,squeeze(nanmean(data_out.shiftMatHat(idx,:,:))));
            
        end
    end
end
figure
trials2show = 18+[-5:10];
subplot(2,1,1)
imagesc(nanmean(allS(trials2show,trials2show,:),3),[-3 3])
axis image
subplot(2,1,2)
imagesc(nanmean(allSHat(trials2show,trials2show,:),3),[-3 3])
axis image

%%

%% across sites
allM=[];
allS=[];
allMHat = [];
allSHat = [];
region = 'VISp';
for iF=1:numel(output)
    if isempty(output{iF})
        continue
    end
    
    for iRep =1:2
        if isempty(output{iF}{iRep})
            continue
        end
        
        data_out = output{iF}{iRep};
        dcorr = data_out.correlation_shifted - data_out.correlation_noshift;
        idx_region  = startsWith(data_out.region,region);
        idx = data_out.similarity>0.4 & startsWith(data_out.region,region)' & dcorr'>.055 & data_out.factors'<0;
        if nnz(idx)<2
            continue
        end

            tmp = squeeze((data_out.corrMat(idx,:,:)));  
            tmpS = squeeze((data_out.shiftMat(idx,:,:)));
            allM=cat(1,allM,tmp);
            allS = cat(1,allS,tmpS);
            allMHat = cat(1,allMHat,squeeze((data_out.corrMatHat(idx,:,:))));
            allSHat = cat(1,allSHat,squeeze((data_out.shiftMatHat(idx,:,:))));
            
        
    end
end

figure
subplot(2,2,1)
imagesc(squeeze(nanmean(allM(:,trials2show,trials2show),1)),[0 .75])
axis image
subplot(2,2,2)
imagesc(squeeze(nanmean(allMHat(:,trials2show,trials2show),1)),[0 .75])
axis image

subplot(2,2,3)
imagesc(squeeze(nanmean(allS(:,trials2show,trials2show),1)),[-2 2])
axis image
subplot(2,2,4)
imagesc(squeeze(nanmean(allSHat(:,trials2show,trials2show),1)),[-2 2])
axis image