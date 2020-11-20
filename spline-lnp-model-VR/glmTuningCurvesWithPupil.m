mm_resp = nanmean(mm(:,opt.time_vecs> 0.1 & opt.time_vecs<1.0),2)-nanmean(mm(:,opt.time_vecs>-1 & opt.time_vecs<-0.1),2);

good_cells = data.sp.cids(data.sp.cgs==2);
opt = load_mismatch_opt;
figure
for cellIDX = 1:numel(good_cells)
    if nnz(glmData(cellIDX).bestModels>=3)
        

spike_t = data.sp.st(data.sp.clu==good_cells(cellIDX));
    %[~,~,spike_idx] = histcounts(spike_t,post);
    
    [spiketrain,~,spikeIDX] = histcounts(spike_t,data.post);
fr = gauss_smoothing(spiketrain,opt.smoothSigma_time/opt.TimeBin);
for iV=1:5
    if ismember(iV,glmData(cellIDX).bestModels)
    subplot(2,5,iV)
    plot(glmData(cellIDX).tuning_curves{iV})
    
    end
    subplot(2,5,6)
    plot(mm(cellIDX,:))
end
pause
clf
    end
end