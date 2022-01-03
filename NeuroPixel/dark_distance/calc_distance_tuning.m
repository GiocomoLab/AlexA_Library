function data_out = calc_distance_tuning(data_in,good_cells,opt)
    data_in.trial = [];
    fr = calcFRVsDist(good_cells,1:max(data_in.trial),data_in,opt);
    fr_zscore = zscore(fr,0,2);
    
    peak_all = nan(numel(good_cells),1);
    peak_loc_all = nan(numel(good_cells),1);
    peak_prom_all = nan(numel(good_cells),1);
    xplot = 0:opt.SpatialBin:opt.max_lag;

    xcs = nan(numel(good_cells),numel(xplot));

    for cIdx = 1:numel(good_cells)
        y = fr_zscore(cIdx,:);
        xc = xcorr(y,opt.max_lag/opt.SpatialBin,"coeff");
        xc = xc(opt.max_lag/opt.SpatialBin+1:end);
        xcs(cIdx,:)=xc;
        [peaks,locs,~,prominence] = findpeaks(xc,'SortStr','descend','NPeaks',1);
        if ~isempty(peaks)
            peak_all(cIdx) = peaks;
            peak_loc_all(cIdx) = xplot(locs);
            peak_prom_all(cIdx) = prominence;
        end
    end
    
    peak_shuf = nan(numel(good_cells),opt.num_shuf);
    for sIdx = 1:opt.num_shuf
        fr_shuf = calcFRVsDist_shuf(good_cells,1:max(data_in.trial),data_in,opt);
        fr_shuf_zscore = zscore(fr_shuf,[],2);
        for cIdx = 1:numel(good_cells)
            y = fr_shuf_zscore(cIdx,:);
            xc = xcorr(y,opt.max_lag/opt.SpatialBin,"coeff");
            xc = xc(opt.max_lag/opt.SpatialBin+1:end);
            peaks = findpeaks(xc,'SortStr','descend','NPeaks',1);
            if ~isempty(peaks)
                peak_shuf(cIdx,sIdx) = max(peaks);
            end
        end
    end
    
    pval = nan(numel(good_cells),1);
    for cIdx = 1:numel(good_cells)
        pval(cIdx) = sum(peak_shuf(cIdx,:)>=peak_all(cIdx))/opt.num_shuf;
    end
    pval(isnan(peak_all) | sum(isnan(peak_shuf),2)>opt.num_shuf/3) = nan;
    
   
    data_out.peak_all = peak_all;
    data_out.peak_loc_all = peak_loc_all;
    data_out.peak_prom_all = peak_prom_all;
    data_out.xcorrs = xcs;
    data_out.pvals = pval;
    data_out.good_cells = good_cells;
    
end