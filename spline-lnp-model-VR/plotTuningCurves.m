function [ fighandle ] = plotTuningCurves(occ_ax,tc_ax,spiketrain,var_name,tuning_var,nbins,minVal,maxVal,TimeBin)
    
    
    [pos_tuning_curve,pos_occupancy,bins] = compute_1d_tuning_curve(tuning_var,spiketrain,nbins,minVal,maxVal);
    xbincent=0.5 * (bins(1:end-1) + bins(2:end));
        axes(occ_ax)

    plot(xbincent,pos_occupancy.*TimeBin)
    title([var_name ' occupancy'])
    axis tight
    
    axes(tc_ax)
    plot(xbincent,pos_tuning_curve./TimeBin,'r')
    title([var_name])
    ylabel('spikes/s')
    axis tight
    
end

