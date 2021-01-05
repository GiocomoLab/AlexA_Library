function count_vec = spikeStruct2countVec(trial_vec, good_cells,opt)
%trial_vec =cat(1,spike_times_struct{:});
count_vec = zeros(numel(good_cells),numel(opt.time_bins)-1);
for iC=1:numel(good_cells)
    idx = trial_vec(:,2)==good_cells(iC);
    [spike_count]=histcounts(trial_vec(idx,1),opt.time_bins);
    count_vec(iC,:)=spike_count;
    
end
