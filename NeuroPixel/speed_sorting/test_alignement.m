function [correlation_shifted,correlation_noshift]=test_alignement(data,ops,factors)

[speed,speed_raw]=calcSpeed(data.posx,ops);
test_trials = ops.trials;
%% test alignement
good_cells = data.sp.cids(data.sp.cgs==2);
%create trial map that only contains numbers for trials to be included
trialMap = nan(1,numel(data.trial_gain));

cntr = 1;
for iT =1:numel(data.trial_gain)
    if ismember(iT,test_trials)
        trialMap(iT)=cntr;
        cntr=cntr+1;
    end
end

%recreate data.trial map that resets numbers of included trials and sets
%all else to nan
trial_sorted = nan(size(data.trial));
for iT=1:numel(trial_sorted)
    trial_sorted(iT)=trialMap(data.trial(iT));
end
nC=numel(good_cells);
correlation_noshift = nan(1,nC);
correlation_shifted = nan(1,nC);
idx = triu(true(numel(test_trials)),1);
for iC=1:nC
    factor_this = factors(iC);
    good_idx = ismember(data.sp.clu,good_cells(iC));
    clu_tmp = data.sp.clu(good_idx);
    st_tmp = data.sp.st(good_idx);
    [uClu,~,clus]=unique(clu_tmp);
    
        
    [spMapShifted]=shiftAllMapsByFactor(ops,clus,st_tmp,1,data.posx,data.post,trial_sorted,speed,factor_this);
    cc=corr(spMapShifted(:,ops.idx)');
    correlation_shifted(iC)=nanmean(cc(idx));
    
    [spMapBL]=shiftAllMapsByFactor(ops,clus,st_tmp,1,data.posx,data.post,trial_sorted,speed,0);
    cc=corr(spMapBL(:,ops.idx)');
    correlation_noshift(iC)=nanmean(cc(idx));
end
end