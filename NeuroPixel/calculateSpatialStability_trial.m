function rho = calculateSpatialStability_trial(idx,posx,trial,p)
% function to calculate spatial stability by correlating firing rate in
% first half of trials with second half of trials
%
% modified 6/6/18 MGC
% all time bins are now equal length
% has not been fully tested yet
% 
% inputs:
%     idx: spike indices
%     posx: position
%     p: params struct (has bin size, track length, etc)
% outputs:
%     rho: spatial stability
    
median_trial = median(unique(trial));   
fr1 = calculateSmoothedFiringRate(idx(trial(idx)<=median_trial),posx(trial<=median_trial),p);
fr2 = calculateSmoothedFiringRate(idx(trial(idx)>median_trial)-sum(trial<=median_trial),posx(trial>median_trial),p);
rho = corr(fr1',fr2');

end