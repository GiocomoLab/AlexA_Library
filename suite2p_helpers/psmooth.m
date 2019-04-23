function [trace]=psmooth(trace,win,fps)
% implementiert tanks prctile smoothing algorithm (dombeck2007)
% subtract 8 percentile level as calculated on +/- 5 seconds window 
% from relative fluorescence signal

prctile_level=8;

window_size=2*floor(min(length(trace),win*fps)/2);

if size(trace,2)>1
    trace=trace';
end
% tot_median=median(trace);
first_win_media=median(trace(1:window_size));
last_win_media=median(trace(end-window_size+1:end));


lt=length(trace);
X = [ones(window_size/2,1)*first_win_media; trace; ones(window_size/2,1)*last_win_media];

indr = (0:window_size-1)';
indc = 1:lt;

ind = indc(ones(1,window_size),1:lt) + indr(:,ones(1,lt));
xx = sort(X(ind),1);
prctile_average=xx(round(window_size*prctile_level/100),:)';

trace=trace-prctile_average+median(prctile_average);