function firing_rate = calculateSmoothedFiringRate(idx,posx,p)
% calculates smoothed firing rate on linear track
% Malcolm Campbell 5/21/15
%
% modified 6/6/18 MGC
% all time bins are now equal length
%
% inputs:
%     idx: spike indices
%     posx: positions
%     p: params (spatial bin size, etc)
% outputs:
%     firing_rate: smoothed firing rate over position

% divide spike counts by occupancy
binedges = p.TrackStart:p.SpatialBin:p.TrackEnd;
time_per_bin = histcounts(posx, binedges);
time_per_bin = time_per_bin * p.TimeBin;
firing_rate = histcounts(posx(idx), binedges);
firing_rate = firing_rate./time_per_bin;

% interpolate missing values
if sum(isnan(firing_rate))>0
    firing_rate = interp1(find(~isnan(firing_rate)),firing_rate(~isnan(firing_rate)),1:numel(firing_rate));
end

% gaussian filter for smoothing
smoothSigma = p.SmoothSigmaFR/p.SpatialBin;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);

% smooth firing rate
firing_rate_smooth = conv(repmat(firing_rate,1,3),gauss_filter,'same');
firing_rate = firing_rate_smooth(numel(firing_rate)+1:numel(firing_rate)*2);

end
