function [acg,spacing,pxx,f]=calculate_spatial_acg(spike_distance,binedges,gauss_filter,time_per_bin)

[aa,~]=histcounts(spike_distance,binedges);

%moothing
firing_rate = aa./time_per_bin;
if sum(isnan(firing_rate))>0
    firing_rate = interp1(find(~isnan(firing_rate)),firing_rate(~isnan(firing_rate)),1:numel(firing_rate));
end
aa = conv(firing_rate,gauss_filter,'same');

[acg,spacing] = xcorr(aa,2000,'coeff');
%[pxx,f]=pwelch(aa,[],[],[],1/2);
pxx=0;
f=0;
end