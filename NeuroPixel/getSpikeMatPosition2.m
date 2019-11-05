function [spMapN,occupancy,good_cells,edges] = getSpikeMatPosition2(spike_times,spike_clusters,posx,post,varargin)

p = inputParser;
defaultEdges = 0:2:400;

addParameter(p,'edges',defaultEdges);
addParameter(p,'trial_vec',[]);
addParameter(p,'trial_selection',[]);
addParameter(p,'max_clust',[]);
parse(p,varargin{:});

edges=p.Results.edges;
bin_centers = edges(1:end-1)+edges(2:end);
bin_centers = bin_centers/2;
nBins = numel(bin_centers);
smooth_fact = 2;
good_cells = unique(spike_clusters);
spMap = zeros(nBins,numel(good_cells));
for iC=1:numel(good_cells)

spike_t = spike_times(spike_clusters==good_cells(iC));
[~,~,spike_idx] = histcounts(spike_t,post);
fr_x=zeros(nBins,1);
for iBin = 1:nBins
    bin_c = bin_centers(iBin);
    dist = posx(spike_idx)-bin_c; %distance of each spike from bin center
    dist(dist>200)=dist(dist>200)-400;
    dist = dist/smooth_fact;
    dist = exp(-1*dist.^2);
    fr_x(iBin) = sum(dist);
end


spMap(:,iC)=(fr_x);

end

occupancy = zeros(1,nBins);
for iBin = 1:nBins
    temp = posx - bin_centers(iBin);
    temp = temp/smooth_fact;
    temp = temp.^2;
    occupancy(iBin) = sum(exp(-1*temp));
end
occupancy = occupancy/50;
spMapN = bsxfun(@rdivide,spMap,occupancy');

