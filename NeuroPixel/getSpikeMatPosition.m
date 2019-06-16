function [av_spike_mat,dwell_time] = getSpikeMatPosition(spike_times,spike_clusters,posx,post,varargin)

p = inputParser;
defaultEdges = 0:401;

addParameter(p,'edges',defaultEdges);
addParameter(p,'trial_vec',[]);
addParameter(p,'trial_selection',[]);
addParameter(p,'max_clust',[]);
parse(p,varargin{:});

edges=p.Results.edges;
max_clust = p.Results.max_clust;
%discretizePosition

discrete_pos=discretize(posx,edges);
dwell_time=histcounts(posx,edges);

if isempty(max_clust)
n_units = max(spike_clusters)+1;
else
    n_units = max_clust+1;
end

%max_t=ceil(max(sp.st)*1000);

av_spike_mat=zeros(n_units,length(edges)-1);
%post_idx=nearestpoint(spike_times,post);
[~,~,post_idx] = histcounts(spike_times,post);
spike_loc=discrete_pos(post_idx);

for iC=1:length(spike_times)
    sp_id=spike_clusters(iC)+1;
    loc=spike_loc(iC);
    
    av_spike_mat(sp_id,loc)=av_spike_mat(sp_id,loc)+1;

end




end


