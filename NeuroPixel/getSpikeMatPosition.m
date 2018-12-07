function [av_spike_mat,dwell_time] = getSpikeMatPosition(spike_times,spike_clusters,posx,post,varargin)

p = inputParser;
defaultEdges = 0:401;

addParameter(p,'edges',defaultEdges);
addParameter(p,'trial_vec',[]);
addParameter(p,'trial_selection',[]);
parse(p,varargin{:});

edges=p.Results.edges;

%discretizePosition

discrete_pos=discretize(posx,edges);
dwell_time=histcounts(posx,edges);

n_units = max(spike_clusters)+1;

%max_t=ceil(max(sp.st)*1000);

av_spike_mat=zeros(n_units,length(edges)-1);
post_idx=nearestpoint(spike_times,post);
spike_loc=discrete_pos(post_idx);

for iC=1:length(spike_times)
    sp_id=spike_clusters(iC)+1;
    loc=spike_loc(iC);
    
    av_spike_mat(sp_id,loc)=av_spike_mat(sp_id,loc)+1;

end




end


