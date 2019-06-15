function hax=raster_map(dataset,cid,varargin)
spike_id = dataset.sp.clu==cid;
    spike_t = dataset.sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,dataset.post);
     
    if numel(varargin)==0
        figure;
        hax=axes();
        col = 'b';
    else
        hax=varargin{1};
        col = varargin{2};
    end
    axes(hax)
    scatter(dataset.posx(spike_idx),dataset.trial(spike_idx),2,col)
end