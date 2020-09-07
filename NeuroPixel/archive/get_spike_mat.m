function [spike_mat,rate_mat,rate_mat_downsampled] = get_spike_mat(sp,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
defaultKernel = reshape(gausswin(401),1,[]);

addParameter(p,'kernel',defaultKernel);

parse(p,varargin{:});

kernel=p.Results.kernel;

n_units = max(sp.clu)+1;

max_t=ceil(max(sp.st)*1000);

spike_mat=zeros(n_units,max_t);
sp_times=round(sp.st*1000);
sp_times(sp_times==0)=1;
IND=sub2ind(size(spike_mat),sp.clu+1,sp_times);
spike_mat(IND)=1;

if nargout >= 2
    rate_mat = conv2(spike_mat,kernel,'same');
end

if nargout==3
    rate_mat_downsampled = rate_mat(:,1:20:end);
end

    

end

