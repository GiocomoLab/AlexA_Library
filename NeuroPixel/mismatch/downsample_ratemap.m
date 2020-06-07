function rate_mat=downsample_ratemap(spike_mat,stepsize)

rate_mat = zeros(size(spike_mat,1),size(spike_mat,2),floor(8001/stepsize));
cntr = 0;
for iS = 1:stepsize:(8001-stepsize)
    cntr = cntr+1;
    idx = iS:iS+stepsize;
    n_spikes =  sum(spike_mat(:,:,idx),3);
    rate_mat(:,:,cntr)=n_spikes;
end