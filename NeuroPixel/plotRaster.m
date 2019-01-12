%%ba = allBA{e}; bins = allBins{e};

binned_array=squeeze(spike_mat(274,:,:));
bins=win(1):0.001:win(2);
    
    % set the new rasters
    [tr,b] = find(binned_array);
    [rasterX,yy] = rasterize(bins(b));
    rasterY = yy*1+reshape(repmat(tr',3,1),1,length(tr)*3)-0.5; % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
figure;plot(rasterX,rasterY)
hold on
%plot(rasterX,rasterY+8,'r')