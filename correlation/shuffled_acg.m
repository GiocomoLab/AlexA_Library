% function [mean_shuffles, quantiles] = shuffled_acg(spike_id,binedges,gauss_filter,time_per_bin,n_it)
% shuffles=randi([-length(spike_id),length(spike_id)],n_it,1);
% tmp_acg=zeros(n_it,4001);
%         for num_it = 1:n_it
%             shuffle_idx = shuffles(num_it);
%             spike_id_shuffled = circshift(spike_id,shuffle_idx);
%             
%         spike_t = sp.st(spike_id_shuffled);
%         [~,~,spike_idx] = histcounts(spike_t,post);
%         spike_distance = distance(spike_idx);
%         
% % smoothSigma = 3;%params.SmoothSigmaFR/params.SpatialBin;
% % smoothWindow = floor(smoothSigma*5/2)*2+1;
% % gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
% 
% 
% 
%         
%         [acg,spacing,pxx,f] = calculate_spatial_acg(spike_distance,binedges,gauss_filter,time_per_bin);
%          tmp_acg(num_it,:)=acg;
%         end
%         mean_shuffles=mean(tmp_acg,1);
%         quantiles = quantile(tmp_acg,[.95 .05]);
%         
% 
% end
% 
function [mean_shuffles, quantiles] = shuffled_acg(sp,post,distance,spike_id,binedges,gauss_filter,time_per_bin,n_it)
spike_t = sp.st(spike_id);
if numel(spike_t)<1
    mean_shuffles = NaN(1,4001);
    quantiles = NaN(2,4001);
    warning('cluster has no spikes')
    return
end
max_t = spike_t(end);
abs_min = min(sp.st);
shuffles=max_t*rand(n_it,1);

tmp_acg=zeros(n_it,4001);
        for num_it = 1:n_it
            shuffle_idx = shuffles(num_it);
            spike_t_shuffled = mod(spike_t+shuffle_idx,max_t)+abs_min;
            
        
        [~,~,spike_idx] = histcounts(spike_t_shuffled,post);
        spike_distance = distance(spike_idx);
        
% smoothSigma = 3;%params.SmoothSigmaFR/params.SpatialBin;
% smoothWindow = floor(smoothSigma*5/2)*2+1;
% gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);



        
        [acg,spacing,pxx,f] = calculate_spatial_acg(spike_distance,binedges,gauss_filter,time_per_bin);
         tmp_acg(num_it,:)=acg;
        end
        mean_shuffles=mean(tmp_acg,1);
        quantiles = quantile(tmp_acg,[.95 .05]);
        

end

