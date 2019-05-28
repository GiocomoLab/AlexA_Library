good_cells = ks.clusters_good;




% sort cells_to_plot by spike_depth (descending)
similarity_map={};
valid_idx=false(size(good_cells));
for k=1:numel(good_cells)
    clusterID = good_cells(k);
    
    post_offset = post+sp.vr_session_offset;
    %trial range pre
    start_pre = post_offset(find((trial-trial_pre(1))==0,1));
    stop_pre = post_offset(find((trial-trial_pre(end))==0,1));
    
    %trial range post
    start_post = post_offset(find((trial-trial_post(1))==0,1));
    stop_post = post_offset(find((trial-trial_post(end))==0,1));
    
    vr_pre = sp.st>(start_pre-sp.vr_session_offset) & sp.st<(stop_pre-sp.vr_session_offset) & sp.clu == clusterID;
    vr_post = sp.st>(start_post-sp.vr_session_offset) & sp.st<(stop_post-sp.vr_session_offset) & sp.clu == clusterID;
    
    start_pre = start_pre*30000;
    stop_pre = stop_pre*30000;
    idxPre = ks.spike_times > start_pre & ks.spike_times <stop_pre & ks.spike_clusters == clusterID;
    
    start_post = start_post*30000;
    stop_post = stop_post*30000;
    idxPost = ks.spike_times > start_post & ks.spike_times <stop_post & ks.spike_clusters == clusterID;
    %time range post
    if nnz(idxPre)>=100 && nnz(idxPost)>=100
        valid_idx(k)=true;
        snippetSetPre = ks.getWaveformsFromRawData('cluster_ids', clusterID,'num_waveforms', 100, 'best_n_channels', 20, 'car', true,'spike_idx',idxPre);
        meanwfPre=mean(snippetSetPre.data,3);
        
        snippetSetPost = ks.getWaveformsFromRawData('cluster_ids', clusterID,'num_waveforms', 100, 'best_n_channels', 20, 'car', true,'spike_idx',idxPost);
        meanwfPost=mean(snippetSetPost.data,3);
        
        %%
        
        
      
    
    catData = cat(2,squeeze(snippetSetPre.data(1,30:80,:)),squeeze(snippetSetPost.data(1,30:80,:)));
    ff=corr(double(catData));
    similarity_map{k}=ff;

    end
end
%%
spatial_map=nanmean(correlation_All(trial_pre(1):trial_post(length(trial_post)),trial_pre(1):trial_post(length(trial_post)),valid_idx),3);
figure
subplot(1,2,1)
imagesc(spatial_map,[0 .4])
colorbar
title('Spatial Similarity')
xlabel('Trial #')
ylabel('Trial #')
set(gca,'XTick',[1:5:40],'XTickLabel',[trial_pre(1):5:trial_post(end)])
set(gca,'YTick',[1:5:40],'YTickLabel',[trial_pre(1):5:trial_post(end)])
axis image
subplot(1,2,2)
cc=cat(3,similarity_map{:});
cc=mean(cc,3);
imagesc(cc,[0 1]);
title('Waveform similarity')
xlabel('100 spikes pre, 100 spikes post')
ylabel('spike #')
colorbar
axis image
%%


h=figure();
for k=1:numel(good_cells)
    if valid_idx(k)
        tmp_spatial = correlation_All(trial_pre(1):trial_post(length(trial_post)),trial_pre(1):trial_post(length(trial_post)),k);
        tmp_spatial = squeeze(tmp_spatial);
        subplot(1,2,1)
        imagesc(tmp_spatial,prctile(tmp_spatial(:),[25 75]))
        xlabel('Trial #')
ylabel('Trial #')
        title('Spatial Similarity')
        colorbar
        axis image
        subplot(1,2,2)
        tmp_sim = similarity_map{k};
        imagesc(tmp_sim,prctile(tmp_sim(:),[5 95]))
        title('Waveform similarity')
        xlabel('100 spikes pre, 100 spikes post')
        colorbar
        axis image
        saveas(h,fullfile(image_save_dir,sprintf('%d.png',k)),'png');
        %pause
        clf
    end

end
% figure
% imagesc(ff)
%st_post = ks.spike_times(idxPost);