good_cells = ks.clusters_good;
h=figure('Position',[680   363   384   615]);


[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
spike_depth = nan(numel(good_cells),1);
for k = 1:numel(good_cells)
    spike_depth(k) = median(spikeDepths(sp.clu==good_cells(k)));
end

% sort cells_to_plot by spike_depth (descending)
[spike_depth,sort_idx] = sort(spike_depth,'descend');
good_cells = good_cells(sort_idx);
similarity_matrix =struct();
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
    if nnz(idxPre)>1 && nnz(idxPost)>1
        snippetSetPre = ks.getWaveformsFromRawData('cluster_ids', clusterID,'num_waveforms', Inf, 'best_n_channels', 20, 'car', true,'spike_idx',idxPre);
        meanwfPre=mean(snippetSetPre.data,3);
        
        snippetSetPost = ks.getWaveformsFromRawData('cluster_ids', clusterID,'num_waveforms', Inf, 'best_n_channels', 20, 'car', true,'spike_idx',idxPost);
        meanwfPost=mean(snippetSetPost.data,3);
        
        %%
        figure(h)
        subplot(3,2,[1 3 5])
        [a,b]=discretize(sp.st(sp.clu==clusterID),post);
        scatter(posx(a),trial(a),2,'k')
        hold on
        [a,b]=discretize(sp.st(vr_pre),post);
        scatter(posx(a),trial(a),2,'b')
        
        [a,b]=discretize(sp.st(vr_post),post);
        scatter(posx(a),trial(a),2,'r')
        axis tight
        
        title(sprintf('cluid = %i',good_cells(k)))
        tvec=snippetSetPost.window(1):snippetSetPost.window(2);
        tvec=double(tvec)/30;
        for ii = 1:3
            subplot(3,2,ii*2)
            plot(tvec,meanwfPre(ii,:),'b')
            hold on
            plot(tvec,meanwfPost(ii,:),'r')
            axis tight
            grid on
            xlabel('time [ms]')
        end
    
    catData = cat(2,squeeze(snippetSetPre.data(1,30:80,:)),squeeze(snippetSetPost.data(1,30:80,:)));
    ff=corr(double(catData));
    similarity_matrix(k).map=ff;
    similarity_matrix(k).n_pre = snippetSetPre.nSnippets;
    similarity_matrix(k).n_post = snippetSetPost.nSnippets;
    saveas(h,fullfile(image_save_dir,sprintf('%d.png',k)),'png');
    clf
    end
end
%%

% figure
% imagesc(ff)
%st_post = ks.spike_times(idxPost);