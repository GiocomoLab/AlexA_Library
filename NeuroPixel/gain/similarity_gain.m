%contrast = load('Z:\giocomo\attialex\NP_DATA\npJ4_0512_contrast_1.mat');

[maps,maps_smoothed] = get_spatial_map(dataset);
%%
% current_map = squeeze(maps_smoothed(94,:,:));
% good_cells = dataset.sp.cids(dataset.sp.cgs==2);
% [W,H]=nnmf(current_map(:,1:234),3);
% raster_map(contrast,good_cells(94));
% axis tight
% 
% figure; 
% subplot(2,1,1);plot(W) %maps
% title('spatial component')
% subplot(2,1,2)
% plot(H') %temporal components
% title('temporal component')
% yyaxis right
% plot(dataset.trial_contrast/100)
bl_maps= maps_smoothed(:,:,dataset.trial_contrast==100 & dataset.trial_gain==1);

mean_bl_maps = nanmean(bl_maps,3);
%% using correlation as similarity
tc_1=mean(bl_maps(:,2:98,1:floor(size(bl_maps,3)/2)),3);
tc_2=mean(bl_maps(:,2:98,floor(size(bl_maps,3)/2)+1:end),3);
tmp = corr(tc_1',tc_2');
stability=diag(tmp);
[s,sidx]=sort(stability,'descend');

similarity=zeros(size(maps,3),size(maps,1));
for iC=1:size(maps,1)
    aa=corr(squeeze(maps_smoothed(iC,:,:)),mean_bl_maps(iC,:)');
    similarity(:,iC)=aa;
end
similarity(isnan(similarity))=0;
figure('Name',Files(iF).name) 
hold on
plot(nanmean(similarity,2))

n=round(numel(sidx)*.2);
plot(nanmean(similarity(:,sidx(1:n)),2))
ylabel('correlation with baseline map')
yyaxis right
plot(dataset.trial_gain,'k')
ylabel('gain')
legend({'all neurons','stable neurons','gain'});

figure('Name',Files(iF).name)  
g_values = unique(dataset.trial_gain);
avgSimilarity=ones(numel(sidx)*numel(g_values),2);
extract_win = [-4:7];
for iC=1:numel(g_values)
    onsets = strfind(dataset.trial_gain' == g_values(iC),[0 1])+1;
    snps = extract_snps(similarity',onsets,'win',[extract_win(1) extract_win(end)]);
    avg = mean(snps,3);
    avg_all= mean(avg);
    avg_stab = mean(avg(sidx(1:n),:));
    start = (iC-1)*numel(sidx)+1;
    stop = iC*numel(sidx);
    avgSimilarity(start:stop,1)=mean(avg,2);
    avgSimilarity(start:stop,2)=g_values(iC)*ones(numel(sidx),1);
    subplot(1,2,1);hold on
        plot(extract_win+1,avg_all)
        ylim([0 1])
        xlim([extract_win(1) extract_win(end)]+1)
        patch([ .9 4.1  4.1 .9],[0 0 1 1],[245/255, 224/255, 66/255],'FaceAlpha',.15,'EdgeColor','None')
    subplot(1,2,2); hold on
    plot(extract_win+1,avg_stab)
            patch([ .9 4.1  4.1 .9],[0 0 1 1],[245/255, 224/255, 66/255],'FaceAlpha',.15,'EdgeColor','None')

    ylim([0 1])
    xlim([extract_win(1) extract_win(end)]+1)
    
    %subplot(1,3,3)
    
    
    
end
%%
avgMaps={};
avgMapSimilarity=zeros(size(maps,1),numel(g_values));
for iC=1:numel(g_values)
    tmp=maps_smoothed(:,:,dataset.trial_gain==g_values(iC));
    tmpMap = nanmean(tmp,3);
    avgMaps{iC}=tmpMap;
    sim = corr(tmpMap',mean_bl_maps');
    sim=diag(sim);
    avgMapSimilarity(:,iC)=sim;
    
end

