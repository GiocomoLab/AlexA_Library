%contrast = load('Z:\giocomo\attialex\NP_DATA\npJ4_0512_contrast_1.mat');

[maps,maps_smoothed] = get_spatial_map(contrast);
%%
% current_map = squeeze(maps_smoothed(94,:,:));
% good_cells = contrast.sp.cids(contrast.sp.cgs==2);
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
% plot(contrast.trial_contrast/100)
bl_maps= maps_smoothed(:,:,contrast.trial_contrast==100 & contrast.trial_gain==1);

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
plot(contrast.trial_contrast,'k')
ylabel('contrast')
legend({'all neurons','stable neurons','contrast'});

figure('Name',Files(iF).name)  
c_values = unique(contrast.trial_contrast);
avgSimilarity=ones(numel(sidx)*numel(c_values),2);
for iC=1:numel(c_values)
    onsets = strfind(contrast.trial_contrast' == c_values(iC),[0 1])+1;
    snps = extract_snps(similarity',onsets,'win',[0 9]);
    avg = mean(snps,3);
    avg_all= mean(avg);
    avg_stab = mean(avg(sidx(1:n),:));
    start = (iC-1)*numel(sidx)+1;
    stop = iC*numel(sidx);
    avgSimilarity(start:stop,1)=mean(avg,2);
    avgSimilarity(start:stop,2)=c_values(iC)*ones(numel(sidx),1);
    subplot(1,2,1);hold on
        plot(1:10,avg_all)
        ylim([0 1])
        xlim([1 10])
    subplot(1,2,2); hold on
    plot(1:10,avg_stab)
    ylim([0 1])
    xlim([1 10])
    
    %subplot(1,3,3)
    
    
    
end
%%
avgMaps={};
avgMapSimilarity=zeros(size(maps,1),numel(c_values));
for iC=1:numel(c_values)
    tmp=maps_smoothed(:,:,contrast.trial_contrast==c_values(iC) & contrast.trial_gain==1);
    tmpMap = nanmean(tmp,3);
    avgMaps{iC}=tmpMap;
    sim = corr(tmpMap',mean_bl_maps');
    sim=diag(sim);
    avgMapSimilarity(:,iC)=sim;
    
end

%% calculate gain reponse
if     numel(unique(contrast.trial_gain))>1

        gain_onset = strfind(contrast.trial_gain'==.5,[0 1])+1;
        gain_onset = gain_onset(1);
        
        pre_trials = gain_onset-5:gain_onset-1;
        gain_trials = gain_onset:gain_onset+4;
        post_trials = gain_onset+5:gain_onset+9;
        
        pre_map = mean(maps_smoothed(:,:,pre_trials),3);
        post_map = mean(maps_smoothed(:,:,post_trials),3);
        gain_map = mean(maps_smoothed(:,:,gain_trials),3);
        
        gain_stability = diag(corr(pre_map',post_map'));
        gain_response = diag(corr(pre_map',gain_map'));
        xcorr_pre_post=zeros(2,size(pre_map,1));
        xcorr_pre_gain=zeros(2,size(pre_map,1));

        for ii=1:size(pre_map,1)
            pre_tmp=pre_map(ii,:);
            pre_tmp = pre_tmp-mean(pre_tmp);
            post_tmp=post_map(ii,:);
            post_tmp = post_tmp-mean(pre_tmp);
            gain_tmp = gain_map(ii,:);
            gain_tmp = gain_tmp-mean(gain_tmp);
            [prepost]=xcorr(pre_tmp,post_tmp,10,'coeff');
            [pregain]=xcorr(pre_tmp,gain_tmp,10,'coeff');
            [aa,jj]=max(prepost);
            xcorr_pre_post(:,ii)=[aa,jj];
            [aa,jj]=max(pregain);
            xcorr_pre_gain(:,ii)=[aa,jj];
        end
else
    xcorr_pre_post = [];
    xcorr_pre_gain =[];
    gain_response = [];
    gain_stability=[];
end


% figure
% 
% for iC=1:numel(c_values)
%     subplot(1,7,iC)
%     imagesc(avgMaps{iC}(sidx,:))
% end
% figure
% imagesc(avgMapSimilarity(sidx,1:6))
% figure('Name',Files(iF).name)  
% subplot(4,1,1)
% plot(contrast.trial_contrast);
% axis tight
% ylim([0 100])
% ylabel('contrast')
% subplot(4,1,[2:4])
% imagesc(similarity(:,sidx)')
% title('similarity to mean baseline')
% ylabel('unit')

%%
% figure
% 
% for iC=1:length(sidx)
%     cc=sidx(iC);
% current_map = squeeze(maps_smoothed(cc,:,:));
% ax1=subplot(2,2,[1 3]);
% [W,H]=nnmf(current_map(13:87,1:234),3);
% spike_id = contrast.sp.clu==good_cells(cc);
%     spike_t = contrast.sp.st(spike_id);
% [~,~,spike_idx] = histcounts(spike_t,contrast.post);
% current_trial = contrast.trial(spike_idx);
% cols = zeros(length(current_trial),3);
% t1=find(H(1,:)>H(2,:) & H(1,:)> H(3,:));
% idx1 = ismember(current_trial,t1);
% cols(idx1,:)=repmat([1 0 0],nnz(idx1),1);
% 
% %find dominant component of 0 contrast trials
% [aa,bb]=max(H);
% [dominant_component,freq] = mode(bb(find(contrast.trial_contrast(1:234)==0)));
% cols = zeros(length(current_trial),3);
% t1=find(bb~=dominant_component);
% idx1 = ismember(current_trial,t1);
% cols(idx1,:)=repmat([1 0 0],nnz(idx1),1);
% raster_map(contrast,good_cells(cc),ax1,cols);
% axis tight
% ylim([1 234])
% 
% 
% subplot(2,2,2);plot(W) %maps
% title('spatial component')
% xlabel(sprintf('Dominant: %d, frac %.2f',dominant_component,freq/nnz(contrast.trial_contrast==0)));
% subplot(2,2,4)
% plot(H') %temporal components
% title('temporal component')
% yyaxis right
% plot(contrast.trial_contrast/100)
% pause
% clf
% end
%% quantification NNMF
% compMap = zeros(size(maps_smoothed,1),nnz(contrast.trial_gain==1));
% zeroMap = compMap;
% for iC=1:size(maps_smoothed,1)
%     current_map = squeeze(maps_smoothed(iC,:,:));
%     [W,H]=nnmf(current_map(:,contrast.trial_gain==1),3);
%     [aa,bb]=max(H);
%     [dominant_component] = mode(bb(contrast.trial_contrast(contrast.trial_gain ==1)==min(contrast.trial_contrast)));
%     compMap(iC,:)=bb;
%     zeroMap(iC,:) = bb==dominant_component;
% end
% figure('Name',Files(iF).name)  
% subplot(2,1,1)
% imagesc(compMap(sidx,:));
% title('component map')
% 
% subplot(2,1,2)
% imagesc(zeroMap(sidx,:))
% title('min contrast map dominance')

%% same with pca

% 
% for iC=1:length(sidx)
%     cc=sidx(iC);
%     spike_id = contrast.sp.clu==good_cells(cc);
%     spike_t = contrast.sp.st(spike_id);
% [~,~,spike_idx] = histcounts(spike_t,contrast.post);
% current_trial = contrast.trial(spike_idx);
% 
% 
% current_map = squeeze(maps_smoothed(cc,:,:));
% ax1=subplot(2,2,[1 3]);
% [coeff,score,lat,tsquared,explained,mu]=pca(current_map(13:87,1:234)');
% 
% [aa,bb]=max(abs(score(:,1:3))');
% [dominant_component,freq] = mode(bb(find(contrast.trial_contrast==0)));
% cols = zeros(length(current_trial),3);
% t1=find(bb~=dominant_component);
% idx1 = ismember(current_trial,t1);
% cols(idx1,:)=repmat([1 0 0],nnz(idx1),1);
% 
% raster_map(contrast,good_cells(cc),ax1,cols);
% axis tight
% 
% 
% subplot(2,2,2);plot(coeff(:,1:3)) %maps
% title('spatial component')
% subplot(2,2,4)
% plot(score(:,1:3)) %temporal components
% title('temporal component')
% yyaxis right
% plot(contrast.trial_contrast/100)
% pause
% clf
% end

%% using cosine as similarity
% tc_1=mean(bl_maps(:,2:98,1:floor(size(bl_maps,3)/2)),3);
% tc_2=mean(bl_maps(:,2:98,floor(size(bl_maps,3)/2)+1:end),3);
% tmp = tc_1*tc_2';
% n1=diag(((tc_1*tc_1')));
% n2=diag(((tc_2*tc_2')));
% 
% ff=diag(tmp);
% 
% sim=ff./(sqrt(n1.*n2));
% 
% [s,sidx]=sort(sim,'ascend');
% 
% similarity=zeros(size(maps,3),size(maps,1));
% for iC=1:205
%     aa=corr(squeeze(maps_smoothed(iC,:,:)),mean_bl_maps(iC,:)');
%     similarity(:,iC)=aa;
% end
% % similarity(isnan(similarity))=0;
% figure
% subplot(2,1,1)
% hold on
% plot(nanmean(similarity,2))
% 
% n=round(numel(sidx)*.2);
% plot(nanmean(similarity(:,sidx(1:n)),2))
% 
% ylabel('correlation with baseline map')
% yyaxis right
% plot(contrast.trial_contrast)
% ylabel('contrast')
% 
% 
