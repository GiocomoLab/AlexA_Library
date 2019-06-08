contrast = load('Z:\giocomo\attialex\NP_DATA\npJ4_0512_contrast_1.mat');

[maps,maps_smoothed] = get_spatial_map(contrast);
%%
current_map = squeeze(maps_smoothed(94,:,:));

[W,H]=nnmf(current_map(:,1:234),3);

figure; plot(W) %maps
plot(H') %temporal components

bl_maps= maps_smoothed(:,:,contrast.trial_contrast==100 & contrast.trial_gain==1);

mean_bl_maps = nanmean(bl_maps,3);
%% using correlation as similarity
tc_1=mean(bl_maps(:,2:98,1:floor(size(bl_maps,3)/2)),3);
tc_2=mean(bl_maps(:,2:98,floor(size(bl_maps,3)/2)+1:end),3);
tmp = corr(tc_1',tc_2');
stability=diag(tmp);
[s,sidx]=sort(stability,'descend');

similarity=zeros(size(maps,3),size(maps,1));
for iC=1:205
    aa=corr(squeeze(maps_smoothed(iC,:,:)),mean_bl_maps(iC,:)');
    similarity(:,iC)=aa;
end
similarity(isnan(similarity))=0;
figure
hold on
plot(nanmean(similarity,2))

n=round(numel(sidx)*.2);
plot(nanmean(similarity(:,sidx(1:n)),2))

ylabel('correlation with baseline map')
yyaxis right
plot(contrast.trial_contrast)
ylabel('contrast')

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
% similarity(isnan(similarity))=0;
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
