files = dir('Z:\giocomo\attialex\NP_DATA\CL_OL_Towers\*_mismatch*.mat');
MERGED = struct();
for iF=1%1:numel(files)
clearvars -except MERGED iF files
mismatch=load(fullfile(files(iF).folder,files(iF).name));
playback=load(fullfile(files(iF).folder,strrep(files(iF).name,'mismatch','playback')));
collect_data_ol_cl
MERGED(iF).difference_by_location = difference_by_location_playback;
MERGED(iF).spatial_map_closedloop = spatialMap_mm;
MERGED(iF).spatial_map_openloop = spatialMap_pb;
MERGED(iF).corrs_playback = corrs_playback;

end

%%
SIMILAR_SCORE = [];
DIFFERENT_SCORE = [];
STABILITY=[];
for iF=1:5
n_cells = size(MERGED(iF).spatial_map_closedloop,1);
n_baseline = size(MERGED(iF).spatial_map_closedloop,3);
n_ol = size(MERGED(iF).spatial_map_openloop,3);
spatialMap_cl=MERGED(iF).spatial_map_closedloop;
spatialMap_ol=MERGED(iF).spatial_map_openloop;


tc_1=mean(spatialMap_cl(:,2:98,1:floor(n_baseline/2)),3);
tc_2=mean(spatialMap_cl(:,2:98,floor(n_baseline/2)+1:end),3);
tc_bl = mean(spatialMap_cl(:,2:98,:),3);

[a,b]=sort(MERGED(iF).corrs_playback,'descend');
frac = .2;
n=floor(frac*n_ol);
similar_idx = b(1:n);
diff_idx = b(end-n+1:end);
t_c_similar = mean(spatialMap_ol(:,:,similar_idx),3);
t_c_different = mean(spatialMap_ol(:,:,diff_idx),3);


tmp = corr(tc_1',tc_2');
stability=diag(tmp);
[~,stability_sort]=sort(stability);

score_similar = zeros(1,n_cells);
score_different = zeros(1,n_cells);
for ii = 1:n_cells
    score_similar(ii)=corr(tc_bl(ii,:)',t_c_similar(ii,2:98)');
    score_different(ii)=corr(tc_bl(ii,:)',t_c_different(ii,2:98)');
end

STABILITY=cat(1,STABILITY,stability);
SIMILAR_SCORE = cat(2,SIMILAR_SCORE,score_similar);
DIFFERENT_SCORE = cat(2,DIFFERENT_SCORE,score_different);
% figure
% subplot(1,2,1)
% scatter(score_similar,score_different,4,stability)
% 
% speed_difference_similar = nanmean(MERGED(iF).difference_by_location(similar_idx,:));
% speed_difference_different = nanmean(MERGED(iF).difference_by_location(diff_idx,:));
% axis image
% 
% xlim([-1 1])
% ylim([-1 1])
% grid on
% xlabel('similar trials')
% ylabel('different trials')
% subplot(1,2,2)
% plot(speed_difference_similar(2:end))
% hold on
% plot(speed_difference_different(2:end))

end
figure
scatter(SIMILAR_SCORE,DIFFERENT_SCORE,4,STABILITY)
axis image
xlim([-1 1])
ylim([-1 1])
grid on
xlabel('similar trials')
ylabel('different trials')
%%
linefig = figure();
avg=zeros(5,50);
for iF=1:5
    n_cells = size(MERGED(iF).spatial_map_closedloop,1);
    n_baseline = size(MERGED(iF).spatial_map_closedloop,3);
    n_ol = size(MERGED(iF).spatial_map_openloop,3);
    spatialMap_cl=MERGED(iF).spatial_map_closedloop;
    spatialMap_ol=MERGED(iF).spatial_map_openloop;
    tc_1=mean(spatialMap_cl(:,2:98,1:floor(n_baseline/2)),3);
tc_2=mean(spatialMap_cl(:,2:98,floor(n_baseline/2)+1:end),3);
tc_bl = mean(spatialMap_cl(:,2:98,:),3);



tmp = corr(tc_1',tc_2');
stability=diag(tmp);
[~,stability_sort]=sort(stability);
sim_map = zeros(n_cells,n_baseline);
[a,similarity_sort]=sort(MERGED(iF).corrs_playback,'descend');

ol_map_sorted = spatialMap_ol(:,:,similarity_sort);
for iC=1:n_cells
    m_cl=squeeze(spatialMap_cl(iC,2:98,:));
    m_ol=squeeze(ol_map_sorted(iC,2:98,:));
    tmp = corr(m_cl,m_ol);
    sim_map(iC,:)=diag(tmp);
end
idx = stability>median(stability);
%idx = stability_sort(end-50:end);
figure(linefig);
hold on
plot(linspace(0,1,n_baseline),nanmean(sim_map(idx,:)))
avg(iF,:)=interp1(linspace(0,1,n_baseline),nanmean(sim_map(idx,:)),linspace(0,1,50));
% figure
% [~,sidx]=sort(nanmean(sim_map(idx,:)));
% imagesc(MERGED(iF).difference_by_location(sidx,2:end))
end
plot(linspace(0,1,50),mean(avg),'k','LineWidth',2)
ylabel('Similarity to Closed Loop')
xlabel('Trials sorted')
%%
%%
figure
for iF=1:5
    n_cells = size(MERGED(iF).spatial_map_closedloop,1);
    n_baseline = size(MERGED(iF).spatial_map_closedloop,3);
    n_ol = size(MERGED(iF).spatial_map_openloop,3);
    spatialMap_cl=MERGED(iF).spatial_map_closedloop;
    spatialMap_ol=MERGED(iF).spatial_map_openloop;
    tc_1=mean(spatialMap_cl(:,2:98,1:floor(n_baseline/2)),3);
tc_2=mean(spatialMap_cl(:,2:98,floor(n_baseline/2)+1:end),3);
tc_bl = mean(spatialMap_cl(:,2:98,:),3);



tmp = corr(tc_1',tc_2');
stability=diag(tmp);
[~,stability_sort]=sort(stability);
sim_map = zeros(n_cells,100);
[a,similarity_sort]=sort(MERGED(iF).corrs_playback,'descend');


for iC=1:n_cells
    m_cl=squeeze(spatialMap_cl(iC,1:100,:));
    m_ol=squeeze(spatialMap_ol(iC,1:100,:));
    tmp = corr(m_cl',m_ol');
    sim_map(iC,:)=diag(tmp);
end
idx = stability>quantile(stability,.50);
%idx = stability_sort(end-50:end);
hold on
plot(linspace(0,400,size(sim_map,2)),nanmean(sim_map(idx,:)))
end
ylabel('Similarity to Closed Loop')
xlabel('Position')

%%
figure
hold on

for iF=1:5
    tmp = nanmean(MERGED(iF).difference_by_location(:,2:end));
    plot(linspace(0,400,numel(tmp)),tmp)
end
xlabel('Position')
ylabel('Vis/Run speed difference')