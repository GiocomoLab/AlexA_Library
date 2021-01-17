opt = load_default_opt;
opt.dark = true;
opt.SpatialBin = 2; % in cm
opt.smoothSigma_dist = 4; % in cm
opt.max_lag = 800; % in cm (for spatial autocorrelation)
opt.num_shuf = 300;
save_path = '/Volumes/T7/attialex/distance_tuning_mm_mec';
%%
mf = dir('/Volumes/T7/attialex/mismatch_mec/*_mismatch*.mat');

needs to be rerun: autocorrelations look different for H3 probes -why? Also, move trigger not always detected correctly

%%

for iF=1:numel(mf)
    
    data_mm = load(fullfile(mf(iF).folder,mf(iF).name));
    
    good_cells = data_mm.sp.cids(data_mm.sp.cgs==2);
    %cgood_cells=good_cells(1:5);
    
    dist_tuning_mm = calc_distance_tuning(data_mm,good_cells,opt);
    mm_opt = load_mismatch_opt;
    n_spikes = zeros(1,numel(good_cells));
    for iC=1:numel(good_cells)
        n_spikes(iC) = nnz(data_mm.sp.clu==good_cells(iC));
    end
    
    n_spikes = n_spikes/max(data_mm.post);
    [cv,cvr]=preview_MM(data_mm,true);
    count_vecN=cv./nanmedian(cv,2);
    count_vec_runN = cvr./nanmedian(cvr,2);
    figure('Name',mf(iF).name)
    [~,sid]=sort(dist_tuning_mm.peak_loc_all);
    subplot(1,3,1)
    imagesc(0:opt.SpatialBin:opt.max_lag,1:numel(sid),dist_tuning_mm.xcorrs(sid,:),[0 0.4])
    subplot(1,3,2)
    imagesc(mm_opt.time_vecs,1:numel(sid),count_vecN(sid,:),[0 3])
    subplot(1,3,3)
    imagesc(mm_opt.time_vecs,1:numel(sid),count_vec_runN(sid,:),[0 3])
    drawnow
    colormap summer
    fn = fullfile(save_path,mf(iF).name);
    matfile_out = matfile(fn,'Writable',true);
    matfile_out.count_vec_run = cvr;
    matfile_out.count_vec = cv;
    matfile_out.cluster_group = data_mm.sp.cgs;
    matfile_out.fr = n_spikes;
    if isfield(data_mm,'anatomy')
        matfile_out.anatomy = data_mm.anatomy;
    end
    matfile_out.dist_tuning_mm = dist_tuning_mm;
    %matfile_out.dist_tuning_dark = dist_tuning_dark;
    
end
%%
mf = dir(fullfile(save_path,'*.mat'));
for iF=1:numel(mf)
    data = load(fullfile(mf(iF).folder,mf(iF).name));
    
    figure('Name',mf(iF).name)
    mm = data.count_vec;
    mm = smoothdata(mm,2,'gaussian',5);
    mmr = mm-nanmean(mm(:,mm_opt.time_vecs<-0.05 & mm_opt.time_vecs>-0.55),2);
    mmr = mmr./data.fr';
    %[~,sid]=sort(data.dist_tuning_dark.peak_prom_all);
    [~,sid]=sort(data.anatomy.tip_distance(data.cluster_group==2),'descend');
    mm_resp = nanmean(mm(:,mm_opt.time_vecs>0.05 & mm_opt.time_vecs<0.6),2)-nanmean(mm(:,mm_opt.time_vecs<-0.05 & mm_opt.time_vecs>-0.55),2);
    [~,sid]= sort(mm_resp);
    %mm=zscore(mm,[],2);
    subplot(1,4,1)
    imagesc(0:opt.SpatialBin:opt.max_lag,1:numel(sid),data.dist_tuning_mm.xcorrs(sid,:),[0 0.4])
    subplot(1,4,2)
    imagesc(0:opt.SpatialBin:opt.max_lag,1:numel(sid),data.dist_tuning_dark.xcorrs(sid,:),[0 0.4])
    subplot(1,4,3)
    imagesc(mm_opt.time_vecs,1:numel(sid),mmr(sid,:),[-1 1])
    subplot(1,4,4)
    imagesc(mm_opt.time_vecs,1:numel(sid),data.count_vec_run(sid,:),[0 3])
    
end
%%
mf = load('/Volumes/T7/attialex/distance_tuning_mm/npI1_0418_mismatch_notowers_2.mat');
figure
scatter(data.dist_tuning_dark.peak_loc_all,data.dist_tuning_mm.peak_loc_all)
hold on
refline(1,0)
xlabel('dark')
ylabel('mm optic flow')
