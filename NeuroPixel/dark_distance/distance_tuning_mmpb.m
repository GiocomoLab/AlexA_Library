opt = load_default_opt;
opt.dark = true;
opt.SpatialBin = 2; % in cm
opt.smoothSigma_dist = 4; % in cm
opt.max_lag = 800; % in cm (for spatial autocorrelation)
opt.num_shuf = 300;
    mm_opt = load_mismatch_opt;

save_path = '/Volumes/T7/attialex/distance_tuning_mm';
%%
mf = dir('/Volumes/T7/attialex/NP_DATA_corrected/np*_mismatch*.mat')

mf_dark = cell(numel(mf),1);
for ii=1:numel(mf)
    tmp = strsplit(mf(ii).name,'_');
    name_first_part = [tmp{1} '_' tmp{2} '_dark_*'];
    potentials = dir(fullfile(mf(ii).folder,name_first_part));
    if numel(potentials)>=1
        for jj=1:numel(potentials)
            if ~ contains(potentials(jj).name,'reward');
                mf_dark{ii} = potentials(jj).name;
            end
        end
    end
end

%%

for iF=1:numel(mf)
    if isempty(mf_dark{iF})
        continue
    end
    data_mm = load(fullfile(mf(iF).folder,mf(iF).name));
    data_dark = load(fullfile(mf(iF).folder,mf_dark{iF}));
    
    good_cells = data_mm.sp.cids(data_mm.sp.cgs==2);
    %good_cells=good_cells(1:5);
    
    %dist_tuning_mm = calc_distance_tuning(data_mm,good_cells,opt);
    %dist_tuning_dark = calc_distance_tuning(data_dark,good_cells,opt);
    n_spikes = zeros(1,numel(good_cells));
for iC=1:numel(good_cells)
n_spikes(iC) = nnz(data_mm.sp.clu==good_cells(iC));
end

n_spikes = n_spikes/max(data_mm.post);
    [cv,cvr]=preview_MM(data_mm,true);
%     count_vecN=cv./nanmedian(cv,2);
%     count_vec_runN = cvr./nanmedian(cvr,2);
%     figure('Name',mf(iF).name)
%     [~,sid]=sort(dist_tuning_dark.peak_loc_all);
%     subplot(1,4,1)
%     imagesc(0:opt.SpatialBin:opt.max_lag,1:numel(sid),dist_tuning_mm.xcorrs(sid,:),[0 0.4])
%     subplot(1,4,2)
%     imagesc(0:opt.SpatialBin:opt.max_lag,1:numel(sid),dist_tuning_dark.xcorrs(sid,:),[0 0.4])
%     subplot(1,4,3)
%     imagesc(mm_opt.time_vecs,1:numel(sid),count_vecN(sid,:),[0 3])
%     subplot(1,4,4)
%     imagesc(mm_opt.time_vecs,1:numel(sid),count_vec_runN(sid,:),[0 3])
%     drawnow
%     colormap summer
    fn = fullfile(save_path,mf(iF).name);
    matfile_out = matfile(fn,'Writable',true);
    matfile_out.count_vec_run = cvr;
    matfile_out.count_vec = cv;
    matfile_out.cluster_group = data_dark.sp.cgs;
    matfile_out.fr = n_spikes;
    %matfile_out.anatomy = data_dark.anatomy;
    %matfile_out.dist_tuning_mm = dist_tuning_mm;
    %matfile_out.dist_tuning_dark = dist_tuning_dark;
    
end
%% load and plot processed data
mf = dir(fullfile(save_path,'*.mat'));
prom = [];
mm_resp_all = [];
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
    %[~,sid]= sort(mm_resp);
    %mm=zscore(mm,[],2);
    subplot(1,3,1)
    imagesc(0:opt.SpatialBin:opt.max_lag,1:numel(sid),data.dist_tuning_mm.xcorrs(sid,:),[0 0.4])
    subplot(1,3,2)
    imagesc(0:opt.SpatialBin:opt.max_lag,1:numel(sid),data.dist_tuning_dark.xcorrs(sid,:),[0 0.4])
    subplot(1,3,3)
    imagesc(mm_opt.time_vecs,1:numel(sid),mmr(sid,:),[-1 1])
    %subplot(1,4,4)
    %imagesc(mm_opt.time_vecs,1:numel(sid),data.count_vec_run(sid,:),[0 3])
    prom=cat(1,prom,[data.dist_tuning_dark.peak_prom_all,data.dist_tuning_dark.pvals]);
    mm_resp_all = cat(1,mm_resp_all,mmr);
end
%%
DT = prom(:,2)<0.01 & prom(:,1)>0.1;
figure
plot(mm_opt.time_vecs,mean(mm_resp_all(DT,:)));
hold on
plot(mm_opt.time_vecs,mean(mm_resp_all(~DT,:)));
legend({'DT','non tuned'})
ylabel('normalized firing rate change')
%% load and plot processed data
mf = dir(fullfile(save_path,'*.mat'));
prom = [];
mm_resp_all = [];
for iF=1:numel(mf)
    data = load(fullfile(mf(iF).folder,mf(iF).name));

        figure('Name',mf(iF).name)
        DT=data.dist_tuning_dark.pvals<0.01 & data.dist_tuning_dark.peak_prom_all>0.1;
        %DT=~DT;
    mm = data.count_vec;
    mm = smoothdata(mm,2,'gaussian',5);
    mmr = mm-nanmean(mm(:,mm_opt.time_vecs<-0.05 & mm_opt.time_vecs>-0.55),2);
    mmr = mmr(DT,:)./data.fr(DT)';
    xcorrs = data.dist_tuning_dark.xcorrs(DT,:);
    xcorrs_mm = data.dist_tuning_mm.xcorrs(DT,:);
        %[~,sid]=sort(data.dist_tuning_dark.peak_prom_all);
        tip_distance = data.anatomy.tip_distance(data.cluster_group==2);
    [~,sid]=sort(tip_distance(DT),'descend');
    mm_resp = nanmean(mm(:,mm_opt.time_vecs>0.05 & mm_opt.time_vecs<0.6),2)-nanmean(mm(:,mm_opt.time_vecs<-0.05 & mm_opt.time_vecs>-0.55),2);
    %[~,sid]= sort(mm_resp);
    %mm=zscore(mm,[],2);
    subplot(1,3,1)
    imagesc(0:opt.SpatialBin:opt.max_lag,1:numel(sid),xcorrs_mm(sid,:),[0 0.4])
    subplot(1,3,2)
    imagesc(0:opt.SpatialBin:opt.max_lag,1:numel(sid),xcorrs(sid,:),[0 0.4])
    subplot(1,3,3)
    imagesc(mm_opt.time_vecs,1:numel(sid),smoothdata(mmr(sid,:),2,'gaussian',5),[-1 1])
    
    prom=cat(1,prom,[data.dist_tuning_dark.peak_prom_all,data.dist_tuning_dark.pvals]);
    mm_resp_all = cat(1,mm_resp_all,mmr);
end

%%

%%
data = load('/Volumes/T7/attialex/distance_tuning_mm/npI1_0418_mismatch_notowers_2.mat');
figure
scatter(data.dist_tuning_dark.peak_loc_all,data.dist_tuning_mm.peak_loc_all)
hold on
refline(1,0)
xlabel('dark')
ylabel('mm optic flow')
