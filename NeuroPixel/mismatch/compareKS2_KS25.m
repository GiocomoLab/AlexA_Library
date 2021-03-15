matfiles = dir('F:\Alex\matfiles_new\AA_2009*MM_*');
path25 = 'F:\Alex\new_2';
%%
MM=[];
MM25=[];
RUN=[];
CHANNEL = [];
FR = [];
FR25=[];
opt = load_mismatch_opt;
for iF=1:numel(matfiles)
file25 = fullfile(path25,strrep(matfiles(iF).name,'.mat','_KS25.mat'));
if ~isfile(file25)
    continue
end
    data = load(fullfile(matfiles(1).folder,matfiles(iF).name));
data25 = load(file25);

good_cells = data.sp.ks_cluster.cluster_id(startsWith(data.sp.ks_cluster.KSLabel,'good'));
good_cells25 = data25.sp.ks_cluster.cluster_id(startsWith(data25.sp.ks_cluster.KSLabel,'good'));

[~,count_vec]=extractMM(data,good_cells,opt);
duration = max(data.post)-min(data.post);
fr = nan(size(good_cells));
for iC=1:numel(good_cells)
    n = nnz(data.sp.clu==good_cells(iC));
    fr(iC)=n/duration;
end
count_vecN = count_vec-mean(count_vec(:,opt.time_vecs<-.1 & opt.time_vecs>-.6),2);
%count_vecN = smoothdata(count_vecN,2,'gaussian',5);
count_vecN = count_vecN./fr;

chan_number =data.sp.waveform_metrics.peak_channel(ismember(data.sp.waveform_metrics.cluster_id,good_cells));
[~,sid]=sort(chan_number,'descend');

subplot(2,2,1)
imagesc(count_vecN(sid,:),[-2 2])
subplot(2,2,3)
hold on
plot(opt.time_vecs,nanmean(count_vecN(fr>1,:)),'b')


[~,count_vec25]=extractMM(data25,good_cells25,opt);
duration = max(data.post)-min(data.post);
fr25 = nan(size(good_cells25));
for iC=1:numel(good_cells25)
    n = nnz(data25.sp.clu==good_cells25(iC));
    fr25(iC)=n/duration;
end
count_vecN = count_vec25-mean(count_vec25(:,opt.time_vecs<-.1 & opt.time_vecs>-.6),2);
%count_vecN = smoothdata(count_vecN,2,'gaussian',5);
count_vecN = count_vecN./fr25;

chan_number =data25.sp.waveform_metrics.peak_channel(ismember(data25.sp.waveform_metrics.cluster_id,good_cells25));
[~,sid]=sort(chan_number,'descend');

subplot(2,2,2)
imagesc(count_vecN(sid,:),[-2 2])
subplot(2,2,3)
hold on
plot(opt.time_vecs,nanmean(count_vecN(fr25>1,:)),'r')



pause
clf
MM=cat(1,MM,count_vec);
MM25=cat(1,MM25,count_vec25);
CHANNEL = cat(1,CHANNEL,chan_number);
FR=cat(1,FR,fr);
FR25= cat(1,FR25,fr25);
end
%%
valid_idx = FR>1;
[~,sid]=sort(CHANNEL(valid_idx),'descend');
MM_smooth = smoothdata(MM(valid_idx,:),2,'gaussian',5);
RUN_smooth = smoothdata(RUN(valid_idx,:),2,'gaussian',5);
%count_vecN=MM_smooth./max(MM_smooth,[],2);
%count_vec_runN = RUN_smooth./max(RUN_smooth,[],2);
MM_smooth = MM_smooth-nanmean(MM_smooth(:,opt.time_vecs<-.1 & opt.time_vecs>-.6),2);
count_vecN = MM_smooth./FR(valid_idx);
RUN_smooth = RUN_smooth-nanmean(RUN_smooth(:,opt.time_vecs<-.1 & opt.time_vecs>-.6),2);
count_vec_runN = RUN_smooth./FR(valid_idx);
figure('Color','w')
subplot(3,2,[1 3])
imagesc(opt.time_vecs,1:nnz(valid_idx),count_vecN(sid,:),[-2 2]);
ylabel('ventral -> dorsal')
xlim([-1 3])
subplot(3,2,[2 4])
imagesc(opt.time_vecs,1:nnz(valid_idx),count_vec_runN(sid,:),[-2 2]);
xlim([-1 3])

subplot(3,2,5)
%plot(opt.time_vecs,nanmean(count_vecN));
boundedline(opt.time_vecs,nanmean(count_vecN),nanstd(count_vecN)/sqrt(nnz(valid_idx)))
xlim([-1 3])

title('MM')
subplot(3,2,6)
boundedline(opt.time_vecs,nanmean(count_vec_runN),nanstd(count_vec_runN)/sqrt(nnz(valid_idx)))
title('Run')
xlim([-1 3])


%%
figure
valid_idx = mean(RUN,2)>1 & mean(MM,2)>1;
[~,sid]=sort(CHANNEL(valid_idx),'descend');
MM_smooth = smoothdata(MM(valid_idx,:),2,'gaussian',5);
tmp = MM_smooth-nanmean(MM_smooth(:,opt.time_vecs<-.1 & opt.time_vecs>-.6),2);
imagesc(opt.time_vecs,1:nnz(valid_idx),tmp(sid,:),[-2 2])

%%
    valid_idx = FR>1;
[~,sid]=sort(CHANNEL(valid_idx),'descend');
MM_smooth = smoothdata(MM(valid_idx,:),2,'gaussian',5);
count_vecN = MM_smooth./FR(valid_idx);
figure
hold on
for ii=1:4
    sta = round((ii-1)*numel(sidx)/4)+1;
    sto = round(ii*numel(sidx)/4);
    take_idx = sta:sto;
    idx = sid(take_idx);
    boundedline(opt.time_vecs,nanmean(count_vecN(idx,:)),nanstd(count_vecN(idx,:))/sqrt(nnz(idx)))
xlim([-1 3])
end
%%
   valid_idx = FR>1;
[~,sid]=sort(CHANNEL(valid_idx),'descend');
MM_smooth = smoothdata(MM(valid_idx,:),2,'gaussian',5);
count_vecN = MM_smooth./FR(valid_idx);
figure
hold on
col = jet(4);
for ii=1:4
    sta = round((ii-1)*numel(sid)/4)+1;
    sto = round(ii*numel(sid)/4);
    take_idx = sta:sto;
    idx = sid(take_idx);
    boundedline(opt.time_vecs,nanmean(count_vecN(idx,:)),nanstd(count_vecN(idx,:))/sqrt(nnz(idx)),'cmap',col(ii,:),'alpha')
xlim([-1 3])
end