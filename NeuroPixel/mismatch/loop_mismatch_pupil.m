%%
%matfiles = dir("F:\Alex\matfiles_new\*_MM_*.mat");
matfiles = dir('F:\Alex\matfiles_new\*_MM_*.mat');

%%
MM=[];
RUN=[];
CHANNEL = [];
FR = [];
PUPIL=[];
PUPIL_RUN=[];
for iF=1:numel(matfiles)
data = load(fullfile(matfiles(1).folder,matfiles(iF).name));
if ~ismember('pupil_area',fieldnames(data))
    continue
end
mismatch_trigger = data.vr_data_resampled.MM>0.5;
good_cells = data.sp.cids(data.sp.cgs==2);
%good_cells = data.sp.ks_cluster.cluster_id(startsWith(data.sp.ks_cluster.KSLabel,'good'));
%good_cells = good_cells(ismember(good_cells,data.sp.waveform_metrics.cluster_id));
%%
opt.speed_t=0.05;
all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
true_speed = data.vr_data_resampled.velM;
if iscolumn(true_speed)
    speed=true_speed';
else
    speed=true_speed;
end
filt = gausswin(61); %61 pretty close to what we use in other
filt = filt/sum(filt);
smooth_speed = conv(speed,filt,'same');
%%
run_periods=smooth_speed>opt.speed_t;
run_window=-30:30;
possibles=strfind(run_periods,ones(1,length(run_window)))+floor(.5*length(run_window));

mm_trigs=all_mm_trigs(ismember(all_mm_trigs,possibles));
%%
opt = load_mismatch_opt;
[spikeTimes,~,aux,~,count_vec]=extract_triggered_spikeTimes(data.sp,data.post(mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data.post' ;smooth_speed;data.pupil_area';data.pupil_com'],'aux_win',opt.aux_win);
duration = max(data.post)-min(data.post);
fr = nan(size(good_cells));
for iC=1:numel(good_cells)
    n = nnz(data.sp.clu==good_cells(iC));
    fr(iC)=n/duration;
end


run_ons = strfind(smooth_speed>opt.speed_t,[zeros(1,30),ones(1,30)])+30;
[spikeTimes,~,aux_run,~,count_vec_run]=extract_triggered_spikeTimes(data.sp,data.post(run_ons),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data.post' ;smooth_speed;data.pupil_area';data.pupil_com'],'aux_win',opt.aux_win);
%%

chan_number =data.sp.waveform_metrics.peak_channel(ismember(data.sp.waveform_metrics.cluster_id,good_cells));
[~,sid]=sort(chan_number,'descend');
%%
RUN=cat(1,RUN,count_vec_run);
MM=cat(1,MM,count_vec);
CHANNEL = cat(1,CHANNEL,chan_number);
FR=cat(1,FR,fr);
PUPIL = cat(1,PUPIL,mean(aux(2,:,:),2));
PUPIL_RUN = cat(1,PUPIL_RUN,mean(aux_run(2,:,:)));
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