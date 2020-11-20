%%
matfiles = dir("F:\Alex\matfiles_new\*_MM_*.mat");
%%
MM=[];
RUN=[];
CHANNEL = [];
for iF=1:13
data = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
mismatch_trigger = data.vr_data_resampled.MM>0.5;
good_cells = data.sp.cids(data.sp.cgs==2);
good_cells = data.sp.ks_cluster.cluster_id(startsWith(data.sp.ks_cluster.KSLabel,'good'));
good_cells = good_cells(ismember(good_cells,data.sp.waveform_metrics.cluster_id));
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
[spikeTimes,~,aux,~,count_vec]=extract_triggered_spikeTimes(data.sp,data.post(mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data.post' ;smooth_speed],'aux_win',opt.aux_win);



run_ons = strfind(smooth_speed>opt.speed_t,[zeros(1,30),ones(1,30)])+30;
[spikeTimes,~,aux,~,count_vec_run]=extract_triggered_spikeTimes(data.sp,data.post(run_ons),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data.post' ;smooth_speed],'aux_win',opt.aux_win);
%%

chan_number =data.sp.waveform_metrics.peak_channel(ismember(data.sp.waveform_metrics.cluster_id,good_cells));
[~,sid]=sort(chan_number,'descend');
%%
RUN=cat(1,RUN,count_vec_run);
MM=cat(1,MM,count_vec);
CHANNEL = cat(1,CHANNEL,chan_number);
end
%%
valid_idx = mean(RUN,2)>1 & mean(MM,2)>1;
[~,sid]=sort(CHANNEL(valid_idx),'descend');
MM_smooth = smoothdata(MM(valid_idx,:),2,'gaussian',10);
RUN_smooth = smoothdata(RUN(valid_idx,:),2,'gaussian',10);
count_vecN=MM_smooth./max(MM_smooth,[],2);
count_vec_runN = RUN_smooth./max(RUN_smooth,[],2);
figure
subplot(3,2,[1 3])
imagesc(opt.time_vecs,1:nnz(valid_idx),count_vecN(sid,:));
ylabel('ventral -> dorsal')

subplot(3,2,[2 4])
imagesc(opt.time_vecs,1:nnz(valid_idx),count_vec_runN(sid,:));

subplot(3,2,5)
plot(opt.time_vecs,nanmean(count_vecN));
title('MM')
subplot(3,2,6)
plot(opt.time_vecs,nanmean(count_vec_runN));
title('Run')