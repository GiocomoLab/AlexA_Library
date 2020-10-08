function preview_MM(data)
mismatch_trigger = data.vr_data_resampled.MM>0.5;
good_cells = data.sp.cids(data.sp.cgs==2);
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
count_vecN=count_vec./max(count_vec,[],2);
count_vec_runN = count_vec_run./max(count_vec_run,[],2);
figure
subplot(3,2,[1 3])
imagesc(opt.time_vecs,1:size(count_vec,1),count_vecN(sid,:));
ylabel('ventral -> dorsal')

subplot(3,2,[2 4])
imagesc(opt.time_vecs,1:size(count_vec,1),count_vec_runN(sid,:));

subplot(3,2,5)
plot(opt.time_vecs,nanmean(count_vecN));
title('MM')
subplot(3,2,6)
plot(opt.time_vecs,nanmean(count_vec_runN));
title('Run')
end
