%load('F:\G5\1207_mismatch_1\1207_mismatch_1.mat')
speed_t=0.05;
% figure('Name',filenames{iF});; plot(speed)
%
if size(mismatch_trigger,1) ~=1
    mismatch_trigger=mismatch_trigger';
end

if nnz(mismatch_trigger>.5)>nnz(mismatch_trigger<.5)
    warning('flipping mismatch trigger');
    all_mm_trigs=strfind(mismatch_trigger<0.1,[0 0 1 1])+2;

else
    
    all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
end

speed=true_speed';

run_periods=smooth(speed,25)>speed_t;
run_window=-30:30;
possibles=strfind(run_periods',ones(1,length(run_window)))+floor(.5*length(run_window));


mm_trigs=all_mm_trigs(ismember(all_mm_trigs,possibles));
possibles=randsample(possibles,100);
if all(size(post)==size(speed))
    speed = speed';
    %posx = posx;
end
%% MM
[spike_mat,win,adata]=extract_triggered_spikeTimes(sp,post(mm_trigs),'win',[-4 4],'aux',[post'; speed; posx'],'aux_win',[-200 200]);
[spike_mat_random,~,adata_random]=extract_triggered_spikeTimes(sp,post(possibles),'win',[-4 4],'aux',[post'; [speed]],'aux_win',[-200 200]);
[spike_mat_all,~,adata_all]=extract_triggered_spikeTimes(sp,post(all_mm_trigs),'win',[-4 4],'aux',[post'; [speed]],'aux_win',[-200 200]);