function [spike_times,count_vec,aux_mat]=extractMM(data,good_cells,opt)

if isfield(data,'mismatch_trigger')
    mismatch_trigger = data.mismatch_trigger;
    if size(mismatch_trigger,1) ~=1
        mismatch_trigger=mismatch_trigger';
    end
    if nnz(mismatch_trigger==1)>nnz(mismatch_trigger==0)
        %in some old files mismatch_trigger was actually the move variable,
        %i.e. move ==0 is mismatch
        mismatch_trigger = mismatch_trigger<0.1;
    end
else
    mismatch_trigger = data.vr_data_resampled.MM>0.7;
end
all_mm_trigs=strfind(mismatch_trigger,[0 0 1 1])+2;
if isempty(good_cells)
    good_cells = data.sp.cids(data.sp.cgs==2);
end

if isfield(data,'true_speed')
    true_speed = data.true_speed;
else
    true_speed = data.vr_data_resampled.velM;
end
if iscolumn(true_speed)
    speed=true_speed';
else
    speed=true_speed;
end
filt = gausswin(opt.speed_filt_win); %61 pretty close to what we use in other
filt = filt/sum(filt);
smooth_speed = conv(speed,filt,'same');

run_periods=smooth_speed>opt.speed_t;

possibles=strfind(run_periods,ones(1,length(opt.run_window)))+floor(.5*length(opt.run_window));

mm_trigs=all_mm_trigs(ismember(all_mm_trigs,possibles));
if numel(mm_trigs)>=5
    [spike_times,win,aux_mat,time_idx_used,count_vec]=extract_triggered_spikeTimes(data.sp,data.post(mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data.post' ;smooth_speed],'aux_win',opt.aux_win,'opt',opt);
else
    warning('less than 5 MM events, returning nans')
    spike_times = nan;
    win = nan;
    aux_mat = nan;
    time_idx_used = nan;
    count_vec = nan(numel(good_cells),250);
end