idx_dt=1;
idx_pt = 1;
vr_timestamps = vr_time(3:3:end);
vr_timestamps= vr_timestamps-vr_timestamps(1);
tmp = zeros(size(framedata.times));
tmp_idx = 1;
vr_idx = 2;
tmp(1)=vr_timestamps(1);
while tmp_idx<numel(framedata.times)
    if (vr_timestamps(vr_idx)-tmp(tmp_idx))>=0.0318
        tmp(tmp_idx+1)=vr_timestamps(vr_idx);
        tmp_idx=tmp_idx+1;
    end
    vr_idx = vr_idx+1;
end
    