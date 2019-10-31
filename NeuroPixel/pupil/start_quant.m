params = readtable('UniversalParams.xlsx');

speed = calcSpeed(posx,params);
trial_speed = zeros(1,max(trial));
trial_pupil = zeros(1,max(trial));
posWindow = [0 400];
for ii = 1:max(trial)
    idx = posx>posWindow(1) & posx<posWindow(2) & trial == ii;
    tmp = mean(speed(idx));
    tmpP=nanmedian(eye_r_upsampled(idx));
    trial_speed(ii)=tmp;
    trial_pupil(ii)=tmpP;
end