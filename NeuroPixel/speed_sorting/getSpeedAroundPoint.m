function trial_speed = getSpeedAroundPoint(speed_raw,posx,trial,ops,point,window)

        trial_speed = zeros(1,numel(ops.trials));
        window = window+point;
        for ii = 1:numel(ops.trials)
            idx = posx>window(1) & posx<window(2) & trial == ops.trials(ii);
            tmp = mean(speed_raw(idx));
    
            trial_speed(ii)=tmp;
        end
end