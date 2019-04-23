function [act_smoothed]=get_smoothed_ROIs(act,win,fps)
% highly useless little function thingy... "it saves space - that's what makes it useful" quote Mr. P.


act_smoothed=zeros(size(act));
nRois=size(act,1);
    for gnd=1:nRois
        cur_act=psmooth(act(gnd,:),win,fps);
        act_smoothed(gnd,:)=cur_act/median(cur_act);
    end

end