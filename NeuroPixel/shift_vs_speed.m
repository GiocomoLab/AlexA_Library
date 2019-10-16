cellIDX = find(good_cells==cells_to_plot(k));

rr=xcorr(mS(60:100,:));
nR=size(mS,2);
li=zeros(nR*(nR-1)/2,2);
idx = 0;
for iT = 1:212
    for jT=(iT+1):212
        idx = idx+1;
        c_idx = (iT-1)*212+jT;
        [mm,midx]=max(rr(:,c_idx));
        if ismember(midx,[1 81])
           midx = nan;
end
        %delays2(iT,jT)=midx-41;
        %delays2(jT,iT)=-midx+41;
        hold on
        li(idx,:)=[midx-41,trial_speed(iT)-trial_speed(jT)];
    %delays(iT,jT)=finddelay(mS(60:100,iT),mS(60:100,jT));
end
end