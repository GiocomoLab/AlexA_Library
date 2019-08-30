
%%
tuning_map_1 = []; 
tuning_map_2 = [];
tuning_map = [];
for iS=1:length(proj_meta)
    dF_valid = proj_meta(iS).dFF(:,proj_meta(iS).vr_data.valid_idx);
    posx_bins = 0:3:300;
    posx_centers = (posx_bins(1:end-1)+posx_bins(2:end))/2;
    posx_bins(1)=-inf;
    posx_bins(end)=inf;
    pos_id = discretize(proj_meta(iS).vr_data.xpos,posx_bins);
    tuning_curve = zeros(size(dF_valid,1),numel(posx_bins)-1);
    t_c_first = tuning_curve;
    t_c_second = tuning_curve;
    first_half = false(size(pos_id));
    first_half(1:round(numel(pos_id)/4))=true;
for iB=1:numel(posx_bins-1)
    tuning_curve(:,iB)=mean(dF_valid(:,pos_id==iB),2);
    t_c_first(:,iB) = mean(dF_valid(:,pos_id == iB & first_half),2);
    t_c_second(:,iB) = mean(dF_valid(:,pos_id == iB & ~first_half),2);
end
    tuning_map = cat(1,tuning_map,tuning_curve);
    tuning_map_1 = cat(1,tuning_map_1,t_c_first);
    tuning_map_2 = cat(1,tuning_map_2,t_c_second);
[~,sid]=max(tuning_curve,[],2);

[~,sid2]=sort(sid);
% figure
% imagesc(tuning_curve(sid2,:))
end
%%

[~,sid]=max(tuning_map_1,[],2);

[~,sid2]=sort(sid);
figure
imagesc(tuning_map_2(sid2,:))
%%
stab = corr(tuning_map_1(:,1:end-1)',tuning_map_2(:,1:end-1)');
stab = diag(stab);
stab_idx = stab>.6;
[~,sid]=max(tuning_map_1(stab_idx,:),[],2);

[~,sid2]=sort(sid);
tmp=tuning_map_2(stab_idx,:);

figure
tmp = tmp(sid2,1:end-1);
win = gausswin(3);
win = win/sum(win);
%tmp = conv2(1,win,tmp);
%tmp = bsxfun(@rdivide,tmp,max(tmp,[],2));
for iR=1:size(tmp,1)
    tmp(iR,:)=smooth(tmp(iR,:),5);
end
ff=max(tmp,[],2);

imagesc(tmp(ff>0,:),[.9 1.5])
ylabel('cell #')
set(gca,'XTick',[1 numel(posx_bins)/2 numel(posx_bins)-1],'XTickLabel',[0 150 300])
xlabel('Distance from start [cm]')
colormap(brewermap([],'Blues'))
colorbar
set(gcf,'Position',[680   368   334   610])
%%
ff=max(tuning_map,[],2);

stab_idx = ff>1.4;
[~,sid]=max(tuning_map_1(stab_idx,:),[],2);

[~,sid2]=sort(sid);
tmp=tuning_map_2(stab_idx,:);
figure
imagesc(tmp(sid2,:),[1 2])