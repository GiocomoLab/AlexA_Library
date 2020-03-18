region = 'MEC';
savepath = '/oak/stanford/groups/giocomo/attialex/tbtxcorr_2cmbinShift';

%shift_dir = sprintf('Z:/giocomo/attialex/images/xcorrv9/%s_0.80_100',region);
matfiles = dir(fullfile(savepath,'*.mat'));
allCorrMat = [];
allM_MEC=[];
allShiftMat = [];
allShiftMatShifted = [];
allCorrMatShifted = [];
STAB = [];
DEPTH = [];
FACT = [];
SITEID = [];
if ismember(region,{'MEC','ECT'})
    mult = -1;
else
    mult = 1;
end

for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    idx_region  = startsWith(data_out.region,region);
    

        allCorrMat=cat(1,allCorrMat,data_out.corrMat(idx_region,:,:));
        allCorrMatShifted=cat(1,allCorrMatShifted,data_out.corrMatShifted(idx_region,:,:));
        
        allShiftMat = cat(1,allShiftMat,data_out.shiftMat(idx_region,:,:));
        allShiftMatShifted = cat(1,allShiftMatShifted,data_out.shiftMatShifted(idx_region,:,:));
        STAB = cat(1,STAB,data_out.stability(idx_region));
        
        DEPTH = cat(1,DEPTH,data_out.depth(idx_region)'*mult);
        FACT = cat(1,FACT,data_out.factors(idx_region)');
        SITEID = cat(1,SITEID,ones(nnz(idx_region),1)*iF);
       
end
%%
savefig{1} =figure();
IDX = STAB>.6;

subplot(2,2,1)
hold on
imagesc(squeeze(nanmean(allCorrMat(IDX,:,:))),[0 0.75])
axis square
hold on
num_tr=16;
trgain = [ones(1,6), 0.8*ones(1,4), ones(1,6)];
for tr = 1:num_tr
    patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
        'EdgeColor',get_color(trgain(tr),100));
    patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
        'EdgeColor',get_color(trgain(tr),100));
end

for tr = 0.5+[6 10]
    plot([0 num_tr+0.5],[tr tr],'w-');
    plot([tr tr],[0 num_tr+0.5],'w-');
end
%xlim([-num_tr/15 num_tr]); ylim([-num_tr/15 num_tr]);
axis image

colorbar
axis image
box off
set(gca,'XTick',[],'YTick',[]);
subplot(2,2,3)
hold on
imagesc(squeeze(nanmean(allCorrMatShifted(IDX,:,:))),[0 0.75])
axis square
hold on
num_tr=16;
trgain = [ones(1,6), 0.8*ones(1,4), ones(1,6)];
for tr = 1:num_tr
    patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
        'EdgeColor',get_color(trgain(tr),100));
    patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
        'EdgeColor',get_color(trgain(tr),100));
end

for tr = 0.5+[6 10]
    plot([0 num_tr+0.5],[tr tr],'w-');
    plot([tr tr],[0 num_tr+0.5],'w-');
end
%xlim([-num_tr/15 num_tr]); ylim([-num_tr/15 num_tr]);
axis image

colorbar
axis image
box off
set(gca,'XTick',[],'YTick',[]);


subplot(2,2,2)
hold on
imagesc(squeeze(nanmean(allShiftMat(IDX,:,:))),[-3 3])
xline(6.5,'w');
xline(10.5,'w');
yline(6.5,'w');
yline(10.5,'w');
axis image
colorbar
box off
set(gca,'XTick',[],'YTick',[]);
title(region)
for tr = 1:num_tr
    patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
        'EdgeColor',get_color(trgain(tr),100));
    patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
        'EdgeColor',get_color(trgain(tr),100));
end

for tr = 0.5+[6 10]
    plot([0 num_tr+0.5],[tr tr],'w-');
    plot([tr tr],[0 num_tr+0.5],'w-');
end

subplot(2,2,4)
hold on
imagesc(squeeze(nanmean(allShiftMatShifted(IDX,:,:))),[-3 3])
xline(6.5,'w');
xline(10.5,'w');
yline(6.5,'w');
yline(10.5,'w');
axis image
colorbar
box off
set(gca,'XTick',[],'YTick',[]);
title(region)
for tr = 1:num_tr
    patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
        'EdgeColor',get_color(trgain(tr),100));
    patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
        'EdgeColor',get_color(trgain(tr),100));
end

for tr = 0.5+[6 10]
    plot([0 num_tr+0.5],[tr tr],'w-');
    plot([tr tr],[0 num_tr+0.5],'w-');
end

%%

%%
figure
IDX = STAB>.4;
subplot(1,2,1)
corr_per_Cell = nanmean(nanmean(allCorrMat(IDX,1:6,7:10),2),3);
corr_per_Cell_corrected = nanmean(nanmean(allCorrMatShifted(IDX,1:6,7:10),2),3);


scatter(corr_per_Cell,corr_per_Cell_corrected,15,'filled')
axis image
grid on
hold on
plot([-15 15],[-15 15],'k')
xlim([-0 1])
ylim([-0 1])
xlabel('shift')
ylabel('speed adjusted correlation')  
set(gca,'CLim',[-.3 .3])
colorbar

subplot(1,2,2)
shifts_per_cell = nanmean(nanmean(allShiftMat(IDX,1:6,7:10),2),3);
shifts_per_cell_corrected = nanmean(nanmean(allShiftMatShifted(IDX,1:6,7:10),2),3);


scatter(shifts_per_cell,shifts_per_cell_corrected,15,'filled')
axis image
grid on
hold on
plot([-15 15],[-15 15],'k')
xlim([-15 15])
ylim([-15 15])
xlabel('shift')
ylabel('speed adjusted shift')  
set(gca,'CLim',[-.3 .3])
colorbar
%%
figure
subplot(1,2,1)
plot([shifts_per_cell,shifts_per_cell_corrected]',[DEPTH(IDX),DEPTH(IDX)]','k')
hold on
scatter(shifts_per_cell,DEPTH(IDX),15,'filled')
set(gca,'YDir','reverse')
xlim([-20 10])
grid on


subplot(1,2,2)
plot([shifts_per_cell,shifts_per_cell_corrected]',[DEPTH(IDX),DEPTH(IDX)]','k')
hold on
scatter(shifts_per_cell_corrected,DEPTH(IDX),15,'filled')
set(gca,'YDir','reverse')
xlim([-20 10])
grid on
%%
figure('Position',[941   206   604   641])
subplot(1,2,1)
scatter(shifts_per_cell,DEPTH(IDX),15,FACT(IDX),'filled')
set(gca,'YDir','reverse')
xlim([-20 10])
grid on
title(region);
ylabel('Depth')
xlabel('shift')
colorbar
aid = unique(SITEID);
midpoints= -2200:200:1500;
half_width = 200;

idx_value = STAB>.4;

SHIFTS = nanmean(nanmean(allShiftMat(:,1:6,7:10),2),3);
SHIFTS_CORRECTED = nanmean(nanmean(allShiftMatShifted(:,1:6,7:10),2),3);

allSites = nan(numel(aid),numel(midpoints));
allSites_corrected = allSites;
for ii=1:numel(aid)
    idx = SITEID==ii;
    valid_idx = idx_value & idx;
    vals = nan(1,numel(midpoints));
    vals_c = vals;
    for iv=1:numel(vals)
    idx = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width & valid_idx;
    vals(iv)=nanmean(SHIFTS(idx));
    vals_c(iv)=nanmean(SHIFTS_CORRECTED(idx));
    end
    %plot(vals,edges(1:end-1))
    allSites(ii,:)=vals;
    allSites_corrected(ii,:)=vals_c;
end




subplot(1,2,2)
%plot(nanmean(allSites),edges(1:end-1))
err = nanstd(allSites)/sqrt(numel(aid));
errorbar(nanmean(allSites),midpoints,[],[],err,err)
hold on
err_c = nanstd(allSites_corrected)/sqrt(numel(aid));
errorbar(nanmean(allSites_corrected),midpoints,[],[],err,err)
ylabel('Depth')
set(gca,'YDir','reverse')
xlabel('Shift')
legend({'Shifts','Corrected'})
