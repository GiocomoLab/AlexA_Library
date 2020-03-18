regions = {'VISp','MEC','RS'};
savepath = '/oak/stanford/groups/giocomo/attialex/tbtxcorr_2cmbinShift_speed5';
main_fig = figure;
for iR = 1:numel(regions)
    region = regions{iR};

%shift_dir = sprintf('Z:/giocomo/attialex/images/xcorrv9/%s_0.80_100',region);
matfiles = dir(fullfile(savepath,'*.mat'));
allCorrMat = [];
allM_MEC=[];
allShiftMat = [];
allShiftMatShifted = [];
allCorrMatShifted = [];
baseline_shifts = [];
similarity_score = [];
similarity_score_MEC=[];
cntr = 0;
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    idx_region  = startsWith(data_out.region,region);
    idx = data_out.stability>.4 & startsWith(data_out.region,region)' & data_out.factors'<100;
    if nnz(idx)>2 && nnz(idx)/nnz(idx_region)>.1 
        %baseline_xcorr = load(fullfile(shift_dir,matfiles(iF).name));
        % tmp_peak = baseline_xcorr.peak;
        %         if any(isnan(tmp_peak))
        %             continue
        %         end
        cntr = cntr+1;
        tmp=squeeze(nanmean(data_out.corrMat(idx,:,:)));
        ff=tmp(7:16,1:6);
        similarity_score(end+1) = mean(ff(:));
        
        tmpS = squeeze(nanmean(data_out.shiftMat(idx,:,:)));
        allCorrMat=cat(3,allCorrMat,tmp);
        allCorrMatShifted=cat(3,allCorrMatShifted,squeeze(nanmean(data_out.corrMatShifted(idx,:,:))));
        
        allShiftMat = cat(3,allShiftMat,tmpS);
        allShiftMatShifted = cat(3,allShiftMatShifted,squeeze(nanmean(data_out.shiftMatShifted(idx,:,:))));
        
        %baseline_shifts = cat(1,baseline_shifts,reshape(tmp_peak',1,[]));
        if startsWith(region,'ECT')
            idx_MEC = data_out.stability>.4 & startsWith(data_out.region,'MEC')';
            tmp_MEC = squeeze(nanmean(data_out.corrMat(idx_MEC,:,:)));
            allM_MEC = cat(3,allM_MEC,tmp_MEC);
            ff=tmp_MEC(7:16,1:6);
            similarity_score_MEC(end+1)=mean(ff(:));
        end
    else
        nstab = nnz(idx);
        ntot = nnz(idx_region);
        sprintf('Kicked %s, stab: %d, tot: %d \n',matfiles(iF).name,nstab,ntot)
    end
end
%%
savefig{1} =figure();

subplot(2,2,1)
hold on
imagesc((nanmean(allCorrMat,3)),[0 0.75])
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
imagesc((nanmean(allCorrMatShifted,3)),[0 0.75])
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
imagesc((nanmean(allShiftMat,3)),[-3 3])
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
imagesc((nanmean(allShiftMatShifted,3)),[-3 3])
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
tf = triu(true(num_tr),1);
Y = [];
for ii=1:size(allCorrMat,3)
    tmp = allCorrMat(:,:,ii);
    Y=cat(1,Y,tmp(tf)');
end
[coeff,score,~,~,expl] = pca(Y);
%%
figure
subplot(1,2,1)
corr_per_Site = nanmean(nanmean(allCorrMat(1:6,7:10,similarity_score>.5),1),2);
corr_per_Site_corrected = nanmean(nanmean(allCorrMatShifted(1:6,7:10,similarity_score>.5),1),2);


scatter(corr_per_Site,corr_per_Site_corrected,15,'filled')
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
shifts_per_Site = nanmean(nanmean(allShiftMat(1:6,7:10,similarity_score>.5),1),2);
shifts_per_Site_corrected = nanmean(nanmean(allShiftMatShifted(1:6,7:10,similarity_score>.5),1),2);


scatter(shifts_per_Site,shifts_per_Site_corrected,15,'filled')
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
figure(main_fig)
subplot(1,2,1)
hold on
scatter(corr_per_Site,corr_per_Site_corrected,15,'filled')
subplot(1,2,2)
hold on
scatter(shifts_per_Site,shifts_per_Site_corrected,15,'filled')

end
%%
figure(main_fig)
subplot(1,2,1)
axis image
grid on
hold on
xlim([0.3 0.8])
ylim([0.3 0.8])
xlabel('similarity')
ylabel('speed adjusted similarity')  
subplot(1,2,2)
axis image
grid on
hold on
xlim([-15 5])
ylim([-15 5])
xlabel('shift')
ylabel('speed adjusted shift') 

