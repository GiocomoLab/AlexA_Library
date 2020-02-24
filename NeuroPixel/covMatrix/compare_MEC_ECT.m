savefig{1} =figure();

subplot(2,2,1)
hold on
imagesc((nanmean(allM,3)),[0 0.75])
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
imagesc((nanmean(allS,3)),[-3 3])
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

subplot(2,2,3)
hold on
imagesc((nanmean(allM_MEC,3)),[0 0.75])
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


subplot(2,2,4)
hold on
imagesc((nanmean(allS_MEC,3)),[-3 3])
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
for ii=1:size(allM,3);
    tmp = allM(:,:,ii);
    Y=cat(1,Y,tmp(tf)');
end
[coeff,score,~,~,expl] = pca(Y);
[~,sort_idx_gc] = sort(score(:,1),'descend');
corrmat_sort = allM(:,:,sort_idx_gc);
corrmat_sort_MEC = allM_MEC(:,:,sort_idx_gc);

savefig{2} = figure('Position',[100 100 1000 800],'Renderer','Painters');
ha = tight_subplot(5,6);
for i = 1:size(corrmat_sort,3)
    axes(ha(i)); hold on;
    imagesc(squeeze(corrmat_sort(:,:,i)));
    caxis([0 0.7]);
    %colorbar;
    %axis square;
    %title(sprintf('session %d',i));
    
    % patches indicating gain value
    for tr = 1:num_tr
        patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
            'EdgeColor',get_color(trgain(tr),100));
        patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
            'EdgeColor',get_color(trgain(tr),100));
    end
    axis image

end

savefig{2} = figure('Position',[100 100 1000 800],'Renderer','Painters');
ha = tight_subplot(5,6);
for i = 1:size(corrmat_sort,3)
    axes(ha(i)); hold on;
    imagesc(squeeze(corrmat_sort_MEC(:,:,i)));
    caxis([0 0.7]);
    %colorbar;
    %axis square;
    %title(sprintf('session %d',i));
    
    % patches indicating gain value
    for tr = 1:num_tr
        patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
            'EdgeColor',get_color(trgain(tr),100));
        patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
            'EdgeColor',get_color(trgain(tr),100));
    end
    axis image

end
%%

savefig{3}=figure('Position',[200   200   706   300]);
subplot(1,2,1)
plot(similarity_score(sort_idx_gc),'k.');
hold on
plot(similarity_score_MEC(sort_idx_gc),'r.');
hold on;
yline(0.5)
xlabel('rep (sorted by PC1)');
ylabel('avg corr: (gc + bl_post) to bl_pre','Interpreter','none');
title('Gain change');
num_stab_reps = sum(similarity_score>0.5);
text(max(xlim()),max(ylim()),sprintf('num stable reps:\n%d/%d (%0.1f%%)',...
    num_stab_reps,numel(similarity_score),100*num_stab_reps/numel(similarity_score)),...
    'HorizontalAlignment','right','VerticalAlignment','top');
set(gca,'box','off');
legend({'ECT','MEC'},'Location','SW')
subplot(1,2,2)
scatter(similarity_score,similarity_score_MEC)
axis image
xlim([0.2 .8])
ylim([.2 .8])
xlabel('stability ECT')
ylabel('Stability MEC')
