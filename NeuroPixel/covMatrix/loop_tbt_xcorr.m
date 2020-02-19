ops.factors = -.25:0.01:.25;
ops.BinWidth = 2;

ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;

%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops.TimeBin = 0.02;
ops.plotfig = false;

ops.smoothSigma = 4;
ops.maxLag = 20; % in cm

%%
[filenames,triggers] = getFilesCriteria('ECT',100,0.8,'F:/NP_DATA');
savepath = 'F:/temp/tbtxcorr';
if ~isfolder(savepath)
    mkdir(savepath)
end
%%
for iF=1:numel(filenames)
    
    
    
    [~,sn]=fileparts(filenames{iF});
    data = load(filenames{iF});
    
    for iRep=1:numel(triggers{iF})
        data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);
        if isfield(data.anatomy,'parent_shifted')
            reg = data.anatomy.parent_shifted;
        else
            reg = data.anatomy.cluster_parent;
        end
        if iscolumn(reg)
            reg = reg';
        end
        
        
        
        ops.trials=triggers{iF}(iRep)+[-6:9];
        [corrMat,shiftMat,stability]=calculateTrialByTrialXCorr(data,ops);
        data_out.corrMat = corrMat;
        data_out.shiftMat = shiftMat;
        stability2 = nanmean(nanmean(corrMat(:,1:6,1:6),3),2);
        data_out.stability = stability2;

        data_out.region = reg(data.sp.cgs==2);
    end
end

%%
region = 'RSP';

%shift_dir = sprintf('Z:/giocomo/attialex/images/xcorrv9/%s_0.80_100',region);
matfiles = dir('F:\temp\tbtxcorr\*.mat');
allM = [];
allS = [];
baseline_shifts = [];
similarity_score = [];
cntr = 0;
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    idx_region  = startsWith(data_out.region,region);
    idx = data_out.stability>.4 & startsWith(data_out.region,region)';
    if nnz(idx)>2 && nnz(idx)/nnz(idx_region)>.2
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
        allM=cat(3,allM,tmp);
        allS = cat(3,allS,tmpS);
        %baseline_shifts = cat(1,baseline_shifts,reshape(tmp_peak',1,[]));
    else
        nstab = nnz(idx);
            ntot = nnz(idx_region);
            sprintf('Kicked %s, stab: %d, tot: %d \n',matfiles(iF).name,nstab,ntot)
    end
end
%%
savefig{1} =figure();

subplot(1,2,1)
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


subplot(1,2,2)
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
%%
% [coeff,score,~,~,expl] = pca(baseline_shifts);
% [~,sort_idx] = sort(score(:,1),'descend');
% savefig{2}=figure;
% plot(similarity_score(sort_idx),'.')
% ylim([0.2 0.8])
% %[~,sort_idx] = sort(similarity_score,'descend');
% savefig{3}=figure('Position',[680   537   423   441]);
% ha = tight_subplot(8,8);
% for ii=1:numel(sort_idx)
%     axes(ha(ii))
%     tmp = allM(:,:,sort_idx(ii));
%     tmp = flipud(tmp);
%         imagesc(tmp,[0 0.75])
%         set(gca,'XTick',[],'YTick',[])
%         axis image
%         box off
% end

% gain change
tf = triu(true(num_tr),1);
Y = [];
for ii=1:size(allM,3);
    tmp = allM(:,:,ii);
    Y=cat(1,Y,tmp(tf)');
end
[coeff,score,~,~,expl] = pca(Y);
[~,sort_idx_gc] = sort(score(:,1),'descend');
corrmat_sort = allM(:,:,sort_idx_gc);

savefig{2} = figure('Position',[100 100 1000 800],'Renderer','Painters');
ha = tight_subplot(8,10);
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
%%
savefig{3}=figure('Position',[200 200 300 300]);
plot(similarity_score(sort_idx_gc),'k.');
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


%%
for iF=1:numel(savefig)
    set(savefig{iF},'Renderer','Painters');
    saveas(savefig{iF},fullfile('F:/temp/figures',sprintf('tbtx_%s_%d.pdf',region,iF)));
end