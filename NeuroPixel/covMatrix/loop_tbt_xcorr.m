ops.BinWidth = 2;

ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.midpoints = ops.edges(1:end-1)*.5 + ops.edges(2:end)*.5;
ops.search_range=[round(20/ops.BinWidth):round(380/ops.BinWidth)];
%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops.TimeBin = 0.02;
ops.plotfig = false;

ops.smoothSigma = 4;
ops.maxLag = 20; % in cm
ops.matched_bl_range=[-10:-1];
ops.trial_range = [-6:9];


smoothSigma = ops.smoothSigma/ops.BinWidth;
ops.filter = gausswin(floor(smoothSigma*5/2)*2+1);
ops.filter = ops.filter/sum(ops.filter);
OAK='/oak/stanford/groups/giocomo/';

%%
gain = 0.8;
contrast = 100;
regions = {'VISp','RS','MEC'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
[tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA'));
filenames=cat(2,filenames,tmp1);
triggers = cat(2,triggers,tmp2);
end
savepath = '/oak/stanford/groups/giocomo/attialex/tbtxcorr_with_shortbaseline2';
shiftDir = fullfile(OAK,'attialex','speed_filtered_new_22binspace_5binspeed2');
if ~isfolder(savepath)
    mkdir(savepath)
end
shift_ops = load(fullfile(shiftDir,'parameters.mat'));
shift_ops = shift_ops.ops;
%%
p = gcp('nocreate');
if isempty(p)
    p = parpool(12);
end
%%
parfor iF=1:numel(filenames)
    
    try
    [~,sn]=fileparts(filenames{iF});
    data = load(filenames{iF});
    
    for iRep=1:numel(triggers{iF})
        if isfield(data.anatomy,'parent_shifted')
            reg = data.anatomy.parent_shifted;
        else
            reg = data.anatomy.cluster_parent;
        end
        if iscolumn(reg)
            reg = reg';
        end
        
        
        ops_here = ops;
        
        ops_here.trials=triggers{iF}(iRep)+ops.matched_bl_range;
        
        if ~all(data.trial_contrast(ops_here.trials)==contrast) || ~all(data.trial_gain(ops_here.trials)==1)
            disp('bl trials violating bl condition, skipping this rep')
            continue
        end
        
        [corrMatBL,shiftMatBL,~,~]=calculateTrialByTrialXCorr(data,ops_here);
        
        
        %calculate trial by trial correlation across all bins
        ops_here.trials=triggers{iF}(iRep)+ops.trial_range;
        
        if ~all(data.trial_contrast(ops_here.trials)==contrast)
            error('gain trials violating contrast condition')
            
        end
        if ~all(data.trial_gain((0:3)+triggers{iF}(iRep))==gain)
            error('gain trials wrong gain')
            
        end
        
        %calculate spatial firing rate
        pp=findPeakFiringRate(data,ops_here);
        
        %find max location for each cell
        nC=numel(pp.region);
        window=[-9:10];
        window_length = numel(window);
        bins2correlate = zeros(nC,window_length);
        for ii=1:nC
            bins2correlate(ii,:)=pp.maxBin_index(ii)+window;
        end
        
        [corrMat,shiftMat,stability,spatialMap]=calculateTrialByTrialXCorr(data,ops_here);
        
        %calculate trial by trial correlation around max locations
        [corrMatPartial,shiftMatPartial,stabilityPartial]=calculateTrialByTrialXCorr(data,ops_here,bins2correlate);
        
        %load shift factors
        shift_data = load(fullfile(shiftDir,[sn '.mat']));
        factors = nanmean(shift_data.all_factors);
        
       
        
        % prepare to shift spatial maps according to factors
        good_idx = ismember(data.sp.clu,data.sp.cids(data.sp.cgs==2));
        clu_tmp = data.sp.clu(good_idx);
        st_tmp = data.sp.st(good_idx);
        [uClu,~,clus]=unique(clu_tmp);
        nClu = numel(uClu);
        [~,sr] = calcSpeed(data.posx,ops_here);
        if ~isfield(shift_ops,'speed_filter')
            speed_filter = shift_ops.filter;
        else
            speed_filter = shift_ops.speed_filter;
        end
        speed = conv(sr,speed_filter,'same');
        trial_sorted = nan(size(data.trial));
        trialMap = nan(1,numel(data.trial_gain));
        
        cntr = 1;
        for iT =1:numel(data.trial_gain)
            if ismember(iT,ops_here.trials)
                trialMap(iT)=cntr;
                cntr=cntr+1;
            end
        end
        for iT=1:numel(trial_sorted)
            trial_sorted(iT)=trialMap(data.trial(iT));
        end
        if nClu ~= size(corrMat,1)
            
            good_cells = data.sp.cids(data.sp.cgs==2);
            
            all_good = ismember(good_cells,uClu);
            factors = factors(all_good);
            corrMat = corrMat(all_good,:,:);
            shiftMat = shiftMat(all_good,:,:);
            corrMatPartial = corrMatPartial(all_good,:,:);
            shiftMatPartial = shiftMatPartial(all_good,:,:);
            bins2correlate = bins2correlate(all_good,:);
            pp.maxBin_index = pp.maxBin_index(all_good);
            pp.maxLoc = pp.maxLoc(all_good);
            pp.depth = pp.depth(all_good);
            pp.region = pp.region(all_good);
            if ~isempty(pp.subregion)
                 pp.subregion = pp.subregion(all_good);
            end
            
            %disp(sprintf('%s has not same number of clusters for correlation and shift, skipping',sn))
            %continue
        end
        
        %
        spMapShifted=shiftAllMapsByFactor(ops_here,clus,st_tmp,nClu,data.posx,data.post,trial_sorted,speed,factors);
        % calculate trial by trial correlation
        [corrMatShifted,shiftMatShifted]=spMapXcorr(spMapShifted,ops_here.maxLag,ops_here.BinWidth);
        data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);

        data_out.corrMat = corrMat;
        data_out.corrMatBL = corrMatBL;
        data_out.shiftMat = shiftMat;
        data_out.shiftMatBL = shiftMatBL;
        data_out.corrMatPartial = corrMatPartial;
        data_out.shiftMatPartial = shiftMatPartial;
        data_out.corrMatShifted = corrMatShifted;
        data_out.shiftMatShifted = shiftMatShifted;
        data_out.bins2correlate = bins2correlate;
        data_out.speedMat = pp.speedMat;
        data_out.maxInd = pp.maxBin_index;
        data_out.maxLoc = pp.maxLoc;
        data_out.depth = pp.depth;
        data_out.region = pp.region;
        data_out.subregion = pp.subregion;
        stability2 = nanmean(nanmean(corrMat(:,1:6,1:6),3),2);
        data_out.stability = stability2;
        data_out.factors = factors;
        data_out.trials = ops_here.trials;
    end
    catch ME
        disp(ME.message)
        disp(sprintf('filenr: %d',iF))
    end
end

%%
region = 'MEC';

%shift_dir = sprintf('Z:/giocomo/attialex/images/xcorrv9/%s_0.80_100',region);
%matfiles = dir(fullfile(savepath,'*.mat'));
matfiles = dir('/Volumes/Samsung_T5/tbtxcorr_2cmbinShift_speed5/*.mat')
allCorrMat = [];
allM_MEC=[];
allShiftMat = [];

baseline_shifts = [];
similarity_score = [];
similarity_score_MEC=[];
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
        allCorrMat=cat(3,allCorrMat,tmp);
        
        allShiftMat = cat(3,allShiftMat,tmpS);
        
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

subplot(1,2,1)
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


subplot(1,2,2)
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


%%


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
for ii=1:size(allCorrMat,3);
    tmp = allCorrMat(:,:,ii);
    Y=cat(1,Y,tmp(tf)');
end
[coeff,score,~,~,expl] = pca(Y);
[~,sort_idx_gc] = sort(score(:,1),'descend');
corrmat_sort = allCorrMat(:,:,sort_idx_gc);

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

figure
imagesc(squeeze(mean(allShiftMat(:,:,similarity_score>.5),3)),[-5 5])

%%
% for iF=1:numel(savefig)
%     set(savefig{iF},'Renderer','Painters');
%     saveas(savefig{iF},fullfile('F:/temp/figures',sprintf('tbtx_%s_%d.pdf',region,iF)));
% end
