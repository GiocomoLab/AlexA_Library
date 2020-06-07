ops.BinWidth = 2;
ops.xbin = ops.BinWidth;
ops.edges = 0:ops.BinWidth:400;
ops.track_length = 400;
ops.nBins = numel(ops.edges)-1;
ops.midpoints = ops.edges(1:end-1)*.5 + ops.edges(2:end)*.5;
ops.search_range=[round(20/ops.BinWidth):round(380/ops.BinWidth)];
%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops.TimeBin = 0.02;
ops.plotfig = false;

ops.smoothSigma = 4;
ops.smoothSigma_time = 0.2;
ops.SpeedCutoff = 2;
ops.maxLag = 20; % in cm
ops.matched_bl_range=[-10:-1];
ops.trial_range = [-6:9];
ops.trials_train = [1:6];
ops.trials_test = [7:16];
ops.xbinedges = 0:ops.xbin:400;
ops.xbincent = ops.xbinedges(1:end-1)+ops.xbin/2;
ops.nBins = numel(ops.xbincent);
ops.edges = ops.xbinedges;

smoothSigma = ops.smoothSigma/ops.BinWidth;
ops.filter = gausswin(floor(smoothSigma*5/2)*2+1);
ops.filter = ops.filter/sum(ops.filter);
%OAK='/oak/stanford/groups/giocomo/';
OAK = '/Volumes/Samsung_T5';
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
savepath = fullfile(OAK,'attialex','tbtxcorr_with_decoder_svm2');
if ~isfolder(savepath)
    mkdir(savepath)
end

%%
% p = gcp('nocreate');
% if isempty(p)
%     p = parpool(12);
% end
%%
for iF=1:numel(filenames)
    
    try
        [~,sn]=fileparts(filenames{iF});
        data = load(filenames{iF});
        
        for iRep=1:numel(triggers{iF})
            if isfield(data.anatomy,'parent_shifted')
                region_this = data.anatomy.parent_shifted;
            else
                region_this = data.anatomy.cluster_parent;
            end
            if iscolumn(region_this)
                region_this = region_this';
            end
            
            
            ops_here = ops;
            
            ops_here.trials=triggers{iF}(iRep)+ops.matched_bl_range;
            
            if ~all(data.trial_contrast(ops_here.trials)==contrast) || ~all(data.trial_gain(ops_here.trials)==1)
                disp('bl trials violating bl condition, skipping this rep')
                continue
            end
            
            
            
            %calculate trial by trial correlation across all bins
            ops_here.trials=triggers{iF}(iRep)+ops.trial_range;
            
            if ~all(data.trial_contrast(ops_here.trials)==contrast)
                error('gain trials violating contrast condition')
                
            end
            if ~all(data.trial_gain((0:3)+triggers{iF}(iRep))==gain)
                error('gain trials wrong gain')
                
            end
            
            
            
            
            
            [corrMat,shiftMat,stability,spatialMap]=calculateTrialByTrialXCorr(data,ops_here);
            
            
            
            idx_m = triu(true(6),1);
            
            baseline_stability = nan(size(corrMat,1),1);
            for iC=1:numel(baseline_stability)
                tmp_m = squeeze(corrMat(iC,1:6,1:6));
                baseline_stability(iC)=nanmean(tmp_m(idx_m));
                
            end
            
            region_good = region_this(data.sp.cgs==2)';
            clus_good = data.sp.cids(data.sp.cgs==2);
            clus_stable = baseline_stability>.5;
            good_cell_idx = clus_stable & startsWith(region_good,regions);
            
            good_cells = clus_good(good_cell_idx);
            
            idxClu = ismember(data.sp.clu,(good_cells));

idxVR=ismember(data.trial,[15:20]) & speed>2;
t_time=data.post(idxVR);
start=min(t_time);
stop=max(t_time);
idxNP=data.sp.st<stop & data.sp.st>=start;
edges = [0:2:400];
[sp,occ,good_cells,track_edges]=getSpikeMatPosition2(data.sp.st(idxClu&idxNP),data.sp.clu(idxClu&idxNP),data.posx(idxVR),data.post(idxVR),'edges',edges);

        sp2=reshape(sp,[size(sp,1) 1 size(sp,2)]);
    
            
            if numel(good_cells)>=5
                corrMat = corrMat(good_cell_idx,:,:);
                shiftMat = shiftMat(good_cell_idx,:,:);
                baseline_stability = baseline_stability(good_cell_idx);
                decode_region = region_good((clus_stable & startsWith(region_good,regions)));
                
                
                %%
                %calculate firing rate
                fr = calcFRVsTime(good_cells,data,ops,ops_here.trials);
                
                % threshold by running speed and trials that we want to look at
                speed = calcSpeed(data.posx,ops);
                %trial_idx = ismember(data.trial,ops_here.trials);
                trial_idx = ismember(data.trial,[21:24]);
                idx = speed>ops.SpeedCutoff & trial_idx;
                X = fr(:,idx);
                trial = data.trial(idx);
                posx = data.posx(idx);
                speed = speed(idx);
                post = decode_calcBayesPost(X, sp2, occ',0.02,false);

                
                % normalize each cell to have FR between 0 and 1
                X = (X-nanmin(X,[],2))./repmat(nanmax(X,[],2)-nanmin(X,[],2),1,size(X,2));
                
                % set nans to 0
                X(isnan(X)) = 0;
                
                % mean subtract
                Xtilde = X - nanmean(X,2);
                
                %prepare y
                num_components = size(Xtilde,1);
                [~,~,posbin] = histcounts(posx,ops.xbinedges);
                posbin(posbin==0) = 1;
                
                %subsample for other classifiers
                sub_sample_idx = false(size(trial));
                sub_sample_idx(1:10:end)=true; %1s p s
                
                
                train_trials = ops_here.trials(ops.trials_train);
                train_trial_idx=ismember(trial,train_trials);
                train_trial_idx = train_trial_idx & sub_sample_idx;
                
                tc = nan(num_components,ops.nBins,1);
                for i = 1:ops.nBins
                    tc(:,i) = mean(Xtilde(:,posbin==i & train_trial_idx),2);
                end
                
                % decode position in gain change trials based on baseline1 trials
                dot_prod = tc' * Xtilde; % predict position
                [~,max_bin] = max(dot_prod);
                pred_pos = ops.xbincent(max_bin);
                
                %Mdl = fitlm(Xtilde(:,train_trial_idx)',posbin(train_trial_idx),'RobustOpts','welsch');
                Mdl = fitcecoc(Xtilde(:,train_trial_idx)',posbin(train_trial_idx));
                %yhat{iFold} = ops.xbincent(predict(Mdl,Xtilde(:,test_trial_idx)'));
                yhat = ops.xbincent(mod(round(predict(Mdl,Xtilde')),ops.track_length/2)+1);
                
                decoder = struct();
                decoder.cluID = good_cells;
                decoder.region = decode_region;
                decoder.pred_pos = pred_pos;
                decoder.true_pos = posx;
                decoder.trial = trial;
                decoder.speed = speed;
                decoder.yhat=yhat;
                
                
                plot_idx = ismember(trial,ops_here.trials(7:10));
                h=figure('Visible','off');
                plot(pred_pos(plot_idx))
                hold on
                plot(yhat(plot_idx));
                saveas(h,fullfile(savepath,sprintf('%s_%d.png',sn,iRep)))
                close(h);
                
                data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);
                data_out.decoder = decoder;
                data_out.corrMat = corrMat;
                data_out.shiftMat = shiftMat;
                data_out.region = decode_region;
                
                data_out.stability_baseline = baseline_stability;
                data_out.trials = ops_here.trials;
            end
        end
    catch ME
        disp(ME.message)
        disp(sprintf('filenr: %d',iF))
    end
end

%%
