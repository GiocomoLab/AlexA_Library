ops = load_default_opt;
ops.trial_range = [0:15];
ops.BinWidth = 2;
ops.edges = 0:ops.BinWidth:400;
ops.xbinedges = ops.edges;
ops.xbincent = .5*ops.edges(1:end-1)+.5*ops.edges(2:end);
ops.SpatialBin = ops.BinWidth;
ops.nBins = numel(ops.edges)-1;
ops.smoothSigma=ops.smoothSigma_dist;
smoothSigma = ops.smoothSigma/ops.BinWidth;
ops.filter = gausswin(floor(smoothSigma*5/2)*2+1);
ops.filter = ops.filter/sum(ops.filter);
ops.max_lag = 30;
ops.maxLag = ops.max_lag;
ops.stab_thresh = 0.5;
ops.trials_train = 1:6;
ops.SpeedCutoff = -1;

%OAK='/oak/stanford/groups/giocomo/';
OAK = '/Volumes/T7';
%%
% gain = 0.5;
% contrast = 100;
regions = {'MEC'};
filenames = {};
triggers = {};
% for iR = 1:numel(regions)
%     
%     [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA_corrected'));
%     filenames=cat(2,filenames,tmp1);
%     triggers = cat(2,triggers,tmp2);
% end
fn = dir(fullfile(OAK,'attialex','NP_DATA_corrected','np*'))
for iF=1:numel(fn)
    if ~contains(fn(iF).name,{'dark','mismatch','playback'})
        filenames=cat(2,filenames,fullfile(fn(iF).folder,fn(iF).name));
        triggers = cat(2,triggers,{5});
    end
end
savepath = fullfile(OAK,'attialex',['tbtxcorr_decoder_baseline']);
if ~isfolder(savepath)
    mkdir(savepath)
end
save(fullfile(OAK,'attialex','parameters.mat'),'ops');
%
% p = gcp('nocreate');
% if isempty(p)
%     p = parpool(12);
% end
%%
for iF=1:numel(filenames)
    data = load(filenames{iF});

    if ~isfield(data,'anatomy')
        continue
    end
    try
        [~,sn]=fileparts(filenames{iF});
        
        for iRep=1:numel(triggers{iF})
            
            if isfield(data.anatomy,'parent_shifted')
                reg = data.anatomy.parent_shifted;
            else
                reg = data.anatomy.cluster_parent;
            end
            if iscolumn(reg)
                reg = reg';
            end
            
            reg=reg(data.sp.cgs==2);
            reg_orig = data.anatomy.cluster_parent((data.sp.cgs==2));
            if iscolumn(reg_orig)
                reg_orig = reg_orig';
            end
            
            
            
            %calculate trial by trial correlation across all bins
            trials=triggers{iF}(iRep)+ops.trial_range;
            ops_here = ops;
            ops_here.trials = trials;
            cellID = data.sp.cids(data.sp.cgs==2);
            reg(startsWith(reg,'VISpm'))={'VISm'};
            reg(startsWith(reg,'RSPa'))={'RA'};
            reg = reg(startsWith(reg,regions));
            cellID = cellID(startsWith(reg,regions));
            
            if numel(cellID)<5
                continue;
            end
            [corrMat,frMat,shiftMat]=trialCorrMat(cellID,trials,data,ops);
            fr = calcFRVsTime(cellID,data,ops);
            speed = calcSpeed(data.posx,ops);
            
            lick_per_bin = histcounts(data.lickt,data.post)/ops.TimeBin;
            lick_per_bin = [0;lick_per_bin'];
            lick_rate = gauss_smoothing(lick_per_bin,10);
            
            [~,~,posx_bin] = histcounts(data.posx,ops.xbincent);
            trial_speed = zeros(numel(trials),numel(ops.xbincent));
            trial_lick = trial_speed;
            for iT=1:numel(trials)
                for iPos = 1:numel(ops.xbincent)
                idx = data.trial ==trials(iT) & posx_bin==iPos;
                trial_speed(iT,iPos)=mean(speed(idx));
                trial_lick(iT,iPos) = mean(lick_rate(idx));
                end
            end
            speed_orig = speed;
            
            fr = fr(:,speed>ops.SpeedCutoff);
            trial = data.trial(speed>ops.SpeedCutoff);
            posx = data.posx(speed>ops.SpeedCutoff);
            lick_rate = lick_rate(speed>ops.SpeedCutoff);

            speed = speed(speed>ops.SpeedCutoff);
            
            % extract firing rate map and position data for trials from this rep
            X = fr(:,ismember(trial,trials));
            mm=max(X,[],2);
            frN=frMat;
            for iT=1:16
                for iP=1:200
                    frN(:,iT,iP)=frMat(:,iT,iP)./mm;
                end
            end
            frMat = frN;
            X=X./mm;
            stab = nanmean(nanmean(corrMat(:,1:6,1:6),2),3);
            if nnz(stab>ops.stab_thresh)<5
                continue
            end
            
            cellID = cellID(stab>ops.stab_thresh);
            region = reg(stab>ops.stab_thresh);
            frMat = frMat(stab>ops.stab_thresh,:,:);
            
            X=X(stab>ops.stab_thresh,:);
            
            Xtilde = X;
            
            
            
            
            corrMat = corrMat(stab>ops.stab_thresh,:,:);
            shiftMat = shiftMat(stab>ops.stab_thresh,:,:);
            score_mat = zeros(3,size(frMat,2),size(frMat,3));
            for iFold = 1:numel(trials)
                take_idx = true(1,numel(trials));
                
            
                    take_idx(iFold)=false;
                    %calculate tuning curve based on 'training trials'
             
                
                tuning_curve = squeeze(mean(frMat(:,take_idx,:),2));
                vn=vecnorm(tuning_curve,2,1);
                tc_n = tuning_curve./vn;
                Xt=squeeze(frMat(:,iFold,:));
                Xt=Xt./vecnorm(Xt,2,1);
                Xt(isnan(Xt))=0;
                %find closest poin on trajectory
                dot_prod = tc_n' * Xt;
                
                [dist,max_bin] = max(dot_prod);
                tmp_e = ops.xbincent - ops.xbincent(max_bin);
                correction_idx = abs(tmp_e)>ops.TrackEnd/2;
                tmp_e(correction_idx) = tmp_e(correction_idx)-ops.TrackEnd*sign(tmp_e(correction_idx));
                dec_position = ops.xbincent(max_bin);
                score_mat(:,iFold,:)=[tmp_e;dist;dec_position];
            end
            
            
            trial_this = trial(ismember(trial,trials));
            posx_this = posx(ismember(trial,trials));
            speed_this = speed(ismember(trial,trials));
            lick_this = lick_rate(ismember(trial,trials));
            posx_shifted = mod(posx_this+200,400);
            [~,~,posbin_orig] = histcounts(posx_this,ops.xbinedges);
            
            X = fr(:,ismember(trial,trials));
            X = zscore(X,[],2);
            X=X(stab>ops.stab_thresh,:);
            frMat = zeros(size(X,1),16,200);
            for iT=1:numel(trials)
                for iPos = 1:numel(ops.xbincent);
                    idx = posbin_orig==iPos & trial_this==trials(iT);
                    frMat(:,iT,iPos)=nanmean(X(:,idx),2);
                end
            end
            
            score_mat2 = zeros(2,size(X,2));
            for iFold = 1:numel(trials)
                take_idx = true(1,numel(trials));
                
            
                take_idx(iFold)=false;
                trial_idx = trial_this ==trials(iFold);
                Xt = X(:,trial_idx)';
                tuning_curve = squeeze(mean(frMat(:,take_idx,:),2))';

                [~,iBin]=pdist2(tuning_curve,Xt,'euclidean','Smallest',1);
                tmp_e2 = posx_this(trial_idx) - ops.xbincent(iBin)';
                correction_idx = abs(tmp_e2)>ops.TrackEnd/2;
                tmp_e2(correction_idx) = tmp_e2(correction_idx)-ops.TrackEnd*sign(tmp_e2(correction_idx));
                dec_position = ops.xbincent(iBin);
                score_mat2(1,trial_idx)=tmp_e2;
                score_mat2(2,trial_idx)=dec_position;
            end
            
          
            data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);
            data_out.corrMat = corrMat;
            
            data_out.scoreMat = score_mat;
            data_out.scoreMat2 = score_mat2;
            data_out.trial_lick = trial_lick;
            data_out.trial_speed = trial_speed;
            
            
            data_out.posx = posx_this';
            data_out.speed = speed_this';
            data_out.lick = lick_this;
            data_out.region = region;
            
            
            data_out.trials = trial_this';
        end
    catch ME
        disp(ME.message)
        disp(sprintf('filenr: %d',iF))
        rethrow(ME)
    end
end

%%
fi = dir(fullfile(savepath,'*.mat'));
idx_win = 110:130;
fastest = [];
slowest = [];
behavior_fast = [];
behavior_slow = [];
for iF=1:numel(fi)
    data = load(fullfile(fi(iF).folder,fi(iF).name));
    %trial_speed = zscore(data.trial_speed,[],'all');
    %trial_lick = zscore(data.trial_lick,[],'all');
av_speed = nanmean(data.trial_speed(:,idx_win),2);%-nanmean(data.trial_speed(:,idx_win-20),2);
av_lick = nanmean(data.trial_lick(:,idx_win),2);
av_error = squeeze(nanmean((data.scoreMat(1,:,idx_win)),3))';
if sum(av_lick)==0
    continue
end
%[~,idx]=min(av_lick);
[~,idx]= min(av_error);
slowest = cat(2,slowest,squeeze(data.scoreMat(3,idx,:)));
tmp = [data.trial_speed(idx,:);data.trial_lick(idx,:)];
behavior_slow = cat(3,behavior_slow,tmp);

%[~,idx]=max(av_lick);
[~,idx]= max(av_error);
fastest = cat(2,fastest,squeeze(data.scoreMat(3,idx,:)));
tmp = [data.trial_speed(idx,:);data.trial_lick(idx,:)];
behavior_fast = cat(3,behavior_fast,tmp);
end
%%
figure
subplot(1,3,1)
plot(ops.xbincent,mean(slowest,2))
hold on
plot(ops.xbincent,mean(fastest,2))
legend({'low lick rate','high lick rate'},'Location','NW')
title('decoded position')
plot([0,400],[0,400],'k--')
subplot(1,3,2)

plot(squeeze(nanmean(behavior_slow(1,:,:),3)))
hold on
plot(squeeze(nanmean(behavior_fast(1,:,:),3)))
title('running speed')
legend({'low lick rate','high lick rate'},'Location','NW')

subplot(1,3,3)

plot(squeeze(nanmean(behavior_slow(2,:,:),3)))
hold on
plot(squeeze(nanmean(behavior_fast(2,:,:),3)))
title('lick rate')
legend({'low lick rate','high lick rate'},'Location','NW')
