ops = load_default_opt;
ops.trial_range = [-6:9];
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

OAK='/oak/stanford/groups/giocomo/';
OAK = '/Volumes/Samsung_T5';
%%
gain = 0.5;
contrast = 100;
regions = {'VISp','RS','MEC'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA_corrected'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end

savepath = fullfile(OAK,'attialex',['tbtxcorr_decoder_' num2str(gain) '_newNorm_nocutoff']);
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
for iF=45:numel(filenames)
    
    try
        [~,sn]=fileparts(filenames{iF});
        
        for iRep=1:numel(triggers{iF})
            data = load(filenames{iF});
            
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
            %occ = histcounts(data.posx(ismember(data.trial,trials)),ops.xbinedges);
            occ = histcounts(data.posx,ops.xbinedges);
            occ_2=gauss_smoothing_no_taper(occ',23,ops.smoothSigma_dist);
            speed = calcSpeed(data.posx,ops);
            speed = speed./data.trial_gain(data.trial);
            
            trial_this = data.trial(ismember(data.trial,trials));
            posx_this = data.posx(ismember(data.trial,trials));
            speed_this = speed(ismember(data.trial,trials));
            [~,~,posbin] = histcounts(posx_this,ops.xbinedges);
            trial_tmp = trial_this-trial_this(1)+1;
            
            
            
            
            stab = nanmean(nanmean(corrMat(:,1:6,1:6),2),3);
            if nnz(stab>ops.stab_thresh)<5
                continue
            end
            
            
            
            cellID = cellID(stab>ops.stab_thresh);
            region = reg(stab>ops.stab_thresh);
            corrMat = corrMat(stab>ops.stab_thresh,:,:);
            shiftMat = shiftMat(stab>ops.stab_thresh,:,:);
            
            frMat_stable = frMat(stab>ops.stab_thresh,:,:);
            rate_mat = squeeze(nanmean(frMat_stable(:,1:6,:),2))';
            rate_mat = reshape(rate_mat,[1 size(rate_mat)]);
            decode_fr = fr(stab>ops.stab_thresh,ismember(data.trial,trials))*ops.TimeBin;
            postProb = decode_calcBayesPost(decode_fr,rate_mat, occ_2',ops.TimeBin,false);
            postProbUniform = decode_calcBayesPost(decode_fr,rate_mat, ones(1,200),ops.TimeBin,false);
            
            [~,dp] = max(postProb,[],2);
            [~,dpUniform]=max(postProbUniform,[],2);
            
            tmp_e = posx_this' - ops.xbincent(dp);
            correction_idx = abs(tmp_e)>ops.TrackEnd/2;
            tmp_e(correction_idx) = tmp_e(correction_idx)-ops.TrackEnd*sign(tmp_e(correction_idx));
            
            tmp_U = posx_this' - ops.xbincent(dpUniform);
            correction_idx = abs(tmp_U)>ops.TrackEnd/2;
            tmp_U(correction_idx) = tmp_U(correction_idx)-ops.TrackEnd*sign(tmp_U(correction_idx));
            
            
            data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);
            data_out.corrMat = corrMat;
            data_out.shiftMat = shiftMat;
            data_out.time_error = tmp_e;
            data_out.time_errorUniform = tmp_U;
            
            data_out.time_posterior = squeeze(postProb)';
            data_out.posx = posx_this;
            
            data_out.region = region;
            
            
            data_out.trials = trial_this;
            data_out.speed = speed(ismember(data.trial,trials));
        end
    catch ME
        disp(ME.message)
        disp(sprintf('filenr: %d',iF))
        rethrow(ME)
    end
end

