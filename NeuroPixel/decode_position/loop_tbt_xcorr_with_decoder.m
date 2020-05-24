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
OAK='/oak/stanford/groups/giocomo/';
%OAK = '/Volumes/Samsung_T5';
%%
gain = 1.0;
contrast = 100;
regions = {'MEC'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA_corrected'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end
triggers_new = triggers;
for iT=1:numel(triggers)
    nT=numel(triggers{iT});
    triggers_new{iT}=[triggers{iT}(1) triggers{iT}(round(nT/2))];
end
triggers=triggers_new;
savepath = fullfile(OAK,'attialex','tbtxcorr_decoder');
if ~isfolder(savepath)
    mkdir(savepath)
end
save(fullfile(OAK,'attialex','parameters.mat'),'ops');
%
p = gcp('nocreate');
if isempty(p)
    p = parpool(12);
end
%%
parfor iF=1:numel(filenames)
    
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
            
            [corrMat,frMat,shiftMat]=trialCorrMat(cellID,trials,data,ops);
            

            stab = nanmean(nanmean(corrMat(:,1:6,1:6),2),3);
            if nnz(stab>ops.stab_thresh)<5
                continue
            end
            
            cellID = cellID(stab>ops.stab_thresh);
            region = reg(stab>ops.stab_thresh);
            frMat = frMat(stab>ops.stab_thresh,:,:);
            
            score_mat = zeros(2,size(frMat,2),size(frMat,3));
            for iFold = 1:numel(trials)
                take_idx = true(1,numel(trials));
                
                if iFold<6
                    take_idx(iFold)=false;
                    %calculate tuning curve based on 'training trials'
                end
                take_idx(7:end)=false;
                
                tuning_curve = squeeze(mean(frMat(:,take_idx,:),2));
                vn=vecnorm(tuning_curve,2,2);
                tc_n = tuning_curve./vn;
                Xt=squeeze(frMat(:,iFold,:));
                Xt=Xt./vecnorm(Xt,2,2);
                %find closest poin on trajectory
                dot_prod = tuning_curve' * Xt;
                
                [dist,max_bin] = max(dot_prod);
                tmp_e = ops.xbincent - ops.xbincent(max_bin);
                correction_idx = abs(tmp_e)>ops.TrackEnd/2;
                tmp_e(correction_idx) = tmp_e(correction_idx)-ops.TrackEnd*sign(tmp_e(correction_idx));
                score_mat(:,iFold,:)=[tmp_e;dist];
            end
            
 
           
            
            
            
            
            data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);
            data_out.corrMat = corrMat;
            data_out.shiftMat = shiftMat;
            
            
            
            data_out.region = reg;
            data_out.region_orig = reg_orig;
            
           
            data_out.trials = trials;
        end
    catch ME
        disp(ME.message)
        disp(sprintf('filenr: %d',iF))
        rethrow(ME)
    end
end

