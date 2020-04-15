ops = load_default_opt;
ops.trial_range = [-6:9];
ops.BinWidth = 2;
ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.smoothSigma=ops.smoothSigma_dist;
smoothSigma = ops.smoothSigma/ops.BinWidth;
ops.filter = gausswin(floor(smoothSigma*5/2)*2+1);
ops.filter = ops.filter/sum(ops.filter);
ops.maxLag = ops.max_lag;
OAK='/oak/stanford/groups/giocomo/';
%OAK = '/Volumes/Samsung_T5';
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
savepath = fullfile(OAK,'attialex','tbtxcorr_v3');
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
            
            reg=reg(data.sp.cgs==2);
            
            
            
            
            %calculate trial by trial correlation across all bins
            trials=triggers{iF}(iRep)+ops.trial_range;
            ops_here = ops;
            ops_here.trials = trials;
            cellID = data.sp.cids(data.sp.cgs==2);
            
            [corrMat,b,shiftMat]=trialCorrMat(cellID,trials,data,ops);
            if ~all(data.trial_contrast(trials)==contrast)
                error('gain trials violating contrast condition')
                
            end
            if ~all(data.trial_gain((0:3)+triggers{iF}(iRep))==gain)
                error('gain trials wrong gain')
                
            end
            
            
            
            
            %load shift factors
            if ~isfile(fullfile(shiftDir,[sn '.mat']))
                sprintf('no speed file for %s',sn)
            end
            shift_data = load(fullfile(shiftDir,[sn '.mat']));
            %tmp = shift_data.all_stability;
            %factors = nanmean(shift_data.all_factors);
            tmp = shift_data.all_stability;
            tmp = tmp>.5;
            tmp_f = shift_data.all_factors;
            tmp_f(tmp)=nan;
            factors = nanmean(tmp_f);
            
            
            % prepare to shift spatial maps according to factors
            good_idx = ismember(data.sp.clu,data.sp.cids(data.sp.cgs==2));
            clu_tmp = data.sp.clu(good_idx);
            st_tmp = data.sp.st(good_idx);
            [uClu,~,clus]=unique(clu_tmp);
            nClu = numel(uClu);
            [~,sr] = calcSpeed(data.posx,ops);
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
                if ismember(iT,trials)
                    trialMap(iT)=cntr;
                    cntr=cntr+1;
                end
            end
            for iT=1:numel(trial_sorted)
                trial_sorted(iT)=trialMap(data.trial(iT));
            end
            
            
            good_cells = data.sp.cids(data.sp.cgs==2);
            
            all_good = ismember(good_cells,uClu);
            factors = factors(all_good);
            corrMat = corrMat(all_good,:,:);
            shiftMat = shiftMat(all_good,:,:);
            reg = reg(all_good);
            
            
            spMapShifted=shiftAllMapsByFactor(ops_here,clus,st_tmp,nClu,data.posx,data.post,trial_sorted,speed,factors);
            % calculate trial by trial correlation
            [corrMatShifted,shiftMatShifted]=spMapXcorr(spMapShifted,ops_here.maxLag,ops_here.BinWidth);
            %v2 add -delay factor to each spike time, then run same code as
            %above
            %         st_adjusted = data.sp.st;
            %         st_old = data.sp.st;
            %         for iFact=1:numel(factors)
            %             if isnan(factors(iFact))
            %                 continue
            %             end
            %             IDX = data.sp.clu == uClu(iFact);
            %             this_fact = factors(iFact);
            %             st_adjusted(IDX)=st_adjusted(IDX)+this_fact;
            %         end
            %         data.sp.st = st_adjusted;
            %         [corrMatS2,b,shiftMatS2]=trialCorrMat(cellID,trials,data,ops);
            %         data.sp.st = st_old;
            corrMatS2 = nan(numel(factors),numel(trials),numel(trials));
            shiftMatS2 = corrMatS2;
            posx_orig = data.posx;
            %%
            for iC = 1:numel(cellID)
                if ~isnan(factors(iC))
                    data.posx=posx_orig;
                    tmp = posx_orig+factors(iC)*speed;
                    tmp = mod(tmp,400);
                    data.posx = tmp;
                    %         for iT = 1:numel(trials)
                    %             t_idx = data.trial == trials(iT);
                    %             tmp = posx_orig(t_idx) + factors(iC)*speed(t_idx);
                    %             tmp(tmp<0)=0;
                    %             tmp(tmp>400)=400;
                    %             data.posx(t_idx)=tmp;
                    %         end
                    
                    [tmpC,b,tmpS]=trialCorrMat(cellID(iC),trials,data,ops);
                    corrMatS2(iC,:,:)=tmpC;
                    shiftMatS2(iC,:,:)=tmpS;
                    
                end
                
            end
            
            corrMatS2=corrMatS2(all_good,:,:);
            shiftMatS2=shiftMatS2(all_good,:,:);
            
            
            data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);
            data_out.corrMat = corrMat;
            data_out.shiftMat = shiftMat;
            
            data_out.corrMatShifted = corrMatShifted;
            data_out.shiftMatShifted = shiftMatShifted;
            data_out.shiftMatS2 = shiftMatS2;
            data_out.corrMatS2 = corrMatS2;
            
            data_out.region = reg;
            
            
            data_out.factors = factors;
            data_out.trials = trials;
        end
    catch ME
        rethrow(ME)
        disp(ME.message)
        disp(sprintf('filenr: %d',iF))
    end
end

