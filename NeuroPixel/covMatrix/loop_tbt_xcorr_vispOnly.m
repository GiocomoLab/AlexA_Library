%new version that uses malcolms approach for shift correction
ops = load_default_opt;
ops.trial_range = [-6:9];
ops.BinWidth = 2;
ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.smoothSigma=ops.smoothSigma_dist;
smoothSigma = ops.smoothSigma/ops.BinWidth;
ops.filter = gausswin(floor(smoothSigma*5/2)*2+1);
ops.filter = ops.filter/sum(ops.filter);
ops.max_lag = 30;
ops.maxLag = ops.max_lag;
%OAK='/oak/stanford/groups/giocomo/';
OAK = '/Volumes/T7';
%OAK = '/Volumes/Crucial X8/';

rec_field_data = load('/Users/attialex/code/AlexA_Library/NeuroPixel/receptive_fields/receptive_field_table.mat');
rec_field_names = unique(rec_field_data.cell_table.spatial_name);
%%
gain = 0.8;
contrast = 100;
regions = {'VISp'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA_corrected'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end
[~,basenames]=fileparts(filenames);
idx = ismember(basenames,rec_field_names);
filenames = filenames(idx);
triggers = triggers(idx);

savepath = fullfile(OAK,'attialex','tbtxcorr_visp');
shiftDir = fullfile(OAK,'attialex','speed_filtered_correctedData_shortidx2');
if ~isfolder(savepath)
    mkdir(savepath)
end
shift_ops = load(fullfile(shiftDir,'parameters.mat'));
shift_ops = shift_ops.ops;
save(fullfile(OAK,'attialex','parameters.mat'),'ops');
%%
% p = gcp('nocreate');
% if isempty(p)
%     p = parpool();
% end
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
            reg_orig = data.anatomy.cluster_parent((data.sp.cgs==2));
            if iscolumn(reg_orig)
                reg_orig = reg_orig';
            end
            
            
            
            %calculate trial by trial correlation across all bins
            trials=triggers{iF}(iRep)+ops.trial_range;
            ops_here = ops;
            ops_here.trials = trials;
            cellID = data.sp.cids(data.sp.cgs==2);
            
            [corrMat,fr_map,shiftMat]=trialCorrMat(cellID,trials,data,ops);
            trial_idx = ismember(data.trial,trials);
            speed = calcSpeed(data.posx(trial_idx),ops);
            
            if ~all(data.trial_contrast(trials)==contrast)
                %error('gain trials violating contrast condition')
                disp('gain trials violating contrast condition')
                continue
            end
            
            
            data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);
            data_out.corrMat = corrMat;
            data_out.shiftMat = shiftMat;
            data_out.trials = trials;
            data_out.region = reg;
            data_out.region_orig = reg_orig;
            data_out.speed = speed;
            data_out.trials = data.trial(trial_idx);
            data_out.good_cells = cellID;
            data_out.baseline_map = squeeze(nanmean(fr_map(:,1:6,:),2));

            
            
            %load shift factors
            if ~isfile(fullfile(shiftDir,[sn '.mat']))
                sprintf('no speed file for %s',sn)
                continue
            end
            shift_data = load(fullfile(shiftDir,[sn '.mat']));
            tmp = shift_data.all_stability;
            stab_idx = tmp>.5;
            tmp_f = shift_data.all_factors;
            tmp_f(~stab_idx)=nan;
            enough_idx = sum(~isnan(tmp_f))>2;
            factors = nanmean(tmp_f);
            factors(~enough_idx)=nan;
            
            
           
            %correct spike times for each unit
            dat2=data;
            for iC=1:numel(cellID)
                if ~isnan(factors(iC))
                idx = dat2.sp.clu==cellID(iC);
                dat2.sp.st(idx)=dat2.sp.st(idx)+factors(iC); %%verify
                end
            end
            [corrMatS,~,shiftMatS]=trialCorrMat(cellID,trials,dat2,ops);

            
            %set to nan where factor was nan
            for iC = 1:numel(cellID)
                if isnan(factors(iC))
                    corrMatS(iC,:,:)=nan;
                    shiftMatS(iC,:,:)=nan;  
                end
                
            end
            
            
            
            
            
            data_out.corrMatShifted = corrMatS;
            data_out.shiftMatShifted = shiftMatS;
            data_out.corrMatS2 = corrMatS;
            data_out.shiftMatS2= shiftMatS;  
            data_out.factors = factors;
        end
    catch ME
        disp(ME.message)
        disp(sprintf('filenr: %d',iF))
        rethrow(ME)
    end
end

