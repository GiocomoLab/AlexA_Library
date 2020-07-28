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
ops.chunksize=100; %in bins,so thats 200 cm
ops.stride_start = 1;%10;
ops.stride = 5;
OAK='/oak/stanford/groups/giocomo/attialex';
OAK = '/Volumes/Samsung_T5/attialex';
%%
gain = 0.5;
contrast = 100;
regions = {'VISp','RS'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'NP_DATA_corrected'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end
savepath = fullfile(OAK,'tbtxcorr_05_reg_200cmChunk');
shiftDir = fullfile(OAK,'attialex','speed_filtered_new_22binspace_5binspeed2');
if ~isfolder(savepath)
    mkdir(savepath)
end

%%
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
            if ~isfield(data.anatomy,'depth')
                depth = data.anatomy.tip_distance(data.sp.cgs==2);
                mec_entry = data.anatomy.z2;
            else
                depth = data.anatomy.depth;
                mec_entry = nan;
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
            
            [PEAKS,SHIFTS,corrMat,shiftMat]=calculatePeakShiftSession_new(data,trials,ops);
            if ~all(data.trial_contrast(trials)==contrast)
                error('gain trials violating contrast condition')
                
            end
            if ~all(data.trial_gain((0:3)+triggers{iF}(iRep))==gain)
                error('gain trials wrong gain')
                
            end
            
            
            
           
            
            data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);
            data_out.corrMat = corrMat;
            data_out.shiftMat = shiftMat;
            data_out.PEAKS = PEAKS;
            
            data_out.SHIFTS = SHIFTS;
            
            data_out.region = reg;
            data_out.depth = depth;
            data_out.mec_entry = mec_entry;
            data_out.good_Cells = cellID;
            data_out.region_orig = reg_orig;
            
            data_out.trials = trials;
        end
    catch ME
        disp(ME.message)
        disp(sprintf('filenr: %d',iF))
        rethrow(ME)
    end
end

