ops = load_default_opt;
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
ops.num_tr_tot = 10;
ops.trial_range = [0:9];

%OAK = '/Volumes/Crucial X8/attialex';
%OAK = '/Volumes/Samsung_T5/attialex';
%OAK = 'F:\attialex';
%%
gain = 1.0;
contrast = 100;
regions = {'VISp','RS','MEC'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'NP_DATA_corrected'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end
savepath = fullfile(OAK,'slidingWindow_baseline_100cm');
if ~isfolder(savepath)
    mkdir(savepath)
end

%%
p = gcp('nocreate');
if isempty(p)
    p = parpool();
end
% filenames={'/Volumes/GoogleDrive/My Drive/alex_data/AA44_190925_gain_1.mat'};
% triggers{1}=[7];
%%
for iF=1:numel(filenames)
    
    try
        [~,sn]=fileparts(filenames{iF});
                    data = load(filenames{iF});
                    
                    
         nrectrials = max(data.trial);
        %find blocks of n_trials of baseline
        bl_trials = 1:nrectrials;
        bl_trials = bl_trials(data.trial_contrast(1:nrectrials)==100 & data.trial_gain(1:nrectrials)==1);
        not_done = true;
        start=1;
        good_starts = [];
        while not_done
            t_vec = start:(start+ops.num_tr_tot-1);
            if all(ismember(t_vec,bl_trials))
                good_starts(end+1)=start;
                start=start+ops.num_tr_tot;
            else
                start=start+1;
            end
            if start>max(bl_trials)
                not_done = false;
            end
            
        end

        for iRep=1:numel(good_starts)
            if isfile(fullfile(savepath,sprintf('%s_%d.mat',sn,iRep)))
                continue
            end

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
            
            ops_temp = ops;
        
         
            %calculate trial by trial correlation across all bins
            trials=good_starts(iRep)+ops.trial_range;
            ops_here = ops;
            ops_here.trials = trials;
            cellID = data.sp.cids(data.sp.cgs==2);
            
            [PEAKS,SHIFTS,corrMat,shiftMat,SPEED,SPEED2]=calculatePeakShiftSession_new(data,trials,ops);

      
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
            data_out.SPEED = SPEED;
            data_out.SPEED2 = SPEED2;
            data_out.trials = trials;
        end
    catch ME
        disp(ME.message)
        disp(sprintf('filenr: %d',iF))
        %rethrow(ME)
    end
end

