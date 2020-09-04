
ops = load_default_opt;
ops.cm_ops = load_default_opt; %this one will get passed to the calcCorrMat (with 2cm spatial bin)
ops.factors = -.3:0.01:.3;

ops.SpatialBin = 1;
ops.idx = [10:ops.SpatialBin:390]/ops.SpatialBin;

ops.n_trials = 16;
ops.plotfig = false;
%OAK='Z:\giocomo';
OAK='/oak/stanford/groups/giocomo/';
%OAK = '/Volumes/T7/'
%pca_data=load('/Users/attialex/Downloads/rep_clusters.mat');

%% savedir =
savedir = fullfile(OAK,'attialex','speed_filtered_correctedData_shortidx3');
%savedir = fullfile('F:/temp/','speed_filtered');
imdir = fullfile(savedir,'images');
if ~isfolder(savedir)
    mkdir(savedir);
    
end
if ~isfolder(imdir)
    mkdir(imdir)
end
mf = matfile(fullfile(savedir,'parameters'),'Writable',true);
mf.ops = ops;


%% find files



gain = 1;
contrast = 100;
%region = 'VISp';
regions = {'VISp','RS','MEC'};
%regions={'CA'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA_corrected'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end



%%
p=gcp('nocreate');
if isempty(p)
    parpool(24);
end

%%
zero_idx = find(ops.factors==0);
fprintf('Progress:\n');
m=numel(filenames);
fprintf(['\n' repmat('.',1,m) '\n\n']);
parfor iF=1:numel(filenames)
    [~,session_name]=fileparts(filenames{iF});
    fprintf('\b|\n');
    if isfile(fullfile(savedir,[session_name '.mat']))
        %disp('exisits')
        continue
    end
    try
        data = load(filenames{iF});
        ops_temp = ops;
        nrectrials = max(data.trial);
        %find blocks of n_trials of baseline
        bl_trials = 1:nrectrials;
        bl_trials = bl_trials(data.trial_contrast(1:nrectrials)==100 & data.trial_gain(1:nrectrials)==1);
        bl=data.trial_contrast==100 & data.trial_gain==1;
        bl_onsets = strfind(bl',[0 1])+3; %+3 so that we skip first two trials after condition change to allow settle
        not_done = true;
        start=3; %skip the first two trials
        good_starts = [];
        
        
        if isfield(data.anatomy,'parent_shifted')
            reg = data.anatomy.parent_shifted;
        else
            reg = data.anatomy.cluster_parent;
        end
        
        if isfield(data.anatomy,'parent_shifted')
            sub_reg = data.anatomy.region_shifted;
        elseif isfield(data.anatomy,'cluster_region')
            
            sub_reg = data.anatomy.cluster_region;
        else
            sub_reg = {};
        end
        
        if iscolumn(reg)
            reg = reg';
            sub_reg = sub_reg';
        end
        region = reg;
        good_cells = data.sp.cids(data.sp.cgs==2 & startsWith(region,{'MEC','ECT','RS','VISp'}));
        
        while not_done
            t_vec = start:(start+ops.n_trials-1);
            if all(ismember(t_vec,bl_trials))
                good_starts(end+1)=start;
                start=start+ops.n_trials;
            else
                %start=start+1;
                %start next closest bl_onset
                bl_onsets(bl_onsets<=start)=[];
                if ~isempty(bl_onsets)
                    start = bl_onsets(1);
                else
                    start = inf;
                end
            end
            if start>max(bl_trials)
                not_done = false;
            end
            
        end
        nC=nnz(data.sp.cgs==2);
        all_factors = nan(numel(good_starts),nC);
        all_stability = all_factors;
        all_firingRate=all_factors;
        %ops_here.trial = find(data.trial_gain ==1 & data.trial_contrast==100);
        %run shift finding for each block
        for iRep=1:numel(good_starts)
            ops_temp.trials = good_starts(iRep)+[0:ops_temp.n_trials-1];
            
            %calculate trial_corr mat
            corrMat = trialCorrMat(good_cells,ops_temp.trials,data,ops_temp.cm_ops);
            bl_stab = nanmean(nanmean(corrMat(:,1:6,1:6),2),3);
            %cluster assignement
            corrMat = squeeze(nanmean(corrMat(bl_stab>=ops.stab_thresh,:,:),1));
            if nanmean(corrMat,'all')>.55
                clu=1;
            else
                clu=0;
            end
            
            if ~(nnz(bl_stab>=ops.stab_thresh)>=5) || clu ~= 1
                skip=true;
                sprintf('skipped nstab: %d, clu %d \n',nnz(bl_stab>=ops.stab_thresh),clu)
            else
                skip = false;
            end
            if ~skip
                [data_out] = findBestShifts_mgc(data,ops_temp);
                [~,mi]=max(data_out.all_stability,[],2);
                factors = ops_temp.factors(mi);
                all_factors(iRep,:)=factors;
                all_stability(iRep,:)=data_out.stability;
            end
            
        end
        
        mf =matfile(fullfile(savedir,session_name),'Writable',true);
        mf.region = data_out.region;
        mf.subregion = data_out.sub_reg;
        mf.depth = data_out.depth;
        mf.CID = data_out.CID;
        mf.start_idx = good_starts;
        mf.all_factors = all_factors;
        mf.all_stability = all_stability;
    catch ME
        %rethrow(ME)
        disp(iF)
        fprintf('%s \nFailed for %s: %d \n',ME.message,filenames{iF},iF)
        %rethrow(ME)
    end
end


