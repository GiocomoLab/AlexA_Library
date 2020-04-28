

ops.factors = -.3:0.01:.3;
ops.BinWidth =1;
ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.TimeBin = 0.02;
ops.idx = [10:ops.BinWidth:390]/ops.BinWidth;% in bins
fi = gausswin(22);
fi=fi'/sum(fi);
ops.filter = fi;
spfi = gausswin(5);
spfi = spfi/sum(spfi);
ops.speed_filter = spfi; 
ops.n_trials = 10;
ops.plotfig = false;
ops.maxLag = 20; % in cm
OAK='/oak/stanford/groups/giocomo/';
%% savedir =
savedir = fullfile(OAK,'attialex','speed_filtered_correctedData');
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



gain = 0.8;
contrast = 100;
region = 'VISp';
regions = {'VISp','RS','MEC'};
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
    parpool(12);
end

%%
zero_idx = find(ops.factors==0);

parfor iF=1:numel(filenames)
    [~,session_name]=fileparts(filenames{iF});

    if isfile(fullfile(savedir,[session_name '.mat']))
        disp('exisits')
        continue
    end
    try
        data = load(filenames{iF});
        ops_temp = ops;
        nrectrials = max(data.trial);
        %find blocks of n_trials of baseline
        bl_trials = 1:nrectrials;
        bl_trials = bl_trials(data.trial_contrast(1:nrectrials)==100 & data.trial_gain(1:nrectrials)==1);
        not_done = true;
        start=1;
        good_starts = [];
        while not_done
            t_vec = start:(start+ops.n_trials-1);
            if all(ismember(t_vec,bl_trials))
                good_starts(end+1)=start;
                start=start+ops.n_trials;
            else
                start=start+1;
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
            [data_out,~] = findBestShifts(data,ops_temp);
            [~,mi]=max(data_out.all_stability,[],2);
            factors = ops_temp.factors(mi);
            all_factors(iRep,:)=factors;
            stability = data_out.all_stability(:,zero_idx);
            all_stability(iRep,:)=stability;
            all_firingRate(iRep,:)=data_out.firing_rate;
        end
        
        mf =matfile(fullfile(savedir,session_name),'Writable',true);
        mf.region = data_out.region;
        mf.subregion = data_out.sub_reg;
        mf.depth = data_out.depth;
        mf.CID = data_out.CID;
        mf.start_idx = good_starts;
        mf.all_factors = all_factors;
        mf.all_stability = all_stability;
        mf.FR = all_firingRate;
    catch ME
        %rethrow(ME)
        disp(iF)
        fprintf('%s \nFailed for %s: %d \n',ME.message,filenames{iF},iF)
        rethrow(ME)
    end
end
