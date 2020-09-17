
ops = load_default_opt;
ops.cm_ops = load_default_opt; %this one will get passed to the calcCorrMat (with 2cm spatial bin)
ops.factors = -.3:0.01:.3;

ops.SpatialBin = 1; 
ops.idx = [10:ops.SpatialBin:390]/ops.SpatialBin;

ops.n_trials = 4;
ops.plotfig = false;
OAK='Z:\giocomo';
%OAK='/oak/stanford/groups/giocomo/';
OAK = '/Volumes/T7/';
%% savedir =
savedir = fullfile(OAK,'attialex','speed_filtered_05');
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



gain = 0.5;
contrast = 100;
%region = 'VISp';
regions = {'MEC'};
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
    parpool();
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
        not_done = true;
        start=1;
        good_starts = [];
        for iT=1:numel(triggers{iF})
            good_starts(end+1)=triggers{iF}(iT)-4;
            good_starts(end+1)=triggers{iF}(iT);
            good_starts(end+1)=triggers{iF}(iT)+4;
        end


        
        nC=nnz(data.sp.cgs==2);
        all_factors = nan(numel(good_starts),nC);
        all_stability = all_factors;
        all_firingRate=all_factors;
        %ops_here.trial = find(data.trial_gain ==1 & data.trial_contrast==100);
        %run shift finding for each block
        for iRep=1:numel(good_starts)
            ops_temp.trials = good_starts(iRep)+[0:ops_temp.n_trials-1];
            
            [data_out] = findBestShifts_mgc(data,ops_temp);
            [a,mi]=max(data_out.all_stability,[],2);
            factors = ops_temp.factors(mi);
            factors(isnan(a))=nan;
            all_factors(iRep,:)=factors;
            all_stability(iRep,:)=data_out.stability;
            
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
