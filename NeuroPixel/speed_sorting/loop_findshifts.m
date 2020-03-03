

ops.factors = -.55:0.01:.55;
ops.BinWidth =2;
ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.n_preceeding = 18;
ops.TimeBin = 0.02;
ops.idx = [10:ops.BinWidth:390]/ops.BinWidth;% in bins
fi = gausswin(5);
fi=fi'/sum(fi);
ops.filter = fi;
ops.n_trials = 10;
ops.plotfig = false;
ops.maxLag = 20; % in cm
OAK='/oak/stanford/groups/giocomo/';
%% savedir =
savedir = fullfile(OAK,'attialex','speed_filtered_5');
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

% data_dir=fullfile(OAK,'attialex','NP_DATA');
% %session_name = {'AA5_190809_gain_1'};
% filenames = {};
% sn = dir(fullfile(data_dir,'*.mat'));
% for iS = 1:numel(sn)
%     if ~(contains(sn(iS).name,'mismatch') || contains(sn(iS).name,'playback') || contains(sn(iS).name,'dark'))
%         filenames{end+1}=sn(iS).name(1:end-4);
%     end
% end

gain = 0.8;
contrast = 100;
region = 'VISp';
[filenames,triggers] = getFilesCriteria(region,contrast,gain,'F:/NP_DATA');

%%
p=gcp('nocreate');
if isempty(p)
    parpool(12);
end

%%


parfor iF=1:numel(filenames)
    try
        data = load(filenames{iF});
        ops_temp = ops;
        bl_trials = 1:numel(data.trial_contrast);
        bl_trials = bl_trials(data.trial_contrast==100 & data.trial_gain==1);
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
        %ops_here.trial = find(data.trial_gain ==1 & data.trial_contrast==100);
        for iRep=1:3%numel(good_starts)
            ops_temp.trials = good_starts(iRep)+[0:ops_temp.n_trials-1];
            [data_out,fighandles] = findBestShifts(data,ops_temp);
            [~,mi]=max(data_out.all_stability,[],2);
            factors = ops_temp.factors(mi);
            all_factors(iRep,:)=factors;
        end
        
        [~,session_name]=fileparts(filenames{iF});
        mf =matfile(fullfile(savedir,session_name));
        mf.stability = data_out.all_stability;
        mf.region = data_out.region;
        mf.subregion = data_out.sub_reg;
        mf.depth = data_out.depth;
        mf.CID = data_out.CID;
        mf.start_idx = good_start;
        mf.all_factors = all_factors;
    catch ME
        fprintf('%s \nFailed for %s: %d \n',ME.message,filenames{iF},iF)
    end
end
