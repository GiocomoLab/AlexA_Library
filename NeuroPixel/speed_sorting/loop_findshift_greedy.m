

ops.factors = -.55:0.01:.55;
ops.BinWidth =1;
ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.TimeBin = 0.02;
ops.idx = [10:ops.BinWidth:390]/ops.BinWidth;% in bins
fi = gausswin(11);
fi=fi'/sum(fi);
ops.filter = fi;
spfi = gausswin(5)'/sum(gausswin(5));
ops.speed_filter = spfi; 
ops.n_trials = 10;
ops.plotfig = false;
ops.maxLag = 20; % in cm
OAK='/oak/stanford/groups/giocomo/';
%OAK='/Volumes/Samsung_T5/';
%% savedir =
savedir = fullfile(OAK,'attialex','speed_filtered_greedy3');
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
regions = {'VISp','MEC','RS'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
[tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA'));
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
        disp('\t \t exisits')
        continue
    end
    try
        data = load(filenames{iF});
        ops_temp = ops;
        
        good_chunks = {};
        chunk_contrast = [];
        current_chunk = 1;
        chunk_size = [];
        tmp = [];
        current_contrast = data.trial_contrast(1);
        current_gain = data.trial_gain(1);
        start_new = false;
        for iT = 1:numel(data.trial_gain)
            bl_gain = data.trial_gain(iT) ==1;
            old_contrast = current_contrast;
            current_contrast = data.trial_contrast(iT);

            contrast_switch = current_contrast ~= old_contrast;
            
            old_gain = current_gain;
            current_gain= data.trial_gain(iT);
            gain_switch = old_gain ~= current_gain;
            
            
            if ~contrast_switch && bl_gain && ~gain_switch % add to list
                tmp(end+1)=iT;
            elseif contrast_switch %end of chunk
                %sprintf('contrast_switch at %d',iT)
                good_chunks{current_chunk}=tmp;
                chunk_size(current_chunk)=numel(tmp);
                chunk_contrast(current_chunk)=old_contrast;
                start_new = true;
            elseif gain_switch && ~bl_gain
                %sprintf('gain_switch at %d',iT)
                
                    good_chunks{current_chunk}=tmp;
                    chunk_size(current_chunk)=numel(tmp);
                    chunk_contrast(current_chunk)=current_contrast;

            elseif gain_switch && bl_gain
                %sprintf('gain_switch at %d',iT)
                start_new = true;
            end
            if start_new
                tmp=[];
                tmp(end+1)=iT;
                current_chunk = current_chunk+1;
                start_new = false;
            end

        end
        keep = false(size(chunk_contrast));
        for iK=1:numel(keep)
            tmp = chunk_contrast(iK)==100 && numel(good_chunks{iK})>=10;
            keep(iK)=tmp;
        end
        good_chunks=good_chunks(keep);
        chunk_contrast = chunk_contrast(keep);
        
        nC=nnz(data.sp.cgs==2);
        all_factors = nan(numel(good_chunks),nC);
        all_stability = all_factors;
        all_firingRate = all_factors;
        %ops_here.trial = find(data.trial_gain ==1 & data.trial_contrast==100);
        for iRep=1:numel(good_chunks)
            ops_temp.trials = good_chunks{iRep};
            [data_out,fighandles] = findBestShifts(data,ops_temp);
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
        mf.start_idx = good_chunks;
        mf.chunk_contrast = chunk_contrast;
        mf.all_factors = all_factors;
        mf.all_stability = all_stability;
        mf.FR = all_firingRate;
    catch ME
        fprintf('%s \nFailed for %s: %d \n',ME.message,filenames{iF},iF)
    end
end
