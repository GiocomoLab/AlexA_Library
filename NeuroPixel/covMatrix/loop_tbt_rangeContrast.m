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
OAK='/oak/stanford/groups/giocomo/';
OAK = '/Volumes/Samsung_T5';
gains = [0.5, 0.6, 0.7, 0.8];
%%
gain = [];
contrast = 50;
regions = {'MEC','VISp','RS'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA_corrected'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end
savepath = fullfile(OAK,'attialex','tbtxcorr_rangeContrast');

if ~isfolder(savepath)
    mkdir(savepath)
end

save(fullfile(savepath,'parameters.mat'),'ops');
%%
% p = gcp('nocreate');
% if isempty(p)
%     p = parpool(12);
% end
%%

%%
for iF=1:numel(filenames)
   
        
        
        data = load(filenames{iF});
        [~,sn]=fileparts(filenames{iF});
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
        trials=1:max(data.trial);
        ops_here = ops;
        ops_here.trials = trials;
        cellID = data.sp.cids(data.sp.cgs==2);
        
        [corrMat,~,~]=trialCorrMat(cellID,trials,data,ops);
        
            
           
         
        
        
        data_out = matfile(fullfile(savepath,sprintf('%s',sn)),'Writable',true);
        data_out.corrMat = corrMat;
        
        
        
        data_out.region = reg;
        
        
        
        data_out.trials = trials;
        data_out.gain = data.trial_gain;
        data_out.trial_contrast = data.trial_contrast;
        
    
end

