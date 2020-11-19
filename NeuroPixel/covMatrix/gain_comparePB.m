fi = dir('/Volumes/T7/attialex/NP_DATA_corrected/AA*back*Towers*.mat');


%%
gain_names={};
pb_names = {};
for iF = 1:numel(fi)
    
    data=load(fullfile(fi(iF).folder,fi(iF).name),'anatomy');
    if numel(fields(data)) ==0
        continue
    end
    idx = strfind(fi(iF).name,'_');
    newStr = [fi(iF).name(1:idx(2)), '*.mat'];
    fa = dir(fullfile(fi(iF).folder,newStr))
    pp = {fa(:).name}
    pot = ~contains(pp,{'mismatch','play','dark'})
    nPot = nnz(pot)
    if nPot>0
        idx = find(nPot,1);
        gain_names{end+1}=pp{idx};
        pb_names{end+1}=fi(iF).name;
    end
    
    
end

%%

for iF=1:numel(gain_names)
    figure
    data = load(fullfile(fi(1).folder,gain_names{iF}));
    plot(data.trial_gain)
    hold on
    plot(data.trial_contrast/100)
end

%%
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
OAK = '/Volumes/T7';
GAIN = 0.5;
savepath = fullfile(OAK,'attialex',['tbtxcorr_', num2str(GAIN), '_withPB']);
if ~isfolder(savepath)
    mkdir(savepath)
end

%%

for iF=1:numel(gain_names)
    data = load(fullfile(fi(1).folder,gain_names{iF}));
    triggers = strfind(data.trial_gain' == GAIN & data.trial_contrast'==100,[0 1])+1;
    data_pb = load(fullfile(fi(1).folder,pb_names{iF}));
    cellID = data.sp.cids(data.sp.cgs==2);
            
    [corrMatPB,fr_map,shiftMat]=trialCorrMat(cellID,1:max(data_pb.trial)-1,data_pb,ops);
    [~,sn]=fileparts(gain_names{iF});
    for iRep=1:numel(triggers)
            
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
            trials=triggers(iRep)+ops.trial_range;
            ops_here = ops;
            ops_here.trials = trials;
            cellID = data.sp.cids(data.sp.cgs==2);
            
            [corrMat,fr_map,shiftMat]=trialCorrMat(cellID,trials,data,ops);
            [corrMatBL,fr_map,shiftMat]=trialCorrMat(cellID,1:20,data,ops);
            
            trial_idx = ismember(data.trial,trials);
            
           
            
            data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);
            data_out.corrMat = corrMat;
            data_out.corrMatBL = corrMatBL;
            data_out.corrMatPB = corrMatPB;
            data_out.shiftMat = shiftMat;
            data_out.trials = trials;
            data_out.region = reg;
            data_out.region_orig = reg_orig;
            data_out.trials = data.trial(trial_idx);

            data_out.baseline_map = squeeze(nanmean(fr_map(:,1:6,:),2));
    end
end

%% next, for each stable neuron, plot shift vs spatial stability pb