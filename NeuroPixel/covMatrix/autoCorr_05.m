%% parameters
ops = load_default_opt;
ops.dark = true; %fore it to use max pos bin
ops.max_lag_auto = 100;
ops.trial_range = [-6:9];

OAK='/oak/stanford/groups/giocomo/';
%OAK = '/Volumes/Samsung_T5';

savepath = fullfile(OAK,'attialex','autocorr');
if ~isfolder(savepath)
    mkdir(savepath)
end
dark_folder = fullfile(OAK,'attialex','distance_tuning_xcorr_only');

%% 0.5 files mec
gain = 0.5;
contrast = 100;
regions = {'MEC'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA_corrected'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end

%%
p = gcp('nocreate');
if isempty(p)
    p = parpool(12);
end
%% main loop



parfor iF=1:numel(filenames)
    try
    [~,sn]=fileparts(filenames{iF});
    
    parts = strsplit(sn,'_');
    session_part = strcat(parts{1},'_',parts{2});
    
    poss=dir(fullfile(dark_folder,strcat(session_part,'*')));
    data = load(filenames{iF});
    nC=nnz(data.sp.cgs==2);
    if isempty(poss)
        disp('no dark data')
        is_tuned = nan(nC,1);
        continue
    else
        dark_data = load(fullfile(dark_folder,poss(1).name));
        
        is_tuned = zeros(nC,1);
        peak_loc_pf = zeros(nC,1);
        MAX_RANGE_BIN = 300;
        upper_limits = squeeze(dark_data.acg_upper_limits(:,:,end));
        for iN = 1:nC
            
            
            a=strfind(diff(dark_data.ACG(iN,:))>0,[0 1]); %find elbow point/location of first minimum
            if numel(a)<1 %no elbow point detected, skip
                continue
            end
            
            if a(1)>MAX_RANGE_BIN % elbow point beyond 500cm, probably noise, skip
                continue
            end
            [pks,loc]=findpeaks(dark_data.ACG(iN,1:MAX_RANGE_BIN),'SortStr','descend','NPeaks',1,'MinPeakProminence',.2);
            if upper_limits(iN,loc)<pks
                peak_loc_pf(iN)=loc;
            end
            
            sub = dark_data.ACG(iN,a(1):MAX_RANGE_BIN); %only search within between elbow point and MAX_RANGE_BIN
            [ma,mi]=max(sub);
            mi = mi+a(1)-1;
            is_tuned(iN) = upper_limits(iN,mi)<ma;
            %         mod_depth(iN)=ma-ACG(iN,a(1));
            %         peak_loc(iN)=mi;
            
            
        end
    end
    
    
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
    
    cellID = data.sp.cids(data.sp.cgs==2);
    
    for iRep=1:numel(triggers{iF})
        trials=triggers{iF}(iRep)+ops.trial_range;
        bl_trials = trials(3:6);
        gain_trials = trials(7:10);
        %         ops_here = ops;
        %         ops_here.trials = trials;
        
        [corrMat,fr_map,shiftMat]=trialCorrMat(cellID,trials,data,ops);
        baseline_fr = calcFRVsDist(cellID,bl_trials,data,ops);
        gain_fr = calcFRVsDist(cellID,gain_trials,data,ops);
        gain_fr_corrected = calcFRVsDist(cellID,gain_trials,data,ops,gain); %use actual distance run on track
        baseline_fr = baseline_fr - mean(baseline_fr,2);
        gain_fr = gain_fr - mean(gain_fr,2);
        gain_fr_corrected = gain_fr_corrected - mean(gain_fr_corrected,2);
        xc_bl = zeros(nC,2*ops.max_lag_auto+1);
        xc_gain = xc_bl;
        xc_gain_corrected = xc_bl;
        for iC=1:nC
            xc_bl(iC,:)=xcorr(baseline_fr(iC,:),'coeff',ops.max_lag_auto);
            xc_gain(iC,:)=xcorr(gain_fr(iC,:),'coeff',ops.max_lag_auto);
            xc_gain_corrected(iC,:)=xcorr(gain_fr_corrected(iC,:),'coeff',ops.max_lag_auto);
        end
        
        data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);
        data_out.corrMat = corrMat;
        data_out.xc_bl = xc_bl;
        data_out.xc_gain = xc_gain;
        data_out.xc_gain_corrected = xc_gain_corrected;
        data_out.is_tuned = is_tuned;
        data_out.region = reg;
    end
    catch ME
        disp(filenames{iF})
    end
end