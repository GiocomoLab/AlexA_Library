
%%
ops.region = 'MEC';
ops.binsize=5;
ops.smoothSigma = 4;
ops.maxlag = 10;
[filenames,triggers] = getFilesCriteria(ops.region,100,0.5,'/oak/stanford/groups/giocomo/attialex/NP_DATA');


%%

%load
for iF = 1:numel(filenames)
    clear data
    if numel(triggers{iF})==1
        continue
    end
    
    
    
    data = load(filenames{iF});
    
    
    %calculate spatial firing rate
    spatialMap=[];
    dwell_time=[];
    edges=[0:ops.binsize:400];
    edges(1)=-.01;
    data.posx(data.posx<0)=0;
    data.posx(data.posx>=400)=399.00;
    trials = 1:max(data.trial);
    for iT=1:length(trials)
        idxVR=data.trial==trials(iT);
        t_time=data.post(idxVR);
        start=min(t_time);
        stop=max(t_time);
        idxNP=data.sp.st<stop & data.sp.st>=start;
        [spM, dT]=getSpikeMatPosition(data.sp.st(idxNP),data.sp.clu(idxNP),data.posx(idxVR),data.post(idxVR),'edges',edges,'max_clust',max(data.sp.clu)+1);
        spatialMap=cat(3,spatialMap,spM);
        dwell_time=cat(1,dwell_time,dT);
    end
    %cellIDX=find(sp.cgs>=1);
    if isfield(data.anatomy,'parent_shifted')
        reg = startsWith(data.anatomy.parent_shifted,ops.region);
    else
        reg = startsWith(data.anatomy.cluster_parent,ops.region);
    end
    if iscolumn(reg)
        reg = reg';
    end
    spatialMap=spatialMap(data.sp.cids+1,:,:);
    
    
    dt=dwell_time';
    dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
    for ii=1:size(spatialMap,1)
        spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
    end
    spatialMap(isnan(spatialMap))=0;
    smoothSigma = ops.smoothSigma/ops.binsize;
    smoothWindow = floor(smoothSigma*5/2)*2+1;
    gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
    filt = reshape(gauss_filter,[1, numel(gauss_filter),1]);
    sPF = repmat(spatialMap,[1,3,1]);
    sPF=convn(sPF,filt,'same');
    iidx = (size(spatialMap,2)+1):(2*size(spatialMap,2));
    sPF = sPF(:,iidx,:);
    % get template trials
    
    spatialMap = sPF;
    
    
    %pp=bsxfun(@rdivide,spatialMap,mm);
    fr_mat = squeeze(nanmean(spatialMap));
    spatialMap=spatialMap(data.sp.cgs==2 & reg,:,:);
    %%
    onsets = triggers{iF}-1; % bc rest of code assumes trigger is last bl trial
    n_pre = 4;
    %calculate map of 4 preceeding trials for first rep
    bl_pre = mean(spatialMap(:,:,onsets(1)-3:onsets(1)),3);
    gain_map = mean(spatialMap(:,:,onsets(2)+1:onsets(2)+4),3);
    
    % for each of the gc trials, similarity (peak xcorr) to bl_pre map
    nCells = size(spatialMap,1);
    xcorrs_bl=nan(nCells,4,2*ops.maxlag+1);
    for iT=1:4
        for iC = 1:size(spatialMap,1)
            tmp_template = bl_pre(iC,:);
            tmp_template = tmp_template-mean(tmp_template);
            tmp_resp = squeeze(spatialMap(iC,:,onsets(1)+iT));
            tmp_resp = tmp_resp-mean(tmp_resp);
            [aa,lags] = xcorr(tmp_template,tmp_resp,ops.maxlag,'coeff');
            xcorrs_bl(iC,iT,:)=aa;
        end
    end
    % similarity to all of .8 2nd rep
    xcorrs_gain=nan(nCells,4,2*ops.maxlag+1);
    for iT=1:4
        for iC = 1:size(spatialMap,1)
            tmp_template = gain_map(iC,:);
            tmp_template = tmp_template-mean(tmp_template);
            tmp_resp = squeeze(spatialMap(iC,:,onsets(1)+iT));
            tmp_resp = tmp_resp-mean(tmp_resp);
            [aa,lags] = xcorr(tmp_template,tmp_resp,ops.maxlag,'coeff');
            xcorrs_gain(iC,iT,:)=aa;
        end
    end
    
    % similarity to 4st 2nd third 4th
    xcorrs_gainAll=nan(4,nCells,4,2*ops.maxlag+1);
    for iTemplate = 1:4
        for iT=1:4
            for iC = 1:size(spatialMap,1)
                tmp_template = squeeze(spatialMap(iC,:,onsets(2)+iTemplate));
                tmp_template = tmp_template-mean(tmp_template);
                tmp_resp = squeeze(spatialMap(iC,:,onsets(1)+iT));
                tmp_resp = tmp_resp-mean(tmp_resp);
                [aa,lags] = xcorr(tmp_template,tmp_resp,ops.maxlag,'coeff');
                xcorrs_gainAll(iTemplate,iC,iT,:)=aa;
            end
        end
    end
    
    %%
    bl_similarity = max(xcorrs_bl,[],3);
    average_gain_similarity = max(xcorrs_gain,[],3);
    trial_gain_similarity={};
    for ii=1:4
        tmp = squeeze(xcorrs_gainAll(ii,:,:,:));
        trial_gain_similarity{ii}=max(tmp,[],3);
    end
    %%
    template_trials = onsets(1)-3:onsets(1);
    idx = find(triu(true(numel(template_trials)),1));
    
    stability = zeros(1,size(spatialMap,1));
    for ii=1:numel(stability)
        tmp = squeeze(spatialMap(ii,:,template_trials));
        tmp = corr(tmp);
        stability(ii)=mean(tmp(idx));
    end
    
    [~,session_name,~]=fileparts(filenames{iF});
    save(sprintf('/oak/stanford/groups/giocomo/attialex/gainmap_comparison/%s',session_name),'bl_similarity','average_gain_similarity','trial_gain_similarity','stability')
    
end