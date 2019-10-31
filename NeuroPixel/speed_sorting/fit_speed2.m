function fit_speed2(filepath,savepath,params,im_save_root)
[~,session_name,~]=fileparts(filepath);

if ~isempty(im_save_root)
    image_save_dir = fullfile(im_save_root,session_name);
    if ~isfolder(image_save_dir)
        mkdir(image_save_dir)
    end
    plot_data=true;
else
    plot_data=false;
end

%load file
data=load(filepath);
fn=fieldnames(data);
for iF=1:numel(fn);eval(sprintf('%s = data.%s;',fn{iF},fn{iF}));end

%process anatomy
if ~isvarname('anatomy')
    return
end
try
    if isfield(anatomy,'parent_shifted')
        region = anatomy.parent_shifted;
    else
        region = anatomy.cluster_parent;
    end
catch
    disp('no anatomy')
    return
end

%calculate spatial maps

trials=[1:max(trial)];
%trials = trials(trial_gain == 1 & trial_contrast == 100);
spatialMap=[];
dwell_time=[];
edges=[0:2:400];
edges(1)=-.01;
posx(posx<0)=0;
posx(posx>400)=400;
for iT=1:length(trials)
    idxVR=trial==trials(iT);
    t_time=post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=sp.st<stop & sp.st>=start;
    [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',edges,'max_clust',max(sp.clu)+1);
    spatialMap=cat(3,spatialMap,spM);
    dwell_time=cat(1,dwell_time,dT);
end
good_cells = sp.cids(sp.cgs==2);

spatialMap=spatialMap(good_cells+1,:,:);

%normalize by dwell time in each bin
dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end

spatialMap(isnan(spatialMap))=0;

%split and flatten maps by region
[regs,~,idx] = unique(region(sp.cgs==2));
flatMaps = cell(numel(regs),1);
nTrials = size(spatialMap,3);

%define window for xcorr
posWindow=[60 140];
posBin = [find(edges == posWindow(1)), find(edges == posWindow(2))];
nBins=posBin(2)-posBin(1)+1;

filt = gausswin(5);
filt = filt/sum(filt);
for iR=1:numel(regs)
    nC=nnz(idx==iR);
    cellList = find(idx==iR);
    if nC<2 %disregard clusters with fewer than 10 cells
        continue
    end
    tmp_flat = zeros(nTrials,nBins*nC);
    for iT = 1:nTrials
        for iC=1:nC
            currentCell = cellList(iC);
            t_idx = (iC-1)*nBins+1:iC*nBins;
            tmp = squeeze(spatialMap(currentCell,:,iT));
            tmp = conv(tmp,filt,'same');
            tmp = tmp(posBin(1):posBin(2));
            tmp_flat(iT,t_idx)=tmp;
        end
    end
    flatMaps{iR}=tmp_flat;
end

%calculate speed difference on a trial by trial basis
speed = calcSpeed(posx,params);
trial_speed = zeros(1,max(trial));
for ii = 1:max(trial)
    idx = posx>posWindow(1) & posx<posWindow(2) & trial == ii;
    tmp = mean(speed(idx));
    
    trial_speed(ii)=tmp;
end

delays = cell(numel(regs),1);
for iR = 1:numel(flatMaps)
    if isempty(flatMaps{iR})
        continue
    else
        delays{iR}=calcXCorr(flatMaps{iR},trial_speed,trial_gain,trial_contrast);
    end
end

data.delays = delays;
data.session = session_name;
data.regions = regs;
data.clu_region = region;
data.posWindow = posWindow;

save(fullfile(savepath,session_name),'data')

if plot_data
    corr = nan(1,numel(delays));
    for iR=1:numel(delays)
        if ~isempty(delays{iR})
            
            bl_idx = delays{iR}(:,3)==1 & delays{iR}(:,5)==1 & delays{iR}(:,4)==100 & delays{iR}(:,6)==100;
            tmp = corrcoef(delays{iR}(bl_idx,1),delays{iR}(bl_idx,2));
            corr(iR)=tmp(1,2);
        end
    end
end
thresh = 0.08;
if plot_data && any(abs(corr)>thresh)
    [~,sid]=sort(trial_speed);
    
    trial_sorted =trial;
    trial_rank = 1:max(trial);
    trial_rank(sid)=trial_rank;
    for ii=1:max(trial)
        idx = trial==ii;
        trial_sorted(idx)=trial_rank(ii);
    end
    data.region = region;
    data.posWindow = posWindow;
    data.trial_sorted = trial_sorted;
    plotRasterSpeedSort(data,params,image_save_dir)
end

end

function delay_list=calcXCorr(flatMap,trial_speed,trial_gain,trial_contrast)
nTrials = size(flatMap,1);
delay_list=zeros(nTrials*(nTrials-1)/2,6); %delay, speed_diff, %i gain, i contrast, jgain j contrast
idx = 0;
max_lag=20;
for iT = 1:nTrials
    for jT=(iT+1):nTrials
        idx = idx+1;
        rr=xcorr(flatMap(iT,:),flatMap(jT,:),max_lag);
        [~,midx]=max(rr);
        if ismember(midx,[1 size(rr,2)])
            midx = nan;
        end
        
        
        tmp_d = max_lag+1;
        delay_list(idx,:)=[midx-tmp_d,trial_speed(iT)-trial_speed(jT), trial_gain(iT),trial_contrast(iT),trial_gain(jT),trial_contrast(jT) ];
    end
end
end