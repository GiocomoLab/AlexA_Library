function data_out=slow_vs_fastTrials_v2(data,image_save_dir,params)
%make it more consistent with xcorr for other shifts

if ~isempty(image_save_dir)
    if exist(image_save_dir,'dir')~=7
        mkdir(image_save_dir);
    end
    save_images = true;
else
    save_images=false;
end


if ~isfield(data,'anatomy')
    return
end
anatomy = data.anatomy;
trial = data.trial;
posx = data.posx;
post = data.post;
sp = data.sp;
trial_gain = data.trial_gain;
trial_contrast = data.trial_contrast;
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

trials=[1:max(trial)];
%trials = trials(trial_gain == 1 & trial_contrast == 100);
spatialMap=[];
dwell_time=[];
binsize=2;
edges=[0:binsize:400];
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
%cellIDX=find(sp.cgs>=1);
good_cells = sp.cids(sp.cgs==2);

spatialMap=spatialMap(good_cells+1,:,:);
%spatialMap=spatialMap(:,1:end-1,:);
%dwell_time=dwell_time(:,1:end-1);
%normalize by dwell time in each bin
dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
smoothSigma = 4/binsize;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
filt = reshape(gauss_filter,[1, numel(gauss_filter),1]);
sPF = repmat(spatialMap,[1,3,1]);
sPF=convn(sPF,filt,'same');
iidx = (size(spatialMap,2)+1):(2*size(spatialMap,2));
sPF = sPF(:,iidx,:);
% get template trials

spatialMap = sPF;



%% calc speed differential
speed = calcSpeed(posx,params);
posWindow= [10 390];
%posWindow= [120 200];
posBin = [find(edges == posWindow(1)), find(edges == posWindow(2))];
trial_speed = zeros(1,max(trial));
for ii = 1:max(trial)
    idx = posx>posWindow(1) & posx<posWindow(2) & trial == ii;
    tmp = mean(speed(idx));
    
    trial_speed(ii)=tmp;
end
[~,sid]=sort(trial_speed);

trial_sorted = trial;
trial_rank = 1:max(trial);
trial_rank(sid)=trial_rank;
for ii=1:max(trial)
    idx = trial==ii;
    trial_sorted(idx)=trial_rank(ii);
end

spMap=shiftdim(spatialMap,1);


bl_trials = find(trial_gain == 1 & trial_contrast == 100);


bl_sorted = sid;
bl_sorted(~ismember(bl_sorted,bl_trials))=[];
nT=numel(bl_sorted);
M=round(.1*nT);

slow_idx = bl_sorted(1:M);
fast_idx = bl_sorted(nT-M:nT);
%h=figure('visible','off');

if save_images
    fig=figure('visible','off');
end
delay_per_cell = zeros(2,numel(good_cells));
trials_fast = zeros(numel(edges)-1,numel(good_cells));
trials_slow = trials_fast;
clu_reg = cell(numel(good_cells),1);
XCORRS = zeros(numel(good_cells),21);
for cellIDX = 1:numel(good_cells)
    cluID = find(sp.cids==good_cells(cellIDX));
    clu_reg{cellIDX}=region{cluID};
    mS=spMap(:,:,cellIDX);
    mS(isnan(mS))=0;
    mS=conv2(mS,filt,'same');
    tmp_slow = nanmean(mS(:,slow_idx),2);
    tmp_fast = nanmean(mS(:,fast_idx),2);
    trials_slow(:,cellIDX)=tmp_slow;
    trials_fast(:,cellIDX)=tmp_fast;
    ts =tmp_slow(posBin(1):posBin(2));
    ts = ts-mean(ts);
    tf = tmp_fast(posBin(1):posBin(2));
    tf = tf-mean(tf);
    [rr,lags]=xcorr(ts,tf,10,'coeff');
    XCORRS(cellIDX,:)=rr;
    [~,tmp]=max(rr);
    delay_per_cell(1,cellIDX)=lags(tmp)*mean(diff(edges));
    delay_per_cell(2,cellIDX)=rr(tmp);
    if save_images
        
        ax=subplot(1,1,1);
        p1=plot(edges(2:end),(tmp_slow));
        hold(ax,'on')
        p2=plot(edges(2:end),tmp_fast);
        title(round(delay_per_cell(1,cellIDX)))
        
        
        for ij=posWindow
            xline(ij,'r');
        end
        legend([p1 p2],{'slow','fast'})
        
        saveas(fig,fullfile(image_save_dir,sprintf('%s_%s_%d.png',clu_reg{cellIDX},session_name,cellIDX)),'png');
        clf
    end
    %%
end

data_out.delay = delay_per_cell;
data_out.XCORR = XCORRS;
data_out.region = clu_reg;
data_out.slow = mean(trial_speed(slow_idx));

data_out.fast = mean(trial_speed(fast_idx));
data_out.slow_trials = trials_slow;
data_out.fast_trials = trials_fast;

end