function [data_out,fighandles] = findBestShifts(data,ops)

%prepare variables
if isfield(data.anatomy,'parent_shifted')
    reg = data.anatomy.parent_shifted;
else
    reg = data.anatomy.cluster_parent;
end

if isfield(data.anatomy,'parent_shifted')
    sub_reg = data.anatomy.region_shifted;
elseif isfield(data.anatomy,'cluster_region')
    
    sub_reg = data.anatomy.cluster_region;
else
    sub_reg = {};
end

if iscolumn(reg)
    reg = reg';
    sub_reg = sub_reg';
end

if isfield(data.anatomy,'depth_shifted')
    depth = data.anatomy.depth_shifted;
elseif isfield(data.anatomy,'depth')
    depth = data.anatomy.depth;
end
if isfield(data.anatomy,'z2') % for MEC cases
    depth = data.anatomy.tip_distance' - data.anatomy.z2;
end

reg = reg(data.sp.cgs==2);
if ~isempty(sub_reg)
    sub_reg=sub_reg(data.sp.cgs==2);
end
depth = depth(data.sp.cgs==2);
good_cells=data.sp.cids(data.sp.cgs==2);
factors = ops.factors;
trials = ops.trials;
nT = numel(trials);
edges = ops.edges;

%create trial map that only contains numbers for trials to be included
trialMap = nan(1,numel(data.trial_gain));
cntr = 1;
for iT =1:numel(data.trial_gain)
    if ismember(iT,trials)
        trialMap(iT)=cntr;
        cntr=cntr+1;
    end
end
%recreate data.trial map that resets numbers of included trials and sets
%all else to nan
trial_sorted = nan(size(data.trial));
for iT=1:numel(trial_sorted)
    trial_sorted(iT)=trialMap(data.trial(iT));
end

[speed,speed_raw]=calcSpeed(data.posx,ops);
if ~isfield(ops,'speed_filter')
    ops.speed_filter = ops.filter;
end
if ~isempty(ops.speed_filter)
    speed_raw = conv(speed_raw,ops.speed_filter,'same');
end
% speed_raw = speed;
% occupancy matrix
OCC=zeros(nT,ops.nBins,numel(factors)); %occupancy matrix, calculates occupancy for each space bin for each shift

for iFactor = 1:numel(factors)
    posxhat = data.posx+factors(iFactor)*speed_raw;
    posxhat = mod(posxhat,max(ops.edges));
    spike_loc_hat = discretize(posxhat,edges);
    for iT=1:numel(spike_loc_hat)
        r=trial_sorted(iT);
        c=spike_loc_hat(iT);
        if ~isnan(r)
            OCC(r,c,iFactor)=OCC(r,c,iFactor)+1;
        end
    end
end

all_stability=zeros(numel(good_cells),numel(factors));
good_idx = ismember(data.sp.clu,good_cells);
clu_tmp = data.sp.clu(good_idx);
st_tmp = data.sp.st(good_idx);
[a,~,clus]=unique(clu_tmp);
nClu = numel(a);

[spMapBL]=shiftAllMapsByFactor(ops,clus,st_tmp,nClu,data.posx,data.post,trial_sorted,speed_raw,0);

firing_rate = nan(numel(good_cells,1));
idxVR = ismember(data.trial,ops.trials);
start_t = min(data.post(idxVR));
stop_t = max(data.post(idxVR));
duration = stop_t-start_t;

for cellIDX=1:numel(good_cells)
    % extract spike times for this cell
    spike_id=data.sp.clu==good_cells(cellIDX);
    spike_t = data.sp.st(spike_id);
    firing_rate(cellIDX)=sum(spike_t<stop_t & spike_t>start_t)/duration;
    % convert to VR idx
    [~,~,spike_idx] = histcounts(spike_t,data.post);
    posx=mod(data.posx,max(ops.edges));
    
    %spike_loc = discretize(posx,edges);
    idx=triu(true(nT),1);
    
    factors = ops.factors;
    VAR=zeros(size(factors));
    STAB=VAR;
    % for each shift, calculate a spatial firing rate map and calculate
    % trial by trial correlation
    % average of this correlation is the STAB for this factor
    for iFactor = 1:numel(factors)
        spMatHat = zeros(nT,ops.nBins);
        posxhat = posx+factors(iFactor)*speed_raw;
        posxhat = mod(posxhat,max(ops.edges));
        spike_loc_hat = discretize(posxhat,edges);
        for ii=1:numel(spike_idx)
            if ~isnan(trial_sorted(spike_idx(ii)))
                spMatHat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))=spMatHat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))+1;
            end
        end
        %spMatHat = medfilt1(spMatHat);
        %divide by occupancy
        spMatHat = spMatHat./squeeze(OCC(:,:,iFactor));
        spMatHat = fillmissing(spMatHat,'pchip',2);
        iidx = (size(spMatHat,2)+1):(2*size(spMatHat,2));
        
        if ~isempty(ops.filter)
            spF = [spMatHat spMatHat spMatHat];
            spF = convn(spF,ops.filter,'same');
            spMatHat = spF(:,iidx);
        end
        
        
        cc=corr(spMatHat(:,ops.idx)');
        stability=nanmean(cc(idx));
        
        STAB(iFactor)=stability;
    end
    
    
    
    all_stability(cellIDX,:)=STAB;
end
fighandles = {};
if ops.plotfig
    tmp_stability = all_stability(:,find(factors==0));
    [a,b]=sort(tmp_stability,'descend','MissingPlacement','last');
    for ii=1:15
        fighandles{ii}=figure('visible','off');
        this_cell=b(ii);
        
        
        subplot(3,1,1)
        plot(factors,all_stability(this_cell,:))
        xlim([min(factors) max(factors)])
        [ma,mi]=max(all_stability(this_cell,:));
        posxhat = posx+factors(mi)*speed_raw;
        posxhat = mod(posxhat,max(ops.edges));
        title(sprintf('%s, facto %.2f',reg{this_cell},factors(mi)))
        spike_id=data.sp.clu==good_cells(this_cell);
        spike_t = data.sp.st(spike_id);
        % convert to VR idx
        	
        
                spMatHat = zeros(nT,ops.nBins);

        for iS=1:numel(spike_idx)
            if ~isnan(trial_sorted(spike_idx(iS)))
                spMatHat(trial_sorted(spike_idx(iS)),spike_loc_hat(spike_idx(iS)))=spMatHat(trial_sorted(spike_idx(iS)),spike_loc_hat(spike_idx(iS)))+1;
            end
        end
        %spMatHat = medfilt1(spMatHat);
        %divide by occupancy
        spMatHat = spMatHat./squeeze(OCC(:,:,iFactor));
        spMatHat(isnan(spMatHat))=0;
        iidx = (size(spMatHat,2)+1):(2*size(spMatHat,2));
        
        if ~isempty(ops.filter)
            spF = [spMatHat spMatHat spMatHat];
            spF = convn(spF,ops.filter,'same');
            spMatHat = spF(:,iidx);
        end
        subplot(3,1,2)
        imagesc(spMapBL(:,:,this_cell))
        
        subplot(3,1,3)
        imagesc(spMatHat);
        

        
%         subplot(3,1,2)
%         scatter(posx(spike_idx),trial_sorted(spike_idx),2);
%         xlim([0 max(ops.edges)])
%         ylim([0 nT])
%         subplot(3,1,3)
%         scatter(posxhat(spike_idx),trial_sorted(spike_idx),2);
%         xlim([0 max(ops.edges)])
%         ylim([0 nT])

    end
end

data_out.all_stability=all_stability;
data_out.region = reg;
data_out.sub_reg = sub_reg;
data_out.depth = depth;
data_out.CID = good_cells;
data_out.firing_rate = firing_rate;
end

