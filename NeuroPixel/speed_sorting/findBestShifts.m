function data_out = findBestShifts(data,ops)


good_cells=data.sp.cids(data.sp.cgs==2);
factors = ops.factors;
trials = ops.trials;

nT = numel(trials);
edges = ops.edges;
trialMap = nan(1,numel(data.trial_gain));
cntr = 1;
for iT =1:numel(data.trial_gain)
    if ismember(iT,trials)
        trialMap(iT)=cntr;
        cntr=cntr+1;
    end
end

trial_sorted = nan(size(data.trial));
for iT=1:numel(trial_sorted)
    trial_sorted(iT)=trialMap(data.trial(iT));
end

[~,speed_raw]=calcSpeed(data.posx,ops);

if ~isempty(ops.filter)
    speed_raw = conv(speed_raw,ops.filter,'same');
end

OCC=zeros(nT,400,numel(factors));

for iFactor = 1:numel(factors)
    posxhat = data.posx+factors(iFactor)*speed_raw;
    posxhat = mod(posxhat,400);
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
if ops.plotfig
    figure
end

for cellIDX=1:numel(good_cells)
    
    spike_id=data.sp.clu==good_cells(cellIDX);
    spike_t = data.sp.st(spike_id);
    [~,~,spike_idx] = histcounts(spike_t,data.post);
    posx=mod(data.posx,400);
    spike_loc = discretize(posx,edges);
    
    idx=triu(true(nT),1);
    
    factors = ops.factors;
    VAR=zeros(size(factors));
    STAB=VAR;
    for iFactor = 1:numel(factors)
        spMatHat = zeros(nT,400);
        posxhat = posx+factors(iFactor)*speed_raw;
        posxhat = mod(posxhat,400);
        spike_loc_hat = discretize(posxhat,edges);
        for ii=1:numel(spike_idx)
            if ~isnan(trial_sorted(spike_idx(ii)))
                spMatHat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))=spMatHat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))+1;
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
        
        
        cc=corr(spMatHat');
        stability=nanmean(cc(idx));
        
        STAB(iFactor)=stability;
    end
    if ops.plotfig
        subplot(3,1,1)
        plot(factors,STAB)
        xlim([min(factors) max(factors)])
        [ma,mi]=max(STAB);
        posxhat = posx+factors(mi)*speed_raw;
        posxhat = mod(posxhat,400);
        title(sprintf('facto %.2f',factors(mi)))
        
        subplot(3,1,2)
        scatter(posx(spike_idx),trial_sorted(spike_idx),2);
        xlim([0 400])
        ylim([0 nT])
        subplot(3,1,3)
        scatter(posxhat(spike_idx),trial_sorted(spike_idx),2);
        xlim([0 400])
        ylim([0 nT])
        pause
        clf
    end
    
    
    all_stability(cellIDX,:)=STAB;
end
data_out.all_stability=all_stability;

end

function occ = getOccpuancy()
occ=[];
end