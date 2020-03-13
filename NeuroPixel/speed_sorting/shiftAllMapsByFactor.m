function spMatHat = shiftAllMapsByFactor(ops,clus,st,nClu,posx,post,trial_sorted,speed_raw,factor)
nTrials = numel(ops.trials);
spMatHat = zeros(nTrials,ops.nBins,nClu);

        
        
        
        if ~isempty(ops.filter)
            if iscolumn(ops.filter)
                ops.filter = ops.filter'; % make sure filtering goes along correct dimension!
            end
        end

if numel(factor)==1
    %one factor for all maps
    % occupancy matrix
    OCC=zeros(nTrials,ops.nBins); %occupancy matrix, calculates occupancy for each space bin for each shift
    
    posxhat = posx+factor*speed_raw;
    posxhat = mod(posxhat,max(ops.edges));
    spike_loc_hat = discretize(posxhat,ops.edges);
    
    for iT=1:numel(spike_loc_hat)
        trial=trial_sorted(iT);
        pos=spike_loc_hat(iT);
        if ~isnan(trial)
            OCC(trial,pos)=OCC(trial,pos)+1;
        end
    end
    
    
    for cellIDX=1:nClu
        % extract spike times for this cell
        spike_id=clus==cellIDX;
        spike_t = st(spike_id);
        % convert to VR idx
        [~,~,spike_idx] = histcounts(spike_t,post);
        
        % for each shift, calculate a spatial firing rate map and calculate
        % trial by trial correlation
        % average of this correlation is the STAB for this factor
        tmpMat = zeros(nTrials,ops.nBins);
        
        
        for ii=1:numel(spike_idx)
            if ~isnan(trial_sorted(spike_idx(ii)))
                tmpMat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))=tmpMat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))+1;
            end
        end
        %spMatHat = medfilt1(spMatHat);
        %divide by occupancy
        
        tmpMat = tmpMat./OCC;
        tmpMat = fillmissing(tmpMat,'pchip',2);
        %tmpMat(isnan(tmpMat))=0;
        iidx = (size(tmpMat,2)+1):(2*size(tmpMat,2));
        
        
        
        if ~isempty(ops.filter)
            spF = [tmpMat tmpMat tmpMat];
            spF = convn(spF,ops.filter,'same');
            tmpMat = spF(:,iidx);
        end
        
        
        spMatHat(:,:,cellIDX)=tmpMat;
    end
    
else
    for cellIDX=1:nClu
        % occupancy matrix
        OCC=zeros(nTrials,ops.nBins); %occupancy matrix, calculates occupancy for each space bin for each shift
        
        posxhat = posx+factor(cellIDX)*speed_raw;
        posxhat = mod(posxhat,max(ops.edges));
        spike_loc_hat = discretize(posxhat,ops.edges);
        
        for iT=1:numel(spike_loc_hat)
            trial=trial_sorted(iT);
            pos=spike_loc_hat(iT);
            if ~isnan(trial)
                OCC(trial,pos)=OCC(trial,pos)+1;
            end
        end
        
        
        
        % extract spike times for this cell
        spike_id=clus==cellIDX;
        spike_t = st(spike_id);
        % convert to VR idx
        [~,~,spike_idx] = histcounts(spike_t,post);
        
        % for each shift, calculate a spatial firing rate map and calculate
        % trial by trial correlation
        % average of this correlation is the STAB for this factor
        tmpMat = zeros(nTrials,ops.nBins);
        
        
        for ii=1:numel(spike_idx)
            if ~isnan(trial_sorted(spike_idx(ii)))
                tmpMat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))=tmpMat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))+1;
            end
        end
        %spMatHat = medfilt1(spMatHat);
        %divide by occupancy
        
        tmpMat = tmpMat./OCC;
        tmpMat = fillmissing(tmpMat,'pchip',2);
        iidx = (size(tmpMat,2)+1):(2*size(tmpMat,2));
        
        if ~isempty(ops.filter)
            spF = [tmpMat tmpMat tmpMat];
            spF = convn(spF,ops.filter,'same');
            tmpMat = spF(:,iidx);
        end
        
        
        spMatHat(:,:,cellIDX)=tmpMat;
    end
    
end
end
