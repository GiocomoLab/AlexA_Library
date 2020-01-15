
spatialMap(isnan(spatialMap))=0;
[regs,~,idx] = unique(region(sp.cgs==2));
flatMaps = {};
nTrials = size(spatialMap,3);

posWindow=[60 140];
posBin = [find(edges == posWindow(1)), find(edges == posWindow(2))];
nBins=posBin(2)-posBin(1)+1;

filt = gausswin(5);
filt = filt/sum(filt);
for iR=1:numel(regs)
    nC=nnz(idx==iR);
    cellList = find(idx==iR);
    if nC<10
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

%ff=xcorr(flatMaps{6}');
speed = calcSpeed(posx,params);
trial_speed = zeros(1,max(trial));
for ii = 1:max(trial)
    idx = posx>posWindow(1) & posx<posWindow(2) & trial == ii;
    tmp = mean(speed(idx));
    
    trial_speed(ii)=tmp;
end


 delay_list=zeros(nTrials*(nTrials-1)/2,6); %delay, speed_diff, %i gain, i contrast, jgain j contrast
    idx = 0;
    max_lag=20;
    for iT = 1:nTrials
        for jT=(iT+1):nTrials
            idx = idx+1;
            rr=xcorr(flatMaps{3}(iT,:),flatMaps{3}(jT,:),max_lag);
            [mm,midx]=max(rr);
            if ismember(midx,[1 size(rr,2)])
                midx = nan;
            end
            
            
            tmp_d = max_lag+1;
            delay_list(idx,:)=[midx-tmp_d,trial_speed(iT)-trial_speed(jT), trial_gain(iT),trial_contrast(iT),trial_gain(jT),trial_contrast(jT) ];
        end
    end
