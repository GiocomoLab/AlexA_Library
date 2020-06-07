function [corrMat,shiftMat]=spMapXcorr(spatialMap,maxLag,binWidth)


nC=size(spatialMap,3);
nT = size(spatialMap,1);
corrMat = nan(nC,nT,nT);
shiftMat = nan(nC,nT,nT);
shift_all = -maxLag:binWidth:maxLag;
ml = maxLag/binWidth;
for cellIDX = 1:size(spatialMap,3)
    
    mS=squeeze(spatialMap(:,:,cellIDX))';
    
    mS=mS-mean(mS);
    %mS=zscore(mS);
    %%
    [xcorr_this,~]=xcorr(mS,'coeff',ml);
    [xcorr_this,max_idx] = max(xcorr_this,[],1); % take max over lags
    invalid_idx = ismember(max_idx,[1 2*ml+1]);
    shift_this = shift_all(max_idx);
    shift_this(invalid_idx)=nan;
    xcorr_this(invalid_idx)=nan;
    
    xcorr_this = reshape(xcorr_this,nT,nT);
    shiftMat(cellIDX,:,:) = reshape(shift_this,nT,nT);
    %corrMat(cellIDX,:,:) = xcorr_this-diag(diag(xcorr_this)); % subtract diags
    corrMat(cellIDX,:,:) = xcorr_this+diag(nan(nT,1));
end
end
