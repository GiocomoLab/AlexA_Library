function [corrMat,shiftMat] = calc_xcorr_snippet(spatialMap,template,startbin,stopbin,maxlag,trial2templateMap)
xcorrs = zeros(size(spatialMap,1),size(spatialMap,2),2*maxlag+1);
corrMat = zeros(size(spatialMap,1),size(spatialMap,2));
shiftMat = corrMat;
n_sets = numel(template);
for iT = 1:size(spatialMap,2)

        template_set = trial2templateMap(iT);
    for iC = 1:size(spatialMap,1)
        tmp_template = template{template_set}(iC,startbin:stopbin);
        tmp_template = tmp_template-mean(tmp_template);
        tmp_resp = squeeze(spatialMap(iC,iT,startbin:stopbin));
        tmp_resp = tmp_resp-mean(tmp_resp);
        [aa,lags] = xcorr(tmp_template,tmp_resp,maxlag,'coeff');
        
        [xcorr_this,max_idx] = max(aa); % take max over lags
        shiftMat(iC,iT)=max_idx;
        corrMat(iC,iT)=xcorr_this;
        %xcorr_this = reshape(xcorr_this,numel(trials),numel(trials));
        %shiftMat(i,:,:) = reshape(shift_this,numel(trials),numel(trials)); 
        %corrMat(i,:,:) = xcorr_this+diag(nan(numel(trials),1)); % make diagonals nan
        
        %xcorrs(iC,iT,:)=aa;
    end
end
