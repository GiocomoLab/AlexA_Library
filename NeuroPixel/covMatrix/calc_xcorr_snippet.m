function [xcorrs,lags] = calc_xcorr_snippet(spatialMap,template,startbin,stopbin,maxlag)
xcorrs = zeros(size(spatialMap,1),size(spatialMap,3),2*maxlag+1);
for iT = 1:size(spatialMap,3)
    for iC = 1:size(spatialMap,1)
        tmp_template = template(iC,startbin:stopbin);
        tmp_resp = squeeze(spatialMap(iC,startbin:stopbin,iT));
        tmp_resp = tmp_resp-mean(tmp_resp);
        [aa,lags] = xcorr(tmp_template,tmp_resp,maxlag,'coeff');
        xcorrs(iC,iT,:)=aa;
    end
end
