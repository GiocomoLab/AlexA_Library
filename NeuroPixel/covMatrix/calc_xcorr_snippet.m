function [xcorrs,lags] = calc_xcorr_snippet(spatialMap,template,startbin,stopbin,maxlag,trial2templateMap)
xcorrs = zeros(size(spatialMap,1),size(spatialMap,3),2*maxlag+1);
n_sets = numel(template);
for iT = 1:size(spatialMap,3)
    if isempty(trial2templateMap)
    if iT<n_sets
        template_set = iT;
    else
        template_set = n_sets;
    end
    else
        template_set = trial2templateMap(iT);
    end
    for iC = 1:size(spatialMap,1)
        tmp_template = template{template_set}(iC,startbin:stopbin);
        tmp_template = tmp_template-mean(tmp_template);
        tmp_resp = squeeze(spatialMap(iC,startbin:stopbin,iT));
        tmp_resp = tmp_resp-mean(tmp_resp);
        [aa,lags] = xcorr(tmp_template,tmp_resp,maxlag,'coeff');
        xcorrs(iC,iT,:)=aa;
    end
end
