function [xcorrs,lags] = calc_xcorr_snippet2(spatialMap,template,startbin,stopbin,maxlag)
xcorrs = zeros(size(spatialMap,1),size(spatialMap,3),2*maxlag+1);
for iT = 1:size(spatialMap,3)
    for iC = 1:size(spatialMap,1)
        %tmp_template = template(iC,startbin:stopbin);
        if iT<5
            %firs 4
            tmp_template = squeeze(spatialMap(iC,startbin:stopbin,iT+4));
        elseif iT>4 && iT<15
            %everything with gain change
            tmp_template = squeeze(spatialMap(iC,startbin:stopbin,iT-4));
            
        elseif iT<20
            tmp_template = squeeze(spatialMap(iC,startbin:stopbin,10));
        else 
            tmp_template = squeeze(spatialMap(iC,startbin:stopbin,10));
        end
        tmp_resp = squeeze(spatialMap(iC,startbin:stopbin,iT));
        [aa,lags] = xcorr(tmp_template,tmp_resp,maxlag,'coeff');
        xcorrs(iC,iT,:)=aa;
    end
end
