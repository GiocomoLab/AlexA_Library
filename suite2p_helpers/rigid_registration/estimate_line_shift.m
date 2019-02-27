function [shift]=estimate_line_shift(template)
% this function attempts to estimate the shift between alternating lines 
% using cross-correlation function.


for ind=1:size(template,1)-1
    if rem(ind,2)==1
        line_xcorr(ind,:) = xcorr(template(ind,:)-mean(template(ind,:)),...
            template(ind+1,:)-mean(template(ind+1,:)),'unbiased');
    else
        line_xcorr(ind,:) = flipdim(xcorr(template(ind,:)-mean(template(ind,:)),...
            template(ind+1,:)-mean(template(ind+1,:)),'unbiased'),2);
    end
end
boundary=floor(size(template,2)/2);
[~,b]=max(mean(line_xcorr(:,boundary+1:end-boundary)));
shift=round(size(line_xcorr,2)/2)-b-boundary;
