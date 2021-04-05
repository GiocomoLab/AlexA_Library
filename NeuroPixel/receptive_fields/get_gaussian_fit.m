%function [rec_field_outline,rec_field_center,rec_field_sign] = fit_receptive_field(spike_triggered_average,z_score_t)


function params = get_gaussian_fit(binaryImage,Z)
[labeledImage, ~] = bwlabel(binaryImage);
blobMeasurements = regionprops(labeledImage, 'area', 'Centroid','PixelList','PixelIdxList');

[ma,mi]=max([blobMeasurements(:).Area]);
params = struct();
if ma>50
    ix = blobMeasurements(mi).PixelList;
    lam = Z(blobMeasurements(mi).PixelIdxList);
    params = fitMVGaus(ix(:,1),ix(:,2),lam,2.5);
    params.PixelList = blobMeasurements(mi).PixelList;
end

end
    
    