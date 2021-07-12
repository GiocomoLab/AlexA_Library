function receptive_fields = fit_receptive_field(spike_triggered_average,z_score_t)

Z = zscore(spike_triggered_average,[],'All');
Z(abs(Z)<z_score_t)=0;
receptive_fields = {};
for iMult = [-1 1]

Z_q = Z*iMult>0;

binaryImage = imfill(Z_q, 'holes');

params = get_gaussian_fit(binaryImage,Z);

if numel(fieldnames(params))>0
    params.field_sign = iMult;
    receptive_fields{end+1}=params;
end
end