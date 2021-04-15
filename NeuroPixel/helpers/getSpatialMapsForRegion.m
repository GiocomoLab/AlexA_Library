function [corrMat,frMat] = getSpatialMapsForRegion(data,region,trial_range)

if isfield(data.anatomy,'parent_shifted')
    reg = data.anatomy.parent_shifted;
else
    reg = data.anatomy.cluster_parent;
end
if iscolumn(reg)
    reg=reg';
end

good_cells_idx = startsWith(reg,region) & data.sp.cgs==2;
good_cells= data.sp.cids(good_cells_idx);

ops = load_default_opt;
[corrMat,frMat,~]=trialCorrMat(good_cells,trial_range,data,ops);

end


