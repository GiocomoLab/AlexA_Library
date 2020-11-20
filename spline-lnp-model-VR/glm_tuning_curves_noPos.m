function [tuning_curves] = glm_tuning_curves_noPos(A,variables,parameters,ctl_pts_all,s,dt)

%% Description
% Given the variables, A, and the parameters,
% this will return the tuning curves for the cell
% if plotfig = 1, this will also plot the tuning curves

% NOTE: I just use A to compute the correct indexes

numVar = numel(A);
num_plot_columns = numVar+1;

variables = sort(variables);
b0 = parameters(1);
param = parameters(2:end);

total_ind = 0;
y_points = {};
% position




total_ind = 0;

for iv=1:max(variables)
    if ismember(iv,variables)
    param_ind = size(A{iv},2);
    total_ind = total_ind + param_ind;
    param1 = param(total_ind - param_ind + 1:total_ind);
    scale(iv) = mean(exp(A{iv}*param1'));
    if iv<6
    [y,x] = spline_1d_plot(param1,ctl_pts_all{iv},s);
    y_points{iv}=y;
    else
        y_points{iv}=nan;
    end
    else
        scale(iv)=NaN;
    end
end

tuning_curves = {};

for iV=1:max(variables)
    
    if ismember(iV,variables)
        scale_factor_ind = setdiff(variables,iV); 
        scale_factor = scale(scale_factor_ind);
        
        tuning_curves{iV} = exp(y_points{iV})*exp(b0)*prod(scale_factor)/dt;
    end
    
end
