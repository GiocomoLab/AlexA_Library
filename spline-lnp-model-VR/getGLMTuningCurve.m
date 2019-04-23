function [tuning_curves,x_vals] = getGLMTuningCurve(A,variables,parameters,ctl_pts_all,s,dt)

%% Description
% Given the variables, A, and the parameters,
% this will return the tuning curves for the cell
% if plotfig = 1, this will also plot the tuning curves

% NOTE: I just use A to compute the correct indexes
numVar = max(variables);

variables = sort(variables);
b0 = parameters(1);
param = parameters(2:end);
scale = zeros(1,numVar);
total_ind = 0;
%extract the marginal values for each model (in case the best model is a
%combination
x_vals={};
y_vals={};
for iv=1:max(variables)
    if ismember(iv,variables)
    param_ind = size(A{iv},2);
    start_ind = total_ind +1;
    stop_ind = total_ind + param_ind;
    total_ind = total_ind+param_ind;
    param_marginal = param(start_ind:stop_ind);
    scale(iv) = mean(exp(A{iv}*param_marginal'));
    [y,x] = spline_1d_plot(param_marginal,ctl_pts_all{iv},s);
    x_vals{iv}=x;
    y_vals{iv}=y;
    else
        scale(iv)=NaN;
    end
end
    
    
    


tuning_curves = {};
for var_k = 1:max(variables)
    
    if ismember(var_k,variables)
        y = y_vals{var_k};
        scale_factor_ind = setdiff(variables,var_k); scale_factor = scale(scale_factor_ind);
        tuning_curves{var_k} = exp(y(:))*exp(b0)*prod(scale_factor)/dt;
    end
    
end
return