function [m] = normr(m)
m=sqrt(m.^2 ./ sum(m.^2,2)) .* sign(m);
end

