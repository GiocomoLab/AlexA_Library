function [m] = normc(m)
m=sqrt(m.^2 ./ sum(m.^2)) .* sign(m);
end

