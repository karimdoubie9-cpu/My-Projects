function [p_star] = Pcrit(p0, k)
% Pcrit computes the choking pressure, for which M = 1
% k is the heat capacity ratio (AKA gamma)
% p0 is the nozzle inlet, tank or stagnation pressure
% p_start is the pressure at which choked flow occurs

p_star = p0 * (2 / (k + 1)) ^ (k / (k - 1)) ;

end

