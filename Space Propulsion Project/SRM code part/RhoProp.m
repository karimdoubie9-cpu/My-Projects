function [rho_p] = RhoProp(mf, rho)
% Computes the average solid density of a solid propellant
% mf is a vector of mass factions (percentage in mass of each component)
% rho is a vector containing the density of each component
% rho and mf must have the same number of elements and be alligned
% the measurement unit used in rho must be consistent within rho and will
% be the output unit for the propellant density

N = length(mf) ;
if N ~= length(rho)
    error("Inconsistent mass fraction vector and the density vector dimensions")
end
if (sum(mf) > 1+1e-3) || (sum(mf) < 1-1e-3)
    error("the sum of all the mass fractions must be 1 (1e-3 toll)")
end

rho_p = sum(mf) / sum(mf./rho) ;

end

