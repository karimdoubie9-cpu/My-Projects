function [epsilon] = ExpCoeff(k, Pc, Pe)
% Computes the expansion coefficient of an optimally expanded nozzle
% k : gas adiabatic index [pure number]
% Pc : combustion chamber pressure [any unit]
% Pe : desired outlet pressure [same unit as Pc]

eps_inv = ((k+1)/2)^(1/(k-1)) * (Pe/Pc)^(1/k) * sqrt((k+1)/(k-1)*(1 - (Pe/Pc)^((k-1)/k))) ;
epsilon = 1/(eps_inv) ;

end

