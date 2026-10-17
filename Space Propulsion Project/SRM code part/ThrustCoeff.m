function [Ct] = ThrustCoeff(k, epsilon, Pc, Pe, Pa)
% Computes the Thrust coefficient of an optimally expanded nozzle
% k : gas adiabatic index [pure number]
% Pc : combustion chamber pressure [any unit]
% Pe : nozzle outlet pressure [same unit as Pc]
% Pa : ambient pressure [same unit as Pc] (default Pa = Pe)
% epsilon : expansion ratio


Ct = sqrt(2*(k^2)/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1 - (Pe/Pc)^((k-1)/k))) ;
Ct = Ct + (Pe-Pa)/Pc*epsilon ;

end