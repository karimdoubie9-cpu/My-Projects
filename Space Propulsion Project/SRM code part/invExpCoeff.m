function [Pe] = invExpCoeff(k, Pc, epsilon)
% Computes the expansion coefficient of an optimally expanded nozzle
% k : gas adiabatic index [pure number]
% Pc : combustion chamber pressure [any unit]
% Pe : desired outlet pressure [same unit as Pc]
% epsilon : expansion coefficient

C = ((k+1)/2)^(1/(k-1)) * sqrt((k+1)/(k-1)) ;

fun = @(Rp) Rp^(1/k) * C * sqrt(1 - Rp^((k-1)/k)) - 1/epsilon ;
dfun = @(Rp) (1/k) * Rp^((1-k)/k) * C * sqrt(1 - Rp^((k-1)/k))...
    - C * (k-1) / 2 / k / sqrt(1 - Rp^((k-1)/k)) ;

Rp = Newton(0.2,100,1e-12,fun,dfun,0.5) ;
%Rp = real(Rp) ;

Pe = Rp * Pc ;

end
