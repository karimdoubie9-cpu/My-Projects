function [out] = EndBurnerSIMP(a, n, dtht, dc, Cstar)
% Simultaion of a solid rocket motor with endburner grain
% Vieille’s law "a" [mm/s/bar^n] [Random Variable]
% Vieille’s law "n" [] [Random Variable]
% dtht: Engine throat diameter [m] [Random Variable]
% dc : combustion chamber diameter [m] [Random Variable]

% Engine parameters
lg = 3e-3 ;                     % Grain lenght [m]
Ae = 8.568294238174935e-08 ;     % Exit area [m2]
k = 1.1855;                     % Specific heat capacity ratio
rho = 1761.9 ;                  % Solid propellant density [kg/m3]
alpha = deg2rad(30) ;       % divergent angle
eta_nozzle = 0.7 ;         % Nozzle correction coefficient due to low Re
 
% Convert "a" to [m/s/Pa^n]
a = 10 .^ (-3 - 5.*n) .* a ;      

% Compute throat area
Atht = 0.25*pi*(dtht .^ 2) ;

% Computes chamber area
Ac = 0.25 * pi * (dc .^ 2) ;

% Chamber pressure 
p = (Ac ./ Atht .* Cstar .* rho .* a) .^ (1 ./ (1 - n)) ;

% Burning rate [m/s]
r = a .* (p .^ n) ;

% Burning time
tb = lg ./ r ;

% Propellant mass
Mp = (rho * lg) .* Ac ;

% Nozzle exit conditions
eps = Ae ./ Atht ;
Pe = invExpCoeff(k, p, eps);

% Thrust
ct = ThrustCoeff(k, eps, p, Pe, 0) ; % Thrust coefficient
lambda = (1 + cos(alpha)) / 2 ;                     % Divergence loss 2d flow
ct = ct * lambda * eta_nozzle ;                     % Corrected Thrust coefficient
F = Atht * p * ct ;

% Total Impulse
I = F .* tb ;

% Specific Impulse
Isp = I ./ Mp ./ 9.81 ;

% Output for parallel processing
out = [tb, F, I, Isp] ;

end

