clear;
close all;
clc;

% Inputs: Propellant/Thermodynamic Properties
alpha = deg2rad(30) ;   % divergent angle
F = 5e-3;               % Thrust (N)
Pc = 100e5;             % Chamber Pressure (Pa)
Tc = 3431;              % Chamber Temperature (K)
MM = 25.832;            % Molar Mass (g/mol)
k = 1.1835;             % Specific Heat Ratio
eps = 100;              % Expansion Ratio
mu = 0.94893e-4;        % Dynamic Viscosity (Pa.s)
a = 1.7270 ;            % Vieille's Law a parameter (mm/s/bar^n)
n = 0.3821 ;            % Vieille's Law n parameter
lc = 3e-3 ;             % Length of a combustion chamber
YieldStress = 34.138e6 ;
rho_mat = 1300 ;        % epoxy density
eta_nozzle = 0.75 ;

% Inputs: Satellite Parameters
M0 = 4;                 % Initial Mass (kg)
dV = 0.5;               % Delta-V (m/s)

% Constants
g0 = 9.81;              % Gravity (m/s^2)
r = 8314.462618 ;               % Universal Gas Constant (J/kmol-K)

% Throat Calculations
lambda = (1 + cos(alpha)) / 2 ;                     % Divergence loss 2d flow
Rstar = r / MM;                                     % Specific Gas Constant (J/kg-K)
Pt = Pcrit(Pc, k);                                  % Throat Pressure (Pa)
[Tt, rho_c, rho_t] = adiabatic_expansion(Pc, Pt, Tc, k, Rstar); % Throat Conditions
wt = gasAdiabaticVelocity(Pc, Pt, Tc, k, Rstar);    % Throat Velocity (m/s)
ct = eta_nozzle * lambda * ThrustCoeff(k, eps, Pc, invExpCoeff(k, Pc, eps), 0) ; % Thrust coefficient
At = F / (Pc * ct);                                 % Throat Area (m^2)
dt = sqrt(4 * At / pi);                             % Throat Diameter (m)

% Exit Calculations
Pe = invExpCoeff(k, Pc, eps);                       % Exit Pressure (Pa)
we = gasAdiabaticVelocity(Pc, Pe, Tc, k, Rstar);    % Exit Velocity (m/s)
mdot = we * rho_c * (Pe/Pc)^(1/k) * eps * At;       % Mass Flow Rate (kg/s)
Isp = F / (mdot * g0);                              % Specific Impulse (s)

% Propellant Mass Calculation
M_fuel = RocketEqn(dV, Isp, M0);                    % Required Propellant Mass (kg)
rho_solid = RhoProp([0.68 0.18 0.14], [1950 2700 920]);     % Propellant Density (kg/m^3)
V_fuel = M_fuel / rho_solid ;                       % Propellant Volume (m^3)

% Reynolds Number Calculation
Re = rho_t * wt * dt / mu;                          % Reynolds Number

% Burning rate
r = a * (Pc*1e-5)^n ;                               % mm/s
r = r / 1000 ;                                      % m/s

% required burning area
Ab = mdot / r / rho_solid ;
db = sqrt(4 * Ab / pi) ;

% Grain length
l = V_fuel / Ab ;

% Burning time 
tb = l / r ;

% Number of thrusters
Ntot = l / lc ;
Nlat = ceil(sqrt(Ntot)) ;

% Minimum Thickness
Ks = 2 ;
twall = Ks * Pc * db / 2 / YieldStress ;

% Nozzle size
Ae = eps * At ;
Ac = Ab ;
% Radii
rc = sqrt(Ac / pi) ;
rt = sqrt(At / pi) ;
re = sqrt(Ae / pi) ;

% Nozzle lenght
Ldiv = (re - rt) / tan(alpha) ;

% array lenght
larray = Nlat * db +  (Nlat + 1) * twall ;

% Array dept
h = (lc + twall + 1.05*Ldiv) ;

% Array Volume 
Varray = larray^2 * h ;

% Nozzle Volume
Vnzls = Nlat^2 * (1/3) * pi * Ldiv * (rt^2 + re^2 + re*rt) ;

% Empty space 
Vempty = Nlat^2 * 0.05*Ldiv * pi * rc^2 ;

% Array mass
Marray = (Varray - V_fuel - Vnzls - Vempty) * rho_mat ;

% Output Results
fprintf('Reynolds Number: %.2f\n', Re);
fprintf('Isp: %.2f s \n', Isp);
fprintf('Propellant Mass: %.2f g\n', M_fuel * 1000);
fprintf('Grain diameter: %.5f mm\n', db*1000);
fprintf('Throat diameter: %.5f mm\n', dt*1000);
fprintf('wall thickness: %.2f mm\n', twall*1000);
fprintf('Divergent length: %.3f mm\n', Ldiv*1000);
fprintf('Array length: %.2f mm\n', larray*1000);
fprintf('Array Depth: %.2f mm\n', h*1000);
fprintf('Array total volume: %.2f mm3\n', Varray*1e9);
fprintf('Array total mass: %.2f g\n', (Marray + M_fuel)*1000);
fprintf('Number of engines: %d \n', Nlat^2);

% Igniter 0.34 w scaled from paper using grain area ratio