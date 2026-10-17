clear;
close all;
clc;

% Inputs: Propellant/Thermodynamic Properties
F = 5e-3;                   % Thrust (N)
Pc = 45.3e5;               % Chamber Pressure (Pa)
Cstar = 1513.13822114321;   % Cstar (m/s)
k = 1.1855;                 % Specific Heat Ratio
mu = 0.92068e-4;            % Dynamic Viscosity (Pa.s)
a = 1.7270 ;                % Vieille's Law a parameter (mm/s/bar^n)
n = 0.3821 ;                % Vieille's Law n parameter

% Mechanical parameters
lc = 3e-3 ;                 % Length of combustion chamber (m)
YieldStress = 34.138e6 ;    % CC material yield stress (Pa)
rho_mat = 1300 ;            % epoxy density (kg/m3)
Lign = 1e-4 ;               % Igniter length (m)
rho_nozzles = 10490 ;        % Nozle material density (kg/m3)

% Nozzle Input Parameters
alpha = deg2rad(30) ;       % divergent angle
eps = 100;                  % Expansion Ratio
eta_nozzle = 0.7 ;         % Nozzle correction coefficient due to low Re

% Inputs: Satellite Parameters
M0 = 4;                     % Initial Mass (kg)
dV = 0.5;                   % Delta-V (m/s)

% Constants
g0 = 9.81;                  % Gravity (m/s^2)

%% Sizing

% Thrust coefficient
Pe = invExpCoeff(k, Pc, eps) ;                      % Exit pressure
Ct = ThrustCoeff(k, eps, Pc, Pe, 0) ;               % Thrust coefficient
lambda = (1 + cos(alpha)) / 2 ;                     % Divergence loss 2d flow
Ct = Ct * lambda * eta_nozzle ;                     % Corrected Thrust coefficient

% Throat Area
At = F / (Pc * Ct);                                 % Throat Area (m^2)
dt = sqrt(4 * At / pi);                             % Throat Diameter (m)

% Mass Flow rate (kg/s)
mdot = Pc * At / Cstar ;

% Specific Impulse (s)
Isp = Ct * Cstar / g0 ;

% Reynolds
Re = 4 * mdot / pi / dt / mu ;

% Propellant Mass
M_fuel = RocketEqn(dV, Isp, M0);                            % Required Propellant Mass (kg)
rho_solid = RhoProp([0.68 0.18 0.14], [1950 2700 920]);     % Propellant Density (kg/m^3)
V_fuel = M_fuel / rho_solid ;                               % Propellant Volume (m^3)

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
twall_min = Ks * Pc * db / 2 / YieldStress ;
twall = 0.5e-3 ;

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
h = (lc + twall + Ldiv + Lign) ;

% Array Volume 
Varray = larray^2 * h ;

% Nozzle Volume
Vnzls = Nlat^2 * (1/3) * pi * Ldiv * (rt^2 + re^2 + re*rt) ;

% Array mass
Marray = (larray^2 * (lc + twall + Lign) - V_fuel) * rho_mat + (larray^2 * Ldiv - Vnzls) * rho_nozzles ;

% Output Results
fprintf('Chamber Pressure: %.3f bar\n', Pc*1e-5);
fprintf('Reynolds Number: %.3f\n', Re);
fprintf('Isp: %.3f s \n', Isp);
fprintf('Propellant Mass: %.3f g\n', M_fuel * 1000);
fprintf('Grain diameter: %.5f mm\n', db*1000);
fprintf('Throat diameter: %.5f mm\n', dt*1000);
fprintf('minimum wall thickness: %.2f mm\n', twall_min*1000);
fprintf('Selected wall thickness: %.2f mm\n', twall*1000);
fprintf('Divergent length: %.4f mm\n', Ldiv*1000);
fprintf('Array length: %.3f mm\n', larray*1000);
fprintf('Array Depth: %.3f mm\n', h*1000);
fprintf('Array total volume: %.3f mm3\n', Varray*1e9);
fprintf('Array total mass: %.3f g\n', (Marray + M_fuel)*1000);
fprintf('Number of engines: %d \n', Nlat^2);

% Igniter 0.34 w scaled from paper using grain area ratio