clear
clc
% --- Input Parameters & Constants ---
gamma = 1.32;        % Specific heat ratio
area_ratio = 7;      % Nozzle area ratio (Ae/At)
P0 = 2.0*10^5;         % Chamber stagnation pressure (Pa)
T_thrust = 5.3*10^-3; % Target thrust (N) - Naming collision with Temperature corrected
R = 461.5;           % Specific gas constant for water vapor (J/kg*K)
g0 = 9.81;           % Standard gravity (m/s^2)
rho = 1000;          % Propellant (water) density (kg/m^3)
Msat = 4;            % Satellite mass (kg)
dV = 0.5;            % Required delta-V (m/s)
Tns = 440*10^6;      % Tank ultimate tensile strength (Pa) - For potential stress calculation (not used below)
Te_K = 24.2 + 5 + 273; % Nozzle exit temperature estimate (K) - Source unclear, might need T0 first
Tin = 60 + 273;        % Propellant inlet temperature (K) - Assumed Celsius input
Tevap = 120 + 273;     % Propellant vaporization temperature (K) - Assumed Celsius input
hevap = 2201000;     % Latent heat of vaporization for water (J/kg)
c = 4030;            % Specific heat capacity of water (J/kg*K)
h_gas_exit = 3596*1000; % Enthalpy of superheated steam at exit state T0 (J/kg) - Approx.
h_gas_vap = 2706*1000;  % Enthalpy of saturated vapor at Tevap (J/kg) - Approx.
alpha_deg = 28;      % Nozzle divergence half-angle (degrees)
inj_coeff = 0.76;    % Injector discharge coefficient
feed_diam = 350*10^-6; % Feed line diameter (m)
tank_safety_factor = 1.1; % Safety factor for mass/volume
pressure_loss_factor = 1.15; % Factor for feed system pressure loss
P_ext = 0;       % External/back pressure estimate (Pa) - Assuming some vacuum loss
kgas = 1.66;         % Specific heat ratio for pressurant gas (e.g., Helium)
Rgas = 4124;         % Specific gas constant for Helium (J/kg*K)
T_pgb_K = 273;       % Pressurant gas initial temperature (K)
L=9*10^-3; %chamber length
w=2*10^-3; %chamber width
h=0.1*10^-3; %chamber height
Vc=L*w*h;
% --- Core Performance Calculations ---
Pc = P0; % Assuming chamber pressure equals stagnation pressure
% T0_c = @(P_c) interp1(deltaH,T1,power/mdot(P_c)+interp1(T1,deltaH,Ttank));
[mach_sup, T_ratio_sup, P_ratio_sup, ~, ~] = flowisentropic(gamma, area_ratio, 'sup'); % Using placeholder variable names from original code

P_sup = P_ratio_sup;   % Exit/Static pressure ratio (Pe/P0)
T0_r = T_ratio_sup;    % Exit/Stagnation temperature ratio (Te/T0) - Name corrected for clarity

P_exit_sup = P_sup * Pc; % Exit pressure (Pa) for supersonic, assuming vacuum exit (Pe = Pa=0) is often used in Cf calc
T0 = Te_K / T_ratio_sup; % Stagnation temperature (K) - Calculated based on estimated Te

Cf = sqrt(((2*gamma^2)/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1 - (P_sup)^((gamma-1)/gamma))) + area_ratio * (P_sup); % Assuming Pe = P_exit_sup and Pa = 0 for this standard Cf definition
% If accounting for non-zero back pressure Pa = P_ext:
% Cf_adj = sqrt(((2*gamma^2)/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1 - (P_sup)^((gamma-1)/gamma))) + area_ratio * (P_sup - P_ext/Pc);

At = T_thrust / (Pc * Cf); % Throat area (m^2)
Ae = area_ratio * At;     % Exit area (m^2)

% Exit velocity
ve = sqrt((2*gamma*R*T0)/(gamma-1) * (1 - P_sup^((gamma-1)/gamma))); % Exit velocity (m/s)


mdot = (T_thrust - (P_exit_sup - P_ext) * Ae) / ve; % Mass flow rate (kg/s)

% Specific Impulse
Isp = T_thrust / (mdot * g0); % Specific impulse (s)

% Characteristic velocity
c_star = Pc * At / mdot; % Characteristic velocity (m/s) - also c* = sqrt(R*T0/gamma) / sqrt((2/(gamma+1))^((gamma+1)/(gamma-1)))

% --- Power Calculation ---
% Heat required to raise liquid temp to boiling point
w1 = mdot * c * (Tevap - Tin); % Power for sensible heat (W)
% Heat required for vaporization
w2 = mdot * hevap;             % Power for latent heat (W)
% Heat required to superheat vapor to stagnation temperature T0
% Using enthalpies if available (h_gas_exit corresponds to T0)
w3 = mdot * (h_gas_exit - h_gas_vap); % Power for superheating (W)
% Alternative using Cp_gas if known: w3_alt = mdot * Cp_gas * (T0 - Tevap);

wtot = w1 + w2 + w3; % Total heating power (W)

% --- Propellant & Tank Sizing ---
% Propellant mass required (Rocket equation + safety margin)
Mp = (Msat * dV / (g0 * Isp)) * tank_safety_factor; % Propellant mass (kg) - Simplified rocket eq. for small dV/ve << 1, better: Mp = Msat * (exp(dV / (g0 * Isp)) - 1);
Vp = Mp / rho * tank_safety_factor; % Propellant volume (m^3)

% Tank sizing (Spherical)
rp = (3/4 * Vp / pi)^(1/3); % Inner radius of spherical tank (m)
t_p = 0.7*10^-3; % Tank wall thickness estimate (m) - Should be calculated based on stress/pressure
M_pt = (4/3) * pi * ((rp + t_p)^3 - rp^3) * 2670; % Tank mass (kg) - 
M_pt2=4*pi*rp^2*t_p*2670;
% --- Feed System & Pressure ---
A_feed = pi * (feed_diam / 2)^2; % Feed line cross-sectional area (m^2)
u_feed = mdot / (A_feed * rho); % Flow velocity in feed line (m/s)

% Total pressure required at start of feed system (tank pressure)
% Includes chamber pressure, feed line losses (factor), external pressure?, velocity head
P_tot = Pc * pressure_loss_factor + 0.001*10^5 + 0.5 * rho * u_feed^2; % Required tank pressure (Pa)

% --- Pressurant Gas System Sizing (Blowdown estimate) ---
P_pgb = 2 * P_tot; % Pressurant gas initial/beginning pressure (Pa) - Assumed factor of 2
P_pgf = P_tot;     % Pressurant gas final/end pressure (Pa) - Equals required tank pressure

% Temperatures during expansion
T_pgb = T_pgb_K; % Pressurant gas initial temperature (K)
T_pgf = T_pgb * (P_pgf / P_pgb)^((kgas - 1) / kgas); % Pressurant gas final temperature (K) - Isentropic expansion

% Volumes (assuming constant total volume V_pgf = V_gas + V_propellant)
% V_gas = V_pgf - V_propellant
V_pgf = Vp / (1 - (P_pgf * T_pgb) / (P_pgb * T_pgf)); % Total tank volume needed (m^3)
V_pgb = V_pgf - Vp; % Initial gas volume (ullage) (m^3)

% Pressurant gas mass
M_pg = P_pgb * V_pgb / (Rgas * T_pgb); % Pressurant gas mass (kg)

% --- Injector Sizing ---
P_inj = 0.15 * Pc; % Pressure drop across injector (Pa) - Assumed % of Pc
A_inj = mdot / (inj_coeff * sqrt(2 * rho * P_inj)); % Injector area (m^2)
D_inj = sqrt(A_inj / pi) * 2 * 1000; % Injector diameter (mm)

% --- Nozzle Divergence Loss ---
alpha = deg2rad(alpha_deg); % Divergence half-angle (radians)
lamda = (1 + cos(alpha)) / 2; % Divergence loss factor
% Thrust accounting for divergence loss (assuming Pa=0 in this calc)
T_2d = mdot * ve * lamda + P_exit_sup * Ae; % Thrust with divergence loss (N)

% rayleigh losses
% consider the chamber to be divided into 3 parts, liquid, mixture and
% steam
%start with the mixture part x goes from 0 to 1

% x=0.267;
cp=1872;
T0_f=w3/(cp*mdot)+Tevap;
v1=mdot/(rho*h*w);
M1=0.00001;%v1/sqrt(gamma*461*Tevap);
P_prev=Pc;
[mach,T_rtemp,P_rtemo,rho,velocity,T0_rtemp,P0_rtemp] = flowrayleigh(gamma,M1);
T0_r2=T0_f/Tevap*T0_rtemp;
 [mach,T_rtemp2,P_rtemo2,rho,velocity,T0_rtemp2,P0_rtemp2] = flowrayleigh(gamma,T0_r2,'totaltsub');
 P0_2=P0_rtemp2/P0_rtemp*Pc

% --- Display Results (Placeholder for your values) ---
fprintf('--- Performance ---\n');
fprintf('Stagnation Temperature (T0): %.2f K\n', T0);
fprintf('Exit Pressure (Pe): %.2f Pa\n', P_exit_sup);
fprintf('Thrust Coefficient (Cf): %.4f\n', Cf);
fprintf('Throat Area (At): %.4e m^2\n', At);
fprintf('Exit Area (Ae): %.4e m^2\n', Ae);
fprintf('Exit Velocity (ve): %.2f m/s\n', ve);
fprintf('Mass Flow Rate (mdot): %.4e kg/s\n', mdot);
fprintf('Specific Impulse (Isp): %.2f s\n', Isp);
fprintf('Characteristic Velocity (c*): %.2f m/s\n', c_star);
fprintf('Thrust w/ Divergence Loss (T_2d): %.4f N\n', T_2d);

fprintf('\n--- Power ---\n');
fprintf('Sensible Heat Power (w1): %.4f W\n', w1);
fprintf('Latent Heat Power (w2): %.4f W\n', w2);
fprintf('Superheat Power (w3): %.4f W\n', w3);
fprintf('Total Power (wtot): %.4f W\n', wtot);

fprintf('\n--- Propellant & Tank ---\n');
fprintf('Propellant Mass (Mp): %.4f kg\n', Mp);
fprintf('Propellant Volume (Vp): %.4e m^3\n', Vp);
fprintf('Tank Inner Radius (rp): %.4f m\n', rp);
fprintf('Tank Mass (M_pt): %.4f kg \n', M_pt);

fprintf('\n--- Pressurant System ---\n');
fprintf('Required Tank Pressure (P_tot): %.2f Pa\n', P_tot);
fprintf('Pressurant Gas Start P (P_pgb): %.2f Pa\n', P_pgb);
fprintf('Pressurant Gas End T (T_pgf): %.2f K\n', T_pgf);
fprintf('Total Tank Volume (V_pgf): %.4e m^3\n', V_pgf);
fprintf('Initial Gas Volume (V_pgb): %.4e m^3\n', V_pgb);
fprintf('Pressurant Gas Mass (M_pg): %.4e kg\n', M_pg);

fprintf('\n--- Feed & Injector ---\n');
fprintf('Feed Line Velocity (u_feed): %.4f m/s\n', u_feed);
fprintf('Injector Diameter (D_inj): %.4f mm\n', D_inj);