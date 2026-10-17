function [results] = VLM_simulator(param)
% Phsyical properties
R = param.R;
g0 = param.g0;
rho = param.rho;
gamma = param.gamma;

% Inputs to be varied
At = param.At;
Ae = param.Ae;
Ainj = param.Ainj;
power = param.wtot;

% Injection conditions
Ptank = param.Ptank;
Ttank = param.Ttank;
Cd = param.Cd;

% Nozzle
eps = Ae/At;
alpha = param.alpha;
%========= SIMULATION ============
load("deltaH_200kPa.mat"); % Load water deltaH tables
Pdrop_feedline = 0.001e5;
Pdrop_dyn_pressure = rho*param.plenumVel^2/2;
Pdrop_rayleigh=3.210e+03;
mdot = @(P_c) Cd*Ainj*sqrt(2*rho*(Ptank-Pdrop_feedline-Pdrop_dyn_pressure-Pdrop_rayleigh-P_c));
T0_c = @(P_c) interp1(deltaH,T1,power/mdot(P_c)+interp1(T1,deltaH,Ttank));
cstar = @(P_c) sqrt(R/gamma*T0_c(P_c))/sqrt((2/(gamma+1))^((gamma+1)/(gamma-1)));

% Solve numerically for chamber pressure
P_c = fsolve(@(P_c) P_c - mdot(P_c)*cstar(P_c)/At,...
    2e5, optimset('Display', 'off')); 

% Compute thrust
[~,T_ratio,P_ratio,~,~] = flowisentropic(gamma, eps, 'sup');

lambda=(1+cos(alpha))/2;
Cf = @(P_c) lambda*sqrt(2*gamma^2/(gamma-1)*(2/(gamma+1))^((gamma+1)/(gamma-1))...
                * (1-P_ratio^((gamma-1)/gamma)))...
                + (P_ratio*P_c)*eps*At/(P_c*At);

results.P_c = P_c;
results.T0_c = T0_c(P_c);
results.mdot = mdot(P_c);
results.cstar = cstar(P_c);
results.Cf = Cf(P_c);
results.Isp = results.cstar*results.Cf/g0;
results.thrust = results.mdot*results.Isp*g0;
results.P_exit = P_c*P_ratio;
results.T_exit = results.T0_c*T_ratio;
