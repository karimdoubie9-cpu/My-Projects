function M_fuel = RocketEqn(DV, Isp, M0)
% Uses the rocket equation to compute the mass of propellant(in kg) needed given:
% Total delta v, DV [m/s
% Specific impulse, Isp [s]
% Dry mass, M0, Kg

% assume M_fuel << M0, so 
M_ratio = exp(DV/(Isp*9.81));
M_fuel = M0*(M_ratio-1);

end