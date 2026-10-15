function UAV = params()
%PARAMS  Nominal hexacopter geometric, inertial, aerodynamic and rotor data.
%
%   UAV = uav.params() returns the nominal plant parameter struct used by
%   Final.slx.  Values are those of the exam script MAIN_exam.m
%   (A3S AY 2025/2026, aerodynamic coefficients per the 04/02/2026 revision).
%
%   The struct is pure data: no gains, no scenario flags, no side effects.
%
%   See also uav.with_package, uav.allocation, harness.config.

%% ---- Rigid body -----------------------------------------------------
UAV.Nr    = 6;                          % [-]      number of rotors
UAV.m     = 1.2;                        % [kg]     vehicle mass
UAV.J     = diag([0.01 0.01 0.018]);    % [kg m^2] inertia about body origin
UAV.r_bg  = [0 0 0];                    % [m]      CoM in body frame (row)
UAV.S     = UAV.m * uav.crossmat(UAV.r_bg);   % [kg m] static moment

%% ---- Linear aerodynamics -------------------------------------------
UAV.D_fv       = diag([0.26 0.26 0]);                % linear velocity damping
UAV.D_tauomega = diag([3.6e-4 3.6e-4 1.20e-2]);      % angular velocity damping
UAV.D = [UAV.D_fv,   zeros(3);
         zeros(3),   UAV.D_tauomega];

%% ---- Generalised mass ----------------------------------------------
UAV.M    = [UAV.m*eye(3),  UAV.S';
            UAV.S,         UAV.J];
UAV.Minv = inv(UAV.M);

%% ---- Propellers -----------------------------------------------------
UAV.b         = 0.215;                  % [m]     arm length (ell in the slides)
UAV.Omega_max = 10300 * 2*pi/60;        % [rad/s] max spin rate (10300 rpm)
UAV.Omega_min =  1260 * 2*pi/60;        % [rad/s] min spin rate ( 1260 rpm)
UAV.k_m       = 1/0.05;                 % [1/s]   inverse motor time constant
UAV.k_f       = 3.65e-6;                % [N s^2/rad^2] thrust coefficient
UAV.sigma     = 0.09;                   % [m]     torque-to-thrust ratio

% Rotor angular placement, evenly spaced in the body xy-plane, and the
% spin-direction pattern (alternating, so yaw torque can be produced).
UAV.gamma = deg2rad([0 60 120 180 240 300]);   % [rad] rotor azimuth
UAV.xi    = [1 -1 1 -1 1 -1];                  % [-]   spin direction

%% ---- Control allocation --------------------------------------------
[UAV.F1, UAV.F2, UAV.F2_pinv] = uav.allocation(UAV);
end
