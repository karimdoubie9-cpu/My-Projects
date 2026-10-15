function P = config(opts)
%CONFIG  Build the complete, self-consistent parameter set for Final.slx.
%
%   P = harness.config()
%   P = harness.config(package_mass=0.3, attitude_gains='designed')
%   P = harness.config(package_mass=0.4, startup='nominal_hover')
%
%   startup='stopped' (default) starts the motors from zero speed.
%   startup='nominal_hover' initializes them at the bare UAV hover speed.
%   Package mass is NEVER used to compute this initial speed. Position,
%   velocity and heading are still taken from each run's reference in input.
%
%   Returns ONE struct whose field names are exactly the workspace variable
%   names the model resolves.  Nothing is written to the base workspace and
%   nothing is read from it -- harness.input pushes P into a
%   Simulink.SimulationInput, so each run carries its own complete
%   configuration and runs cannot contaminate each other.
%
%   Every value is assigned exactly once.  Where the old main.m assigned a
%   name twice, the SECOND (effective) value is the one reproduced here.
%
%   See also harness.input, harness.scenario, harness.model_vars.

arguments
    opts.package_mass   (1,1) double = 0.4        % [kg] as in the old main.m
    opts.package_r_c    (3,1) double = [0.02; 0.02; -0.15]
    opts.package_L      (1,1) double = 0.2
    opts.attitude_gains (1,:) char {mustBeMember(opts.attitude_gains, ...
                                    {'as_flown','designed'})} = 'designed'
    opts.adapt_mass     (1,1) double = NaN        % NaN -> nominal UAV mass
    opts.lambda_sensor  (3,1) double = [0.017; 0.017; 0.008]
    opts.startup        (1,:) char {mustBeMember(opts.startup, ...
                                    {'stopped','nominal_hover'})} = 'stopped'
end

%% ---- Plant ---------------------------------------------------------
UAV = uav.params();
UAV_pkg = uav.with_package(UAV, opts.package_mass, ...
                           r_c = opts.package_r_c, L = opts.package_L);

%% ---- Environment / shared constants --------------------------------
par.g   = 9.81;                 % [m/s^2]
par.e_3 = [0 0 1]';

%% ---- Attitude loop --------------------------------------------------
att = design.attitude(UAV, opts.attitude_gains);
par.k_R1 = att.k_R(1);
par.k_R2 = att.k_R(2);
par.k_R3 = att.k_R(3);
ctrl.Komega = att.Komega;

%% ---- Position loop --------------------------------------------------
pos = design.position(UAV);

% Horizontal channel (x and y share a gain set: the two axes are dynamically
% identical, same mass and same D_fv entry, and both are actuated by tilting).
ctrl.K_p = pos.K_p;
ctrl.K_v = pos.K_v;
ctrl.K_i = pos.K_i;

% Vertical channel, tuned separately. z is NOT dynamically equivalent to x/y:
% it is driven by the collective directly instead of through the attitude
% cascade, it carries the whole weight rather than a tilt component, and its
% drag coefficient D_fv(3,3) is zero. A separate gain set is therefore the
% physically honest choice rather than a convenience.
%
% Copied only when design.position supplies them, so an older or reduced
% version of that function does not break the configuration.
if isfield(pos, 'K_p_Z'), ctrl.K_p_Z = pos.K_p_Z; end
if isfield(pos, 'K_v_Z'), ctrl.K_v_Z = pos.K_v_Z; end
if isfield(pos, 'K_i_Z'), ctrl.K_i_Z = pos.K_i_Z; end

%% ---- Adaptive augmentation -----------------------------------------
adapt = design.adaptive(UAV, par.g, ...
                        adapt_mass   = opts.adapt_mass);

%% ---- Assemble the workspace contract -------------------------------
P.UAV     = UAV;
P.UAV_pkg = UAV_pkg;
P.ctrl    = ctrl;
P.par     = par;
P.adapt   = adapt;
P.k_f     = UAV.k_f;

% Default flags and initial conditions; harness.input overrides the ones that
% depend on the scenario (MRAC_on, package_exist, p_i_0, v_i_0, q_0, ...).
P.MRAC_on       = 0;
P.robust_on     = 0;
P.package_exist = 0;
P.p_i_0    = [0 0 0]';
P.v_i_0    = [0 0 0]';
P.q_0      = [1 0 0 0]';
P.omegab_0 = [0 0 0]';
P.vb_0     = [0 0 0]';

% Actual motor states and consistent prehistory for the thrust delay.
% These are plant initial conditions, not additional controller parameters.
% In a loaded run, nominal_hover still leaves the unknown payload imbalance.
switch opts.startup
    case 'stopped'
        P.Omega_prop_0 = zeros(UAV.Nr, 1);
    case 'nominal_hover'
        Omega_hover = sqrt(UAV.m * par.g / (UAV.Nr * UAV.k_f));
        P.Omega_prop_0 = Omega_hover * ones(UAV.Nr, 1);
end
P.T_prop_0 = P.k_f .* P.Omega_prop_0.^2;
% harness.input updates T_prop_0 using the initial rotor-effectiveness signal.

% Task 4 switches of the single model (harness.build_single_model). These
% defaults leave Tasks 1-3 exactly as before: no floor, and PBMRAC enabled
% directly by MRAC_on (the discrete flight gate is compiled out).
P.task4_ground_z = -Inf;               % [m] floor height; -Inf = no floor
P.pb_flight_gate = 0;                  % 0 direct MRAC_on, 1 measured flight gate
P.pb_gate_Ts     = 0.01;               % [s] gate sample time
P.pb_gate_cfg    = [0; 0.15; 0.02; 0.10]; % [ground_z; enable h; enable vz; confirm t]

% Design metadata: not read by the model, carried along for the reports.
P.design.attitude = att;
P.design.position = pos;
P.design.startup  = opts.startup;
end
