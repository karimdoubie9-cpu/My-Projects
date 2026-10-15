function in = input(ref, opts)
%INPUT  Build a fully specified Simulink.SimulationInput for one run.
%
%   in = harness.input(ref)
%   in = harness.input(ref, mrac=1, package=1, wind=W, lambda=L)
%   P = harness.config(startup='nominal_hover');
%   in = harness.input(ref, config=P, package=1)
%
%   Rotor initial speed comes from config; the thrust-delay prehistory is
%   matched to that speed and the initial effectiveness of the simulated plant.
%
%   Every variable the model needs is set on the SimulationInput object, so:
%     * nothing is written to the base workspace,
%     * a run cannot inherit leftovers from the previous one,
%     * the returned objects can be handed straight to parsim (Task 3.4).
%
%   The initial conditions are derived FROM the reference so the vehicle always
%   starts on the trajectory; that is what the old push_ref helper did in five
%   separate copies.
%
%   Inputs
%     ref      reference struct from traj.helix / traj.setpoint / traj.mission
%     mrac     0 = baseline, 1 = baseline + PBMRAC augmentation
%     package  0 = bare vehicle, 1 = package attached to the plant
%     wind     3-channel wind VELOCITY signal   (default: none)
%     force    3-channel wind FORCE signal      (default: none)
%     lambda   6-channel rotor effectiveness    (default: all ones)
%     disturbance  3-channel force in the INERTIAL frame, added to the
%              actuator wrench inside UAV (default: none).  Distinct from
%              'force': that one is a wind force handled by the aerodynamic
%              path, this one is an external force on the airframe.
%     config   parameter bundle from harness.config
%
%   See also harness.scenario, harness.config, harness.model_vars.

arguments
    ref            struct
    opts.mrac      (1,1) double {mustBeMember(opts.mrac,[0 1])}    = 0
    opts.package   (1,1) double {mustBeMember(opts.package,[0 1])} = 0
    opts.robust    (1,1) double {mustBeMember(opts.robust,[0 1])} = 0
    opts.wind        struct = struct([])
    opts.force       struct = struct([])
    opts.lambda      struct = struct([])
    opts.disturbance struct = struct([])
    opts.config      struct = struct([])
    opts.model     (1,:) char = 'Final'
    opts.stop_time (1,1) double = NaN
    opts.sync_body_velocity (1,1) logical = true
end

P = opts.config;
if isempty(fieldnames(P)), P = harness.config(); end

t = ref.time;
if isempty(fieldnames(opts.wind)),   opts.wind   = disturb.none(t, 3);        end
if isempty(fieldnames(opts.force)),  opts.force  = disturb.none(t, 3);        end
if isempty(fieldnames(opts.lambda)), opts.lambda = disturb.none(t, 6, 1);     end
if isempty(fieldnames(opts.disturbance)), opts.disturbance = disturb.none(t, 3); end

% Scenario-dependent entries
P.MRAC_on       = opts.mrac;
P.package_exist = opts.package;
P.robust_on     = opts.robust;
P.traj_ref      = ref;
P.wind_vel      = opts.wind;
P.wind_force    = opts.force;
P.lambda_sig    = opts.lambda;
P.Disturbance_i = opts.disturbance;

% A saved parameter bundle from before selectable startup was introduced
% retains the previous zero-speed startup when it has no motor IC field.
if ~isfield(P, 'Omega_prop_0')
    P.Omega_prop_0 = zeros(P.UAV.Nr, 1);
end
P.Omega_prop_0 = P.Omega_prop_0(:);
validateattributes(P.Omega_prop_0, {'numeric'}, ...
    {'real','finite','nonnegative','numel',P.UAV.Nr}, ...
    mfilename, 'Omega_prop_0');

% The Transport Delay is AFTER multiplication by rotor effectiveness.
% From-Workspace signals use rows for time and columns for rotor channels.
% Their first sample defines the constant prehistory assumed at startup.
lambda_0 = P.lambda_sig.signals.values(1,:).';
validateattributes(lambda_0, {'numeric'}, ...
    {'real','finite','nonnegative','numel',P.UAV.Nr}, ...
    mfilename, 'initial rotor effectiveness');
P.T_prop_0 = lambda_0 .* P.k_f .* P.Omega_prop_0.^2;

% Start on the trajectory, level, heading aligned with the initial yaw
P.p_i_0    = traj.pos(ref, 1);
P.v_i_0    = traj.vel(ref, 1);
P.q_0      = eul2quat([traj.yaw(ref, 1), 0, 0], 'ZYX')';
P.omegab_0 = [0 0 0]';

% Body-frame initial velocity.
%
% BEHAVIOUR CHANGE vs. the old scripts, which always used vb_0 = [0;0;0] while
% setting p_i_0/v_i_0 to the first trajectory sample.  On any moving reference
% (the helix starts at 0.71 m/s) that is an inconsistent initial state: the
% inertial velocity says the vehicle is moving and the body velocity says it is
% not, which injects a spurious startup transient into every run.
%
% Default true makes the two agree: vb_0 = R(q_0)' * v_i_0.
% Pass sync_body_velocity=false to reproduce the old, inconsistent behaviour.
if opts.sync_body_velocity
    P.vb_0 = quat2rotm(P.q_0.')' * P.v_i_0;
else
    P.vb_0 = [0 0 0]';
end

stop_time = opts.stop_time;
if isnan(stop_time), stop_time = t(end); end

%% ---- Validate the contract before touching Simulink ----------------
missing = setdiff(harness.model_vars(), fieldnames(P));
if ~isempty(missing)
    error('harness:input:missingVars', ...
          'Configuration is missing required model variable(s): %s', ...
          strjoin(missing, ', '));
end

%% ---- Push into the SimulationInput ---------------------------------
in = Simulink.SimulationInput(opts.model);
names = setdiff(fieldnames(P), {'design'});   % 'design' is report metadata
for i = 1:numel(names)
    in = in.setVariable(names{i}, P.(names{i}));
end
% Preserve the trajectory endpoint instead of rounding it for display.
% A rounded-down StopTime can leave the last reporting interval uncovered.
in = in.setModelParameter('StopTime', sprintf('%.17g',stop_time));
end
