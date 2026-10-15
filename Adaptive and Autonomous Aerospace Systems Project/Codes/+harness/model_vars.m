function [required, dead] = model_vars()
%MODEL_VARS  The workspace contract of Final.slx, stated explicitly.
%
%   [required, dead] = harness.model_vars()
%
%   REQUIRED lists every base-workspace variable that a LIVE (non-commented)
%   block of Final.slx resolves.  harness.input validates against this list, so
%   a scenario that forgets to define one fails loudly instead of silently
%   inheriting whatever the previous run left behind.
%
%   DEAD lists names that the old scripts assigned but that no live block reads.
%   They are kept here as documentation -- deliberately NOT set by
%   harness.config -- so nobody re-adds them thinking they matter.
%
%   Derived by walking the unpacked .slx block tree and skipping every subtree
%   whose Commented flag is 'on' (notably the whole 'ADAPTIVE PREDICTOR'
%   subsystem, which is a disabled earlier version of the augmentation).

required = { ...
    ... % --- parameter structs -------------------------------------------
    'UAV'              % .m .S .D .M .Minv .F1 .F2 .F2_pinv .Omega_max .k_m
    'UAV_pkg'          % .m .S .M .Minv
    'ctrl'             % .K_p .K_v .K_i (horizontal) .Komega
                       %   plus .K_p_Z .K_v_Z .K_i_Z for the separate
                       %   vertical channel, when design.position defines
                       %   them. Not in the required list: the model only
                       %   reads them once its z Gain blocks point at them.
    'par'              % .g .e_3 .k_R1 .k_R2 .k_R3
    'adapt'            % .m .g .L_p .L_v .Gamma_theta .Gamma_lambda .sigma_e
                       %   .lambda_proj_center/_r0/_r1 (sigma_lambda: unused)
    ... % --- scalars ------------------------------------------------------
    'k_f'              % thrust coefficient (actuator + allocation blocks)
    'MRAC_on'          % request to enable the adaptive controller (through the
                       %   'PBMRAC enable' variant, see pb_flight_gate)
    'package_exist'    % selects the 'package' / 'no package' plant variant
    'robust_on'        % e-modification on theta_hat + projection on Lambda_hat
    ... % --- initial conditions -------------------------------------------
    'p_i_0'            % also seeds the predictor integrator
    'v_i_0'
    'q_0'
    'omegab_0'
    'vb_0'
    'Omega_prop_0'     % actuator speed Integrator initial condition [rad/s]
    'T_prop_0'         % thrust Transport Delay initial output, per rotor [N]
    ... % --- Task 4 switches (defaults in harness.config) --------------------
    'task4_ground_z'   % floor height [m]; -Inf = no floor (Tasks 1-3)
    'pb_flight_gate'   % 'PBMRAC enable' variant: 0 direct MRAC_on, 1 flight gate
    'pb_gate_Ts'       % flight-gate sample time [s] (used only when gated)
    'pb_gate_cfg'      % [ground_z; enable h; enable vz; confirm t]
    ... % --- From-Workspace signals ---------------------------------------
    'traj_ref'         % Switch1(642) selects this over the internal generator
    'wind_vel'
    'wind_force'
    'lambda_sig'
    'Disturbance_i'    % 3-channel force in the INERTIAL frame; UAV/MatrixMultiply
                       %   rotates it into the body (R') and UAV/Sum adds it to
                       %   the actuator wrench, padded with zero torque
    };

dead = { ...
    ... % Gains and constants inside the commented-out ADAPTIVE PREDICTOR
    'A_p'
    'B_p'
    'adapt.L'
    'theta_bar'
    'FLAG_adapt'
    'proj_FLAG'
    ... % Projection-operator parameters: same disabled subsystem, and several
    ... % of these were never defined anywhere, which is why the model only
    ... % compiles with that subsystem commented out.
    'par.center_lambda'
    'par.center_theta_p'
    'par.theta_M_lambda'
    'par.eps_proj'
    'par.theta_M_theta_p'
    'par.Gamma_lambda'
    'par.Gamma_theta_p'
    ... % The live integrators hard-code IC 0 and 1, so these do nothing
    'theta_hat_0'
    'Lambda_hat_diag_0'
    ... % Never referenced by any block, live or dead
    'x_hat_0'
    'p_0'
    'stop_time'          % StopTime is a model parameter, not a variable
    'UAV.lambda_act'     % effectiveness comes from the lambda_sig signal
    'par.K_R'            % the zeta/omega_n design was wired to nothing;
    'par.K_omega'        %   see design.attitude for what the model uses
    'par.eul_d'
    'adapt.lambda_min'
    'adapt.lambda_max'
    };
end
