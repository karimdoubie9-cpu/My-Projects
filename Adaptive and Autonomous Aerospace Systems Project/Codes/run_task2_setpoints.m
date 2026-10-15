%% ========================================================================
%  TASK 2.3 - BASELINE CONTROLLER: TIME-DOMAIN TESTING
%  ------------------------------------------------------------------------
%  THE CONTROL ARCHITECTURE UNDER TEST  (P2CS1, cascaded / inner-outer)
%
%     trajectory --> POSITION ctrl --> ATTITUDE PLANNER --> ATTITUDE ctrl
%                                                              |
%                                              CONTROL ALLOCATION --> rotors
%

%  OUTER LOOP (position). PID on the position error, with the integral term
%  present specifically to reject a CONSTANT force bias - the unmodelled
%  weight of the package in Task 3, or a steady headwind:
%
%     f_c = m*( a_d - K_p*e_p - K_v*e_v - K_i*INT(e_p) ) + m*g*e_3
%     e_p = p - p_d,   e_v = v - v_d
%
%  ATTITUDE PLANNER. A coplanar rotor set makes force only along body z, so
%  the desired thrust direction IS the desired body z axis: b3_d = f_c/|f_c|.
%  That plus the commanded heading fixes R_d completely.
%
%  INNER LOOP (attitude). Geometric control on SO(3)
%
%     e_R  = sum_i k_Ri * ( R_e'*e_i  x  e_i ),   R_e = R_d'*R
%     tau_c = -e_R - K_omega*(omega - omega_d)
%
%  ALLOCATION. T_r = pinv(F2) * [T; tau], minimum-norm across the 6 rotors,
%  then Omega_r = sqrt(T_r/k_f).
%
%  ------------------------------------------------------------------------
%  WHAT EACH SCENARIO IS FOR (all BASELINE, MRAC_on = 0, no package)
%     A. lateral step  - setpoint case used to tune, horizontal axis
%     B. vertical step - same, but the vertical axis has no attitude
%                        cascade in front of it, so it is the faster loop
%     C. helix         - trajectory tracking: persistent, non-zero a_d, so
%                        it exposes tracking lag that a step never shows
%     D. C + reduced propeller effectiveness            [Task 2.3.2]
%     E. C + wind gust
%
%  Each scenario declares its OWN reference and its OWN disturbances, so no
%  run can inherit anything left over from the one before it.
%  ========================================================================

clear; close all; clc;
viz.defaults();   % Times New Roman, heavier lines (see +viz/defaults.m)

P = harness.config();          % one parameter set, shared by every run below

fprintf('Baseline gains: K_p = %.3f, K_v = %.3f, K_i = %.3f\n', ...
        P.ctrl.K_p, P.ctrl.K_v, P.ctrl.K_i);
fprintf('Attitude gains (%s): k_R = [%.3g %.3g %.3g], effective K_R = [%.3g %.3g %.3g]\n\n', ...
        P.design.attitude.mode, P.design.attitude.k_R, P.design.attitude.Keff_R);

%% ---- A. Lateral step -----------------------------------------------
%  A 1 m step in x from hover. This is the SETPOINT case the gains were tuned
%  against, so it is the direct check of the design objectives (overshoot and
%  settling time) on the full nonlinear model rather than the linear one.

ref_lat = traj.setpoint([0;0;3], [3;0;3], 35);
R_lat   = harness.scenario(ref_lat, config=P, label='Lateral step');

viz.tracking({R_lat}, {'Baseline'}, show='both', name='#A Lateral step');
viz.path3d({R_lat},  {'Baseline'}, name='#A Lateral step 3D');
harness.report({R_lat}, {'Lateral step'}, title='Lateral step (1 m in x)');

%  THE TUNING OBJECTIVES, CHECKED ON THE NONLINEAR MODEL.
%  Targets are the same ones the position loop was shaped to in
%  run_task2_frequency.m, so the linear design and the nonlinear check are
%  being judged by identical criteria.
OS_TARGET = 35;    % [%] maximum overshoot
TS_TARGET = 15;    % [s] maximum settling time, +-2% band

[~, S_lat] = viz.step_response(R_lat, 1, os_target=OS_TARGET, ts_target=TS_TARGET, ...
                               label='Lateral step (1 m in x)', ...
                               name='#A Lateral step metrics');

viz.actuators({R_lat}, {'Lateral step'}, P.UAV, name='#A Lateral step actuators');
viz.attitude({R_lat},  {'Lateral step'},        name='#A Lateral step attitude');

%% ---- B. Vertical step ----------------------------------------------
%  1 m step in z. .
ref_ver = traj.setpoint([0;0;1], [0;0;2], 100);
R_ver   = harness.scenario(ref_ver, config=P, label='Vertical step');

viz.tracking({R_ver}, {'Baseline'}, show='both', name='#B Vertical step');
harness.report({R_ver}, {'Vertical step'}, title='Vertical step (1 m in z)');


[~, S_ver] = viz.step_response(R_ver, 3, os_target=OS_TARGET, ts_target=TS_TARGET, ...
                               label='Vertical step (1 m in z)', ...
                               name='#B Vertical step metrics');

%  Side-by-side, since the cascade cost is the point being made
fprintf('\n=== Setpoint objectives: lateral vs vertical ===\n');
fprintf('%-26s %10s %10s %10s %10s\n', ...
        'case','rise [s]','OS [%]','Ts [s]','e_ss [m]');
fprintf('%s\n', repmat('-',1,70));
fprintf('%-26s %10.2f %10.2f %10.2f %10.4f\n', 'Lateral (x, via attitude)', ...
        S_lat.RiseTime, S_lat.Overshoot, S_lat.SettlingTime, S_lat.SSError);
fprintf('%-26s %10.2f %10.2f %10.2f %10.4f\n', 'Vertical (z, collective)', ...
        S_ver.RiseTime, S_ver.Overshoot, S_ver.SettlingTime, S_ver.SSError);
fprintf('targets: OS <= %g %%, Ts <= %g s\n', OS_TARGET, TS_TARGET);

%% ---- C. Helicoidal tracking, nominal --------------------------------
%  Trajectory tracking, the case Task 2.3.1 asks for. A helix is chosen
%  because it excites all three translational axes SIMULTANEOUSLY and keeps
%  the heading turning, so the yaw loop is exercised too - while staying at
%  0.71 m/s, well inside the 4-5 m/s the exam allows.
%
%  Unlike a step, the reference has persistent non-zero velocity and
%  acceleration, so any residual error here is TRACKING LAG (a bandwidth
%  limitation), not a transient. That is why the metric of merit is the RMS
%  error over the whole run rather than settling time.
ref_hel = traj.helix(150);
viz.kinematics(ref_hel, name='Helicoidal reference');

R_hel = harness.scenario(ref_hel, config=P, label='Helix nominal');

viz.tracking({R_hel}, {'Baseline'}, show='both', name='#C Helix tracking');
viz.axis_errors({R_hel}, {'Baseline'}, name='#C Helix per-axis errors');
viz.errors({R_hel},   {'Baseline'}, name='#C Helix error norms');
viz.path3d({R_hel},   {'Baseline'}, name='#C Helix 3D');
harness.report({R_hel}, {'Helix nominal'}, title='Helicoidal tracking (nominal)');

%% ---- D. Reduced propeller effectiveness  [Task 2.3.2] ---------------
%  Each rotor r now delivers only lambda_r of its commanded thrust. The
%  baseline controller is NOT told: its allocation still inverts the nominal
%  F2, so the wrench it actually produces is
%
%     w_actual = F1 * F2 * diag(lambda) * pinv(F2) * u   ~=   F1*u
%

lam_deg = disturb.fault(ref_hel.time, t0=0, lambda=[0.8; 0.6; 0.75; 0.78; 0.5; 0.9]);
R_deg   = harness.scenario(ref_hel, config=P, lambda=lam_deg, ...
                           label='Helix, degraded rotors');

viz.axis_errors({R_hel, R_deg}, {'Nominal','Reduced effectiveness'}, ...
                name='#D Reduced effectiveness, per-axis errors');
viz.errors({R_hel, R_deg}, {'Nominal','Reduced effectiveness'}, ...
           name='#D Effect of reduced propeller effectiveness');
viz.path3d({R_hel, R_deg}, {'Nominal','Reduced effectiveness'}, ...
           name='#D Degraded rotors 3D');

viz.actuators({R_hel, R_deg}, {'Nominal','Reduced effectiveness'}, P.UAV, ...
              name='#D Degraded rotors actuator usage');
harness.report({R_hel, R_deg}, {'Nominal','Degraded'}, ...
               title='Impact of reduced propeller effectiveness');

%% ---- E. Wind gust --------------------------------------------s-------
%  A 1-cosine ("1 minus cosine") gust, the standard discrete gust shape from
%  airworthiness practice:
%
%     V(tau) = 0.5*V_peak*(1 - cos(2*pi*tau/T_gust)),   0 <= tau <= T_gust
%
%  It rises smoothly, peaks at mid-window and returns to zero, so the response
%  is a property of the CONTROLLER and not an artefact of a discontinuous
%  input. It enters as a wind VELOCITY, so the force the airframe feels is
%  built from the relative airspeed through D_fv - physically meaningful, and
%  it means faster flight through the same gust produces a different load.
GUST = struct('t0', 60, 'duration', 20, 'peak', [8;0;0]);
wind = disturb.gust(ref_hel.time, t0=GUST.t0, duration=GUST.duration, peak=GUST.peak);

R_gust = harness.scenario(ref_hel, config=P, wind=wind, label='Helix + gust');

gust_window = [GUST.t0, GUST.t0 + GUST.duration];
viz.errors({R_hel, R_gust}, {'No gust','With gust'}, ...
           events=gust_window, eventNames={'Wind gust'}, name='#E Gust response');
viz.tracking({R_hel, R_gust}, {'No gust','With gust'}, ...
             events=gust_window, name='#E Gust tracking');
viz.tracking({R_hel, R_gust}, {'No gust','With gust'}, ...
             overlay=wind, overlayLabel='V_{gust}  [m/s]', ...
             name='#E Gust tracking with gust overlay');
viz.axis_errors({R_hel, R_gust}, {'No gust','With gust'}, events=gust_window, ...
                name='#E Gust per-axis errors');
harness.report({R_hel, R_gust}, {'No gust','With gust'}, ...
               window=gust_window, title='Wind gust response');

viz.attitude({R_hel, R_gust}, {'No gust','With gust'}, events=gust_window, ...
             name='#E Gust attitude response');
viz.actuators({R_hel, R_gust}, {'No gust','With gust'}, P.UAV, events=gust_window, ...
              name='#E Gust actuator usage');

% Gust shape, for the slide that documents the disturbance model
figure('Name','#E Gust profile','Color','w','Position',[140 140 700 400]);
plot(wind.time, vecnorm(wind.signals.values,2,2), 'LineWidth', 2.6); grid on
xlabel('t  [s]'); ylabel('|V_{wind}|  [m/s]');
title(sprintf('1-cosine gust: %.0f m/s peak over %.0f s', ...
              norm(GUST.peak), GUST.duration));
xlim([GUST.t0-10, GUST.t0+GUST.duration+10]);

% Parameter table for the slide that documents the disturbance model
fprintf('\n=== Gust model parameters ===\n');
fprintf('  %-22s %s\n',            'shape',         '1-cosine');
fprintf('  %-22s %.1f m/s\n',      'peak velocity', norm(GUST.peak));
fprintf('  %-22s [%.0f %.0f %.0f]\n', 'direction',  GUST.peak/max(norm(GUST.peak),eps));
fprintf('  %-22s %.1f s\n',        'duration',      GUST.duration);
fprintf('  %-22s %.1f s\n',        'onset',         GUST.t0);

%% ---- F. External force disturbance  [disturbance rejection] ---------
%  A sustained force applied directly to the airframe, entering the model on
%  the Disturbance_i channel: given in the INERTIAL frame, rotated into the
%  body by R' and added to the actuator wrench (zero torque), so it is an
%  external FORCE and not a wind velocity.
%

DIST = struct('t0', 60, 'value', [10;0;0]);   % [N] in the inertial frame
dist = disturb.step(ref_hel.time, t0=DIST.t0, value=DIST.value);

R_dist = harness.scenario(ref_hel, config=P, disturbance=dist, ...
                          label='Helix + force disturbance');

dist_window = [DIST.t0, ref_hel.time(end)];
viz.errors({R_hel, R_dist}, {'No disturbance','With force disturbance'}, ...
           events=dist_window, eventNames={'Force step'}, ...
           name='#F Force disturbance, error norms');
viz.tracking({R_hel, R_dist}, {'No disturbance','With force disturbance'}, ...
             events=dist_window, show='both', name='#F Force disturbance tracking');
viz.axis_errors({R_hel, R_dist}, {'No disturbance','With force disturbance'}, ...
                events=dist_window, name='#F Force disturbance, per-axis errors');

%  The tilt the vehicle has to hold, and what that tilt costs in rotor speed.
viz.attitude({R_hel, R_dist}, {'No disturbance','With force disturbance'}, ...
             events=dist_window, name='#F Force disturbance attitude response');
viz.actuators({R_hel, R_dist}, {'No disturbance','With force disturbance'}, P.UAV, ...
              events=dist_window, name='#F Force disturbance actuator usage');

harness.report({R_hel, R_dist}, {'No disturbance','With disturbance'}, ...
               window=dist_window, title='External force disturbance rejection');

%  Steady-state check, read AFTER the transient has died out.

settle = DIST.t0 + 30;
i_ep   = R_dist.t      >= settle;
i_ep0  = R_hel.t       >= settle;
i_fc   = R_dist.f_c_t  >= settle;
i_fc0  = R_hel.f_c_t   >= settle;

e_ss    = mean(vecnorm(R_dist.p(i_ep,:)  - R_dist.p_d(i_ep,:),  2, 2));
e_ss0   = mean(vecnorm(R_hel.p(i_ep0,:)  - R_hel.p_d(i_ep0,:),  2, 2));
df_c    = mean(R_dist.f_c(i_fc,:), 1) - mean(R_hel.f_c(i_fc0,:), 1);

fprintf('\n=== External force disturbance (%.1f N, inertial) ===\n', norm(DIST.value));
fprintf('  %-30s [%.3f %.3f %.3f] N\n', 'applied disturbance',    DIST.value);
fprintf('  %-30s [%.3f %.3f %.3f] N\n', 'change in commanded f_c', df_c);
fprintf('  %-30s %.4f m  (was %.4f m undisturbed)\n', 'mean |e_p| after settling', e_ss, e_ss0);

%  Tilt, as the physical picture: the vehicle leans into the disturbance, and
%  the lean is close to atan(F/mg) once the turn's own tilt is accounted for.
if isfield(R_dist, 'euler')
    tiltOf   = @(R, idx) rad2deg(mean(acos(cos(R.euler(idx,2)).*cos(R.euler(idx,3)))));
    tilt_act = tiltOf(R_dist, R_dist.euler_t >= settle);
    tilt_nom = tiltOf(R_hel,  R_hel.euler_t  >= settle);
    fprintf('  %-30s %.2f deg  (was %.2f deg undisturbed)\n', 'mean tilt after settling', tilt_act, tilt_nom);
    fprintf('  %-30s %.2f deg\n', 'atan(F/mg) for reference', ...
            rad2deg(atan2(norm(DIST.value), P.UAV.m*P.par.g)));
end

% Disturbance profile, for the slide that documents the disturbance model
figure('Name','#F Force disturbance profile','Color','w','Position',[140 140 700 400]);
stairs(dist.time, dist.signals.values(:,1), 'LineWidth', 2.6); grid on
xlabel('t  [s]'); ylabel('F_{d,x}  [N]');
title(sprintf('Inertial force step: %.1f N from t = %.0f s', ...
              norm(DIST.value), DIST.t0));
ylim([-0.5, norm(DIST.value)+0.5]);

%% ---- G. Hover at the origin, force step  [recovery time] -------------
%  The cleanest possible disturbance-rejection test: the setpoint is a fixed
%  point, [0 0 0]
%

HOV = struct('t0', 20, 'T', 120, 'amps', [3 6 10]);  % [N] step magnitudes, x
TOL = 0.05;                                          % [m] recovery band

ref_hov  = traj.setpoint([0;0;0], [0;0;0], HOV.T);
runs_hov = cell(1, numel(HOV.amps));
tags_hov = cell(1, numel(HOV.amps));
rec      = nan(numel(HOV.amps), 6);

for k = 1:numel(HOV.amps)
    F = HOV.amps(k);
    d = disturb.step(ref_hov.time, t0=HOV.t0, value=[F;0;0]);
    R = harness.scenario(ref_hov, config=P, disturbance=d, ...
                         label=sprintf('F_d = %g N', F));

    runs_hov{k} = R;
    tags_hov{k} = sprintf('F_d = %g N', F);

    nep   = vecnorm(R.p - R.p_d, 2, 2);
    after = find(R.t >= HOV.t0);

    [pk, ipk] = max(nep(after));
    t_pk      = R.t(after(ipk)) - HOV.t0;

    %  Recovery is the LAST exit from the band, not the first entry: a response
    %  that dips inside the band and swings back out has not recovered yet.
    last_out = find(nep(after) > TOL, 1, 'last');
    if isempty(last_out)
        t_rec = 0;                                  % never left the band
    elseif last_out == numel(after)
        t_rec = NaN;                                % still outside at the end
    else
        t_rec = R.t(after(last_out+1)) - HOV.t0;
    end

    tail = R.t >= HOV.t0 + 0.8*(R.t(end) - HOV.t0);
    e_ss = mean(nep(tail));

    tilt = NaN;
    if isfield(R, 'euler')
        e_tail = R.euler_t >= HOV.t0 + 0.8*(R.t(end) - HOV.t0);
        tilt   = rad2deg(mean(acos(cos(R.euler(e_tail,2)).*cos(R.euler(e_tail,3)))));
    end

    rec(k,:) = [F, pk, t_pk, t_rec, e_ss, tilt];
end

hov_window = [HOV.t0, HOV.T];
viz.errors(runs_hov, tags_hov, events=hov_window, eventNames={'Force step'}, ...
           name='#G Hover recovery, error norms');
viz.tracking(runs_hov, tags_hov, events=hov_window, ...
             name='#G Hover recovery tracking');
viz.axis_errors(runs_hov, tags_hov, events=hov_window, ...
                name='#G Hover recovery, per-axis errors');
viz.attitude(runs_hov, tags_hov, events=hov_window, ...
             name='#G Hover recovery attitude response');
viz.actuators(runs_hov, tags_hov, P.UAV, events=hov_window, ...
              name='#G Hover recovery actuator usage');
harness.report(runs_hov, tags_hov, window=hov_window, ...
               title='Hover recovery from an external force step');


s_hov = viz.style();
figure('Name','#G Hover recovery from force step','Color','w','Position',[140 140 860 440]);
hold on; grid on
for k = 1:numel(runs_hov)
    R = runs_hov{k};
    plot(R.t - HOV.t0, vecnorm(R.p - R.p_d, 2, 2), 'LineWidth', 2.1, ...
         'Color', s_hov.colors{k}, 'DisplayName', tags_hov{k});
    if isfinite(rec(k,4))
        xline(rec(k,4), '--', sprintf('%.1f s', rec(k,4)), ...
              'Color', s_hov.colors{k}, 'HandleVisibility','off');
    end
end
yline(TOL, ':k', sprintf('%.0f cm band', 100*TOL), 'HandleVisibility','off');
xlabel('time since the force step  [s]'); ylabel('|e_p|  [m]');
title('Recovery from a step force while hovering at the origin');
legend('Location','best'); xlim([-5, 60]);

fprintf('\n=== Hover at [0 0 0], recovery from a step force in x ===\n');
fprintf('%-8s %10s %10s %12s %11s %9s\n', ...
        'F_d [N]','peak [m]','t_peak [s]','t_rec [s]','e_ss [m]','tilt [deg]');
fprintf('%s\n', repmat('-',1,64));
for k = 1:size(rec,1)
    fprintf('%-8.1f %10.4f %10.2f %12.2f %11.4f %9.2f\n', rec(k,:));
end
fprintf('recovery = back inside a %.0f cm band and staying there; ', 100*TOL);
fprintf('reference tilt atan(F/mg): ');
fprintf('%.1f ', rad2deg(atan2(HOV.amps, P.UAV.m*P.par.g)));
fprintf('deg\n');

fprintf('\nTask 2 time-domain scenarios complete.\n');

%% ---- Save the figures ------------------------------------------------
% Written to Figures/Task2_setpoints/ , named from each figure's Name property and
% numbered so the folder sorts in slide order.
viz.export('Task2_setpoints');
