%% ========================================================================
%  TASK 3.1 / 3.3 - PACKAGE EFFECT AND PREDICTOR-BASED MRAC
%  ------------------------------------------------------------------------
%  3.1  THE PACKAGE AS A PLANT PERTURBATION
%  A rigid package of mass m_p is attached at r_c from the body origin. It
%  changes the plant in three ways at once, and only the first is obvious:
%
%    (a) mass:     m -> m + m_p
%    (b) CoM:      r_bg moves to (m*r_bg + m_p*r_c)/(m+m_p), so the static
%                  moment S = m*skew(r_bg) is no longer ZERO. The generalized
%                  mass matrix M = [m*I, S'; S, J] therefore stops being block
%                  diagonal: translation and rotation are now COUPLED.
%    (c) inertia:  J -> J + J_pkg + m_p*((r_c'r_c)I - r_c r_c')  (parallel axis)
%
%  The baseline controller knows none of this - the exam says the package is
%  "completely unknown". So it flies with the wrong mass, ignores the coupling
%  and carries an unmodelled weight m_p*g. See uav.with_package.
%
%  3.2/3.3  THE AUGMENTATION
%  Approximate uncertain translational dynamics used for the design:
%
%     p_dot = v
%     v_dot = -g*e_3 + (1/m)*( Lambda_f * u + f_e )
%
%  with Lambda_f a diagonal CONTROL EFFECTIVENESS matrix and f_e a constant
%  exogenous force. Note what each one absorbs physically:
%     f_e       <- the package WEIGHT and any steady wind (a force bias)
%     Lambda_f  <- the extra INERTIA, which scales commanded force into
%                  actual acceleration, and any rotor thrust loss
%
%  PREDICTOR (a copy of the assumed dynamics, driven by the state error):
%     p_hat_dot = v_hat + L_p*(p - p_hat)
%     v_hat_dot = -g*e_3 + (1/m)*( Lambda_hat*u + theta_hat ) + L_v*(v - v_hat)
%
%  UPDATE LAWS (driven by the PREDICTION error, not the tracking error - this
%  is what makes the augmentation predictor-based and separates the adaptation
%  dynamics from the reference-tracking error):
%
%     s = (v - v_hat)/m
%
%     theta_hat_dot =
%         Gamma_theta * ( s - sigma_e*norm(s)*theta_hat )
%
%     Lambda_hat_dot =
%         Proj( Lambda_hat, Gamma_lambda .* u_h .* s )
%
%  u_h = u*min(1, T_MAX/||u||) is the ACHIEVABLE force (saturation hedging):
%  both the predictor and the Lambda law use it instead of the command u.
%
%  ROBUSTIFICATION:
%    - theta_hat uses e-modification. The additional damping is proportional
%      to the prediction-error magnitude, limiting parameter drift in the
%      presence of bounded disturbances while vanishing as the prediction
%      error tends to zero.
%
%    - Lambda_hat uses a projection operator. The nominal update is retained
%      inside the admissible set and modified only near its boundary so that
%      the effectiveness estimate remains bounded. This is particularly
%      important because Lambda_hat is inverted in the augmented control law.
%
%  No sigma-modification is used in the final implementation. Neither
%  e-modification nor projection should be interpreted as guaranteeing
%  convergence of the adaptive states to physical parameters in the absence
%  of sufficient excitation.
%
%  AUGMENTED INPUT (cancel what has been estimated):
%     u_a = (Lambda_hat^-1 - I)*u_b - Lambda_hat^-1 * theta_hat
%
%  Why only position is augmented (the exam's Remark): the attitude loop is
%  faster (5.16 vs 1.66 rad/s bandwidth, ~3x, and it settles much faster) and
%  the geometric attitude controller is robust to the package's inertia change,
%  so the attitude error decays before it perturbs the outer loop.
%
%  ------------------------------------------------------------------------
%  SECTIONS
%    A [3.1] baseline, with vs without package - how bad is the problem
%    B [3.3] with package, baseline vs PBMRAC - does the augmentation fix it
%    C [3.3] payload sweep 0.5/0.6/0.7 kg + the Lambda_hat estimate
%    D [3.3] package + persistent asymmetric rotor fault (500 s)
%    E [3.3] package + rotor fault + a short wind gust
%  All fly the SAME helix, so any difference is the payload or the controller.
%  ========================================================================

clear; close all; clc;

%ref = traj.helix(500);
viz.defaults();   % Times New Roman, heavier lines (see +viz/defaults.m)

T_TEST = 300;      % final Task 3 simulation horizon [s]

HELIX_RADIUS   = 10;
HELIX_RATE     = 0.2;
HELIX_CLIMB    = 2;

HELIX_CENTER   = [10; 10];
% HELIX_PHASE0 removed: traj.helix now searches the entry phase itself (free-entry minimum snap).
APPROACH_TIME  = 12;

ref = traj.helix(T_TEST, ...
    radius        = HELIX_RADIUS, ...
    rate          = HELIX_RATE, ...
    climb         = HELIX_CLIMB, ...
    center        = HELIX_CENTER, ...
    start_p       = [0;0;0], ...
    approach_time = APPROACH_TIME);
%% ---- A. Effect of the package on the baseline  [Task 3.1] ----------
%  Same controller, same trajectory, package toggled on the PLANT only. The
%  gap between the two runs is the cost of the unmodelled payload.
%
%  Expect the steady-state error to grow least: the baseline's integral term
%  does eventually absorb the constant m_p*g bias. The transient and RMS
%  errors grow more, because integral action is slow and cannot compensate the
%  changed inertia at all.
P = harness.config(package_mass = 0.4);

R_nopkg = harness.scenario(ref, config=P, package=0, mrac=0, label='No package');
R_pkg   = harness.scenario(ref, config=P, package=1, mrac=0, label='With package');

viz.axis_errors({R_nopkg, R_pkg}, {'No package','With package'}, ...
                name='#A Package effect, per-axis errors');
viz.errors({R_nopkg, R_pkg}, {'No package','With package'}, ...
           name='#A Package effect, error norms');
viz.tracking({R_nopkg, R_pkg}, {'No package','With package'}, ...
             show='both', name='#A Package effect, tracking');
viz.path3d({R_nopkg, R_pkg}, {'No package','With package'}, ...
           name='#A Package effect, 3D');
harness.report({R_nopkg, R_pkg}, {'No package','With package'}, ...
               title='Task 3.1 - effect of an unknown package on the baseline');

%% Initial vertical-transient diagnostic

runs = {R_nopkg, R_pkg};
names = {'No package','With package'};
masses = [P.UAV.m, P.UAV.m + 0.4];

fprintf('\n=== INITIAL VERTICAL TRANSIENT ===\n');

for k = 1:2
    R = runs{k};

    idx = R.t <= 3;

    [vzmin, j] = min(R.v(idx,3));
    tt = R.t(idx);

    fprintf('\n%s\n', names{k});
    fprintf('  min vz in first 3 s = %+8.4f m/s at t = %.3f s\n', ...
        vzmin, tt(j));

    [zmin, jz] = min(R.p(idx,3));
    fprintf('  min z  in first 3 s = %+8.4f m   at t = %.3f s\n', ...
        zmin, tt(jz));

    if isfield(R,'T_out')
        it = R.T_out_t <= 3;

        % Signed rotor thrusts are negative for upward collective
        Tup = -sum(R.T_out(it,:),2);

        fprintf('  required hover thrust = %.4f N\n', ...
            masses(k)*9.81);
        fprintf('  actual upward thrust: min %.4f N, max %.4f N\n', ...
            min(Tup), max(Tup));

        fprintf('  initial upward thrust = %.4f N\n', Tup(1));
    end
end
%% ---- B. Adaptive augmentation with the package  [Task 3.3] ---------
%  Identical plant and reference; the only change is MRAC_on, which switches
%  the augmented input u_b + u_a in ahead of the allocation. The improvement
%  here is attributable to the adaptation alone.
R_bl = harness.scenario(ref, config=P, package=1, mrac=0, label='Baseline');
R_ad = harness.scenario(ref, config=P, package=1, mrac=1, robust=1, label='PBMRAC');

viz.axis_errors({R_bl, R_ad}, {'Baseline','PBMRAC'}, ...
                name='#B Adaptive augmentation, per-axis errors');
viz.errors({R_bl, R_ad}, {'Baseline','PBMRAC'}, ...
           name='#B Adaptive augmentation, error norms');
viz.tracking({R_bl, R_ad}, {'Baseline','PBMRAC'}, show='both', ...
             name='#B Adaptive augmentation, tracking');
viz.path3d({R_bl, R_ad}, {'Baseline','PBMRAC'}, ...
           name='#B Adaptive augmentation, 3D');
harness.report({R_bl, R_ad}, {'Baseline','PBMRAC'}, ...
               title='Task 3.3 - baseline vs PBMRAC with the package attached');

%  Where the augmentation's extra force actually goes. u_a raises the
%  collective to carry the unmodelled package weight, so T_c should sit closer
%  to the TRUE weight (m + m_p)*g than the baseline's does; the attitude
%  figure confirms the augmentation is not buying its accuracy with a more
%  aggressive tilt.
viz.actuators({R_bl, R_ad}, {'Baseline','PBMRAC'}, harness.config().UAV, ...
              name='#B Adaptive augmentation, actuator usage');
viz.attitude({R_bl, R_ad}, {'Baseline','PBMRAC'}, ...
             name='#B Adaptive augmentation, attitude');
%% IMPORTANT CHECK
p = traj.pos(ref);
v = traj.vel(ref);
a = traj.acc(ref);
t = ref.time(:);

Ta = APPROACH_TIME;

iL = find(t <  Ta, 1, 'last');
iR = find(t >= Ta, 1, 'first');

fprintf('\n=== HELIX JOIN CONTINUITY ===\n');
fprintf('tL = %.6f s,  tR = %.6f s\n', t(iL), t(iR));

fprintf('\nPosition jump [m]:\n');
disp(p(:,iR) - p(:,iL))

fprintf('Velocity jump [m/s]:\n');
disp(v(:,iR) - v(:,iL))

fprintf('Acceleration jump [m/s^2]:\n');
disp(a(:,iR) - a(:,iL))

fprintf('Norms:\n');
fprintf('|Dp| = %.6g m\n',     norm(p(:,iR)-p(:,iL)));
fprintf('|Dv| = %.6g m/s\n',   norm(v(:,iR)-v(:,iL)));
fprintf('|Da| = %.6g m/s^2\n', norm(a(:,iR)-a(:,iL)));
%%
mb = harness.metrics(R_bl);
ma = harness.metrics(R_ad);

windows = [ ...
      0    APPROACH_TIME;
      APPROACH_TIME  APPROACH_TIME+40;
    100    T_TEST];

names = {'approach','join transient','mature helix'};

fprintf('\n=== PBMRAC PER-AXIS CHECK ===\n');

for q = 1:size(windows,1)

    ib = mb.t >= windows(q,1) & mb.t <= windows(q,2);
    ia = ma.t >= windows(q,1) & ma.t <= windows(q,2);

    rb = sqrt(mean(mb.ep(ib,:).^2,1));
    ra = sqrt(mean(ma.ep(ia,:).^2,1));

    fprintf('\n%s  [%.0f %.0f] s\n', ...
        names{q}, windows(q,1), windows(q,2));

    fprintf('             x          y          z\n');
    fprintf('Baseline  %8.4f   %8.4f   %8.4f\n', rb);
    fprintf('PBMRAC    %8.4f   %8.4f   %8.4f\n', ra);
end
%%
p = traj.pos(ref);
v = traj.vel(ref);
a = traj.acc(ref);

fprintf('Initial reference:\n');
fprintf('|v(0)| = %.6g m/s\n', norm(v(:,1)));
fprintf('|a(0)| = %.6g m/s^2\n', norm(a(:,1)));
disp('v(0) ='); disp(v(:,1))
disp('a(0) ='); disp(a(:,1))

%% Lift-off diagnostic
runs  = {R_nopkg, R_pkg, R_bl, R_ad};
names = {'No package','Package baseline','Baseline','PBMRAC'};

z_eps  = 1e-3;    % 1 mm
vz_eps = 1e-3;

fprintf('\n=== LIFT-OFF DIAGNOSTIC ===\n');

for k = 1:numel(runs)
    R = runs{k};

    % first point clearly above the floor and moving upward
    idx = find(R.p(:,3) > z_eps & R.v(:,3) > vz_eps, 1, 'first');

    fprintf('\n%s\n', names{k});

    if isempty(idx)
        fprintf('  No lift-off detected.\n');
    else
        fprintf('  lift-off time : %.4f s\n', R.t(idx));
        fprintf('  z at lift-off : %.6f m\n', R.p(idx,3));
        fprintf('  vz at lift-off: %.6f m/s\n', R.v(idx,3));

        % tracking error at that instant
        ep = R.p(idx,:) - R.p_d(idx,:);
        fprintf('  position error at lift-off = %.4f m\n', norm(ep));
        fprintf('  desired z at lift-off      = %.4f m\n', R.p_d(idx,3));
    end
end
%% ---- C. Payload sweep ----------------------------------------------
%  Three payloads, each flown by both controllers, on a LONGER horizon (250 s)
%  because the estimates need time to move.
%
%  The exam asks: what is the effect of the added package mass on the control
%  effectiveness? Read it off Lambda_hat. The augmentation's model says
%  acceleration = Lambda_f*u/m with m the NOMINAL mass; the true vehicle has
%  m + m_p, so the same commanded force produces less acceleration and the
%  estimator explains it as Lambda_f < 1. The heavier the payload, the deeper
%  the dip - which is what the printed table checks quantitatively.
masses  = [0.5 0.6 0.7];
T_sweep = 500;                       % longer horizon: the estimates need time
%ref_s   = traj.helix(T_sweep);
ref_s = traj.helix(T_sweep, ...
    radius        = HELIX_RADIUS, ...
    rate          = HELIX_RATE, ...
    climb         = HELIX_CLIMB, ...
    center        = HELIX_CENTER, ...
    start_p       = [0;0;0], ...
    approach_time = APPROACH_TIME);
% Six runs back to back: keep the model compiled between them (~2.3x faster,
% and verified to give bit-identical results -- see harness.fastrestart).
% Deliberately unread: onCleanup handle, released by the clear below.
frGuard = harness.fastrestart('Final'); %#ok<NASGU>

sweep = struct('m', {}, 'ad', {}); %, 'bl', {}
for j = 1:numel(masses)
    Pj = harness.config(package_mass = masses(j));
    sweep(j).m  = masses(j);
    sweep(j).bl = harness.scenario(ref_s, config=Pj, package=1, mrac=0);
    sweep(j).ad = harness.scenario(ref_s, config=Pj, package=1, mrac=1,robust =1);
    fprintf('payload %.1f kg done\n', masses(j));
end
clear frGuard      % release fast restart before plotting

% Position error norm vs payload
figure('Name','#C Error vs payload','Color','w','Position',[100 100 1200 380]);
for j = 1:numel(masses)
    subplot(1, numel(masses), j)
    mb = harness.metrics(sweep(j).bl);
    ma = harness.metrics(sweep(j).ad);
    semilogy(mb.t, mb.nep, 'LineWidth',1.5); hold on
    semilogy(ma.t, ma.nep, 'LineWidth',1.5); grid on
    xlabel('t  [s]'); ylabel('|e_p|  [m]');
    title(sprintf('m_p = %.1f kg', masses(j)));
    if j == 1, legend('Baseline','PBMRAC','Location','best'); end
end

% Estimated control effectiveness vs payload
if isfield(sweep(1).ad, 'lambda_hat')
    figure('Name','#C Effectiveness vs payload','Color','w','Position',[120 120 1200 380]);
    for j = 1:numel(masses)
        subplot(1, numel(masses), j)
        L = sweep(j).ad.lambda_hat;
        plot(sweep(j).ad.lambda_hat_t, L, 'LineWidth',1.5); grid on
        xlabel('t  [s]'); ylabel('$\hat\Lambda_f$','Interpreter','latex');
        title(sprintf('m_p = %.1f kg', masses(j)));
        if j == 1
            legend('\Lambda_{f,x}','\Lambda_{f,y}','\Lambda_{f,z}','Location','best');
        end
    end

    % Report the minimum, maximum and final values without interpreting
    % Lambda_hat as a direct physical identification of the payload. In the
    % final robustness implementation Lambda_hat is bounded by projection only;
    % no leakage term drives it back towards the nominal value. Since the helix
    % is only weakly exciting, convergence of the estimate to a unique physical
    % effectiveness value is not expected. The observed excursions are therefore
    % used to compare the adaptive response across payload masses.
    % fprintf('\n%-10s %12s %12s %12s %14s\n', ...
    %         'payload','min Lam_x','min Lam_y','min Lam_z','Lam_z(end)');
    % for j = 1:numel(masses)
    %     L = sweep(j).ad.lambda_hat;
    %     fprintf('%-10.1f %12.4f %12.4f %12.4f %14.4f\n', masses(j), min(L,[],1), L(end,3));
    % end
    fprintf('\n%-8s  %-27s  %-27s  %-27s\n', ...
    'payload', 'Lambda_x [min max final]', ...
               'Lambda_y [min max final]', ...
               'Lambda_z [min max final]');

    for j = 1:numel(masses)
    
        L = sweep(j).ad.lambda_hat;
    
        fprintf('%-8.1f  [%6.4f %6.4f %6.4f]  ', ...
            masses(j), min(L(:,1)), max(L(:,1)), L(end,1));
    
        fprintf('[%6.4f %6.4f %6.4f]  ', ...
            min(L(:,2)), max(L(:,2)), L(end,2));
    
        fprintf('[%6.4f %6.4f %6.4f]\n', ...
            min(L(:,3)), max(L(:,3)), L(end,3));
    end
    fprintf([ ...
      '\nReading: the added mass reaches the augmentation through TWO channels.\n' ...
      '  * Its weight is a constant force offset, which theta_hat absorbs -- that\n' ...
      '    is the channel that removes the steady-state error, and it is why the\n' ...
      '    steady-state metric improves most.\n' ...
      '  * The extra inertia scales the commanded-force-to-acceleration gain, so\n' ...
      '    Lambda_hat dips below 1 while the vehicle is accelerating, more deeply\n' ...
      '    for a heavier payload.\n' ...
      'Because the helix is a gentle, weakly exciting manoeuvre, the adaptive\n' ...
      'states are not interpreted as identified physical parameters. Lambda_hat\n' ...
      'is kept inside its admissible set by projection; no sigma leakage towards\n' ...
      'the nominal value is used. The additive theta_hat estimate instead uses\n' ...
      'e-modification, whose damping scales with the prediction error. The trends\n' ...
      'are therefore used to assess adaptive behaviour, not parameter convergence.\n']);
else
    warning(['lambda_hat was not logged. Enable signal logging on the ' ...
             '''lambda hat'' output in Final.slx to produce the Task 3.3 figure.']);
end

% Summary across the sweep
fprintf('\n%-10s %12s %12s %12s\n','payload','RMS baseline','RMS PBMRAC','improvement');
for j = 1:numel(masses)
    mb = harness.metrics(sweep(j).bl);
    ma = harness.metrics(sweep(j).ad);
    fprintf('%-10.1f %12.4f %12.4f %11.1f%%\n', masses(j), mb.rms_p, ma.rms_p, ...
            100*(1 - ma.rms_p/mb.rms_p));
end

clear sweep frGuard

%% ---- D. Package AND degraded rotors, together ----------------------
%  The hardest case in Task 3, and the one that motivates the augmentation's
%  Lambda_f term rather than just its f_e term:
%
%     * the package adds unknown MASS  -> a constant force bias (weight) and a
%       gain change (more inertia per unit commanded force);
%     * the rotors lose effectiveness  -> a further gain change, this time
%       asymmetric across the six arms.
%
%  Both degrade the same channel, so the baseline's integral action - which can
%  eventually absorb a constant bias - cannot fix either gain error. The
%  adaptive law can, because Lambda_hat multiplies the input rather than adding
%  to it.
%
%  NOTE ON SHAPE: lambda_sig is a 6-channel SIGNAL (one per rotor), not a
%  matrix, so the fault is given as the vector [0.8 0.6 1 1 0.4 1]. Rotors 2
%  and 5 sit on opposing arms, so the loss is asymmetric and produces a
%  parasitic moment as well as a thrust deficit.
%
%  The fault is applied as a STEP part-way through rather than from t = 0. That
%  is deliberate: it lets the same run show Lambda_hat before the fault (where
%  it reflects the payload alone) and after (payload + rotor loss), which is
%  what makes the estimate legible.

T_DE = 500;
FAULT_AFTER_JOIN = 50;


T_FAULT = APPROACH_TIME + FAULT_AFTER_JOIN;
%ref_DE = traj.helix(T_DE);
ref_DE = traj.helix(T_DE, ...
    radius        = HELIX_RADIUS, ...
    rate          = HELIX_RATE, ...
    climb         = HELIX_CLIMB, ...
    center        = HELIX_CENTER, ...
    start_p       = [0;0;0], ...
    approach_time = APPROACH_TIME);
LAM_FAULT = [0.8; 0.6; 0.9; 1; 0.8; 0.9];
%T_FAULT   = 60;                       % [s] onset, mid-run

lam_deg = disturb.fault(ref_DE.time, t0 = T_FAULT, lambda = LAM_FAULT);

R_bl_deg = harness.scenario(ref_DE, config=P, package=1, mrac=0, lambda=lam_deg, ...
                            label='Baseline, pkg + fault');
R_ad_deg = harness.scenario(ref_DE, config=P, package=1, mrac=1, robust =1, lambda=lam_deg, ...
                            label='PBMRAC, pkg + fault');

runs_D = {R_bl_deg, R_ad_deg};
tags_D = {'Baseline','PBMRAC'};
evt_D  = [T_FAULT, ref_DE.time(end)];

% THE headline figure for this scenario: what the augmentation estimates
viz.effectiveness({R_ad_deg}, {'PBMRAC'}, lambda = lam_deg, t_fault = T_FAULT, ...
                  name = '#D Control effectiveness, package + rotor fault');

viz.errors(runs_D, tags_D, events = evt_D, eventNames = {'Rotor fault'}, ...
           name = '#D Package + fault, error norms');
viz.axis_errors(runs_D, tags_D, events = evt_D, ...
                name = '#D Package + fault, per-axis errors');
viz.actuators(runs_D, tags_D, P.UAV, events = evt_D, ...
              name = '#D Package + fault, actuator usage');

harness.report(runs_D, tags_D, window = [T_FAULT, ref_DE.time(end)], ...
               title = 'Task 3 - package AND reduced propeller effectiveness');

% How much of the damage is the fault, as opposed to the package alone?
fprintf('\n%-34s %12s %12s\n', 'case', 'RMS |e_p|', 'ss |e_p|');
fprintf('%s\n', repmat('-',1,60));
cases = {R_bl,     'Baseline, package only'
         R_bl_deg, 'Baseline, package + fault'
         R_ad,     'PBMRAC,   package only'
         R_ad_deg, 'PBMRAC,   package + fault'};
for c = 1:size(cases,1)
    mc = harness.metrics(cases{c,1});
    fprintf('%-34s %12.4f %12.4f\n', cases{c,2}, mc.rms_p, mc.ss_p);
end

L = R_ad_deg.lambda_hat;

fprintf('\n=== Lambda_hat observed range ===\n');
fprintf('axis          min          max        final\n');
fprintf('------------------------------------------------\n');

names = {'x','y','z'};

for k = 1:3
    fprintf('%s        %8.4f     %8.4f     %8.4f\n', ...
        names{k}, min(L(:,k)), max(L(:,k)), L(end,k));
end

%% ---- E. Package + asymmetric rotor degradation + short high-magnitude gust
%  Robustness stress case built directly on section D.
%
%  Nothing in the plant/controller is changed:
%    * same 0.4 kg package from P;
%    * same asymmetric rotor-effectiveness step from section D;
%    * same helix;
%    * one short 1-cosine gust is superimposed later in the run.
%
%  The fault remains active permanently after T_FAULT. The gust is temporary,
%  so the pre- and post-gust plant conditions are identical. This makes it
%  possible to see whether a transient unmodelled disturbance merely causes a
%  bounded adaptive excursion or leaves the adaptive states permanently shifted.


GUST_AFTER_JOIN  = 125;

T_GUST  = APPROACH_TIME + GUST_AFTER_JOIN;

%T_GUST = 140;                         % [s] well after the rotor fault
D_GUST = 20;                           % [s] short burst
V_GUST = [-10; 0;0 ];                   % [m/s] peak inertial wind velocity

wind_burst = disturb.gust(ref_DE.time, ...
    t0       = T_GUST, ...
    duration = D_GUST, ...
    peak     = V_GUST);

R_bl_burst = harness.scenario(ref_DE, ...
    config  = P, ...
    package = 1, ...
    mrac    = 0, ...
    robust  = 0,...
    lambda  = lam_deg, ...
    wind    = wind_burst, ...
    label   = 'Baseline, pkg + fault + gust');

R_ad_burst = harness.scenario(ref_DE, ...
    config  = P, ...
    package = 1, ...
    mrac    = 1, ...
    robust  = 1, ...
    lambda  = lam_deg, ...
    wind    = wind_burst, ...
    label   = 'PBMRAC, pkg + fault + gust');

runs_E = {R_bl_burst, R_ad_burst};
tags_E = {'Baseline','PBMRAC'};

evt_E = [T_FAULT, ref_DE.time(end);
         T_GUST,  T_GUST + D_GUST];

evt_names_E = {'Rotor fault active','Short gust'};

% ---- E1. Disturbance profile -----------------------------------------
figure('Name','#E Short high-magnitude gust','Color','w', ...
       'Position',[120 120 1000 430]);
plot(wind_burst.time, wind_burst.signals.values, 'LineWidth',1.4);
grid on
xline(T_GUST,'--','Gust starts','HandleVisibility','off');
xline(T_GUST + D_GUST,'--','Gust ends','HandleVisibility','off');
xlabel('t  [s]');
ylabel('V_{gust}  [m/s]');
title(sprintf('Short 1-cosine gust: %.0f m/s peak, %.0f s duration', ...
              norm(V_GUST), D_GUST));
legend('V_{w,x}','V_{w,y}','V_{w,z}','Location','best');

% ---- E2. Same graphical scheme as the previous Task-3 cases ----------
viz.tracking(runs_E, tags_E, ...
    events = evt_E, ...
    show = 'both', ...
    overlay = wind_burst, ...
    overlayLabel = 'V_{gust}  [m/s]', ...
    name = '#E Package + fault + gust, tracking');

viz.axis_errors(runs_E, tags_E, ...
    events = evt_E, ...
    name = '#E Package + fault + gust, per-axis errors');

viz.errors(runs_E, tags_E, ...
    events = evt_E, ...
    eventNames = evt_names_E, ...
    name = '#E Package + fault + gust, error norms');

viz.path3d(runs_E, tags_E, ...
    name = '#E Package + fault + gust, 3D');

viz.actuators(runs_E, tags_E, P.UAV, ...
    events = evt_E, ...
    name = '#E Package + fault + gust, actuator usage');

viz.attitude(runs_E, tags_E, ...
    events = evt_E, ...
    name = '#E Package + fault + gust, attitude');

viz.effectiveness({R_ad_burst}, {'PBMRAC'}, ...
    lambda = lam_deg, ...
    t_fault = T_FAULT, ...
    name = '#E Control effectiveness, package + fault + gust');

T_RECOVERY = 25;
report_end = min(ref_DE.time(end), T_GUST + D_GUST + T_RECOVERY);

harness.report(runs_E, tags_E, ...
    window = [T_GUST, report_end], ...
    title = 'Task 3E - package + asymmetric rotor fault + short gust');

% ---- E3. Adaptive-state response: section D versus section E ----------
axis_name = {'x','y','z'};

if isfield(R_ad_deg,'lambda_hat') && isfield(R_ad_burst,'lambda_hat')
    figure('Name','#E Lambda response to gust','Color','w', ...
           'Position',[140 100 1050 720]);
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    for i = 1:3
        nexttile
        plot(R_ad_deg.lambda_hat_t, R_ad_deg.lambda_hat(:,i), ...
             '--','Color',[0.45 0.45 0.45],'LineWidth',1.25, ...
             'DisplayName','Package + fault'); hold on
        plot(R_ad_burst.lambda_hat_t, R_ad_burst.lambda_hat(:,i), ...
             'm','LineWidth',1.45, ...
             'DisplayName','Package + fault + gust');

        xline(T_GUST,'--','HandleVisibility','off');
        xline(T_GUST + D_GUST,'--','HandleVisibility','off');
        grid on
        ylabel(sprintf('\hat{\lambda}_{%s}',axis_name{i}), ...
               'Interpreter','tex');

        if i == 1, legend('Location','best'); end
        if i == 3, xlabel('t  [s]'); end
    end

    title(tl,'Effect of the temporary gust on the effectiveness adaptive states');
end

theta_field = '';
if isfield(R_ad_burst,'theta_a_hat_pred') && isfield(R_ad_deg,'theta_a_hat_pred')
    theta_field = 'theta_a_hat_pred';
elseif isfield(R_ad_burst,'theta_hat') && isfield(R_ad_deg,'theta_hat')
    theta_field = 'theta_hat';
end

if ~isempty(theta_field)
    th_D = R_ad_deg.(theta_field);
    tt_D = R_ad_deg.([theta_field '_t']);
    th_E = R_ad_burst.(theta_field);
    tt_E = R_ad_burst.([theta_field '_t']);

    figure('Name','#E Additive adaptive-state response to gust','Color','w', ...
           'Position',[160 100 1050 720]);
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    for i = 1:3
        nexttile
        plot(tt_D, th_D(:,i), ...
             '--','Color',[0.45 0.45 0.45],'LineWidth',1.25, ...
             'DisplayName','Package + fault'); hold on
        plot(tt_E, th_E(:,i), ...
             'm','LineWidth',1.45, ...
             'DisplayName','Package + fault + gust');

        xline(T_GUST,'--','HandleVisibility','off');
        xline(T_GUST + D_GUST,'--','HandleVisibility','off');
        grid on
        ylabel(sprintf('\hat{\theta}_{%s}  [N]',axis_name{i}), ...
               'Interpreter','tex');

        if i == 1, legend('Location','best'); end
        if i == 3, xlabel('t  [s]'); end
    end

    title(tl,'Effect of the temporary gust on the additive adaptive states');
else
    warning(['theta adaptive state is not available in both D and E runs. ' ...
             'Log theta_a_hat_PRED or theta_hat to compare the gust response.']);
end

% ---- E4. Pre-gust / post-gust adaptive-state displacement ------------
PRE_WINDOW  = [T_GUST - 30, T_GUST - 5];
POST_WINDOW = [T_GUST + D_GUST + 20, T_GUST + D_GUST + 45];

fprintf('\n=== Task 3E - adaptive-state displacement across the gust ===\n');
fprintf('pre window : %.1f -> %.1f s\n', PRE_WINDOW);
fprintf('post window: %.1f -> %.1f s\n', POST_WINDOW);

if isfield(R_ad_burst,'lambda_hat')
    tL = R_ad_burst.lambda_hat_t;
    L  = R_ad_burst.lambda_hat;

    k_pre  = tL >= PRE_WINDOW(1)  & tL <= PRE_WINDOW(2);
    k_post = tL >= POST_WINDOW(1) & tL <= POST_WINDOW(2);

    L_pre  = mean(L(k_pre,:),1);
    L_post = mean(L(k_post,:),1);

    fprintf('\nLambda_hat mean shift:\n');
    fprintf('axis          pre         post        post-pre\n');
    for i = 1:3
        fprintf('%s        %9.4f    %9.4f    %+9.4f\n', ...
                axis_name{i}, L_pre(i), L_post(i), L_post(i)-L_pre(i));
    end
end

if ~isempty(theta_field)
    tth = R_ad_burst.([theta_field '_t']);
    th  = R_ad_burst.(theta_field);

    k_pre  = tth >= PRE_WINDOW(1)  & tth <= PRE_WINDOW(2);
    k_post = tth >= POST_WINDOW(1) & tth <= POST_WINDOW(2);

    th_pre  = mean(th(k_pre,:),1);
    th_post = mean(th(k_post,:),1);

    fprintf('\ntheta_hat mean shift [N]:\n');
    fprintf('axis          pre         post        post-pre\n');
    for i = 1:3
        fprintf('%s        %9.4f    %9.4f    %+9.4f\n', ...
                axis_name{i}, th_pre(i), th_post(i), th_post(i)-th_pre(i));
    end
end

fprintf('\nTask 3 scenarios complete.\n');

%% ---- Save the Task 3 scenario figures (A-E) ----------------------------
% Exported before the Monte Carlo, then closed, so Figures/Task3_final holds
% only the scenario figures and Figures/Task3_final_montecarlo only the
% campaign's.
viz.export('Task3_final');
close all;

%% ========================================================================
%  FINAL TASK 3.4 MONTE CARLO
%  ========================================================================
%  The refinement diagnostics above are intentionally disabled. The final
%  Monte Carlo campaign uses the selected controller configuration and the
%  final uncertainty model.
%  MC = run_montecarlo_final(n=500);
% % %% ---- Save the figures ------------------------------------------------
% % % Figures are exported using the figure Name properties and numbered in
% % % presentation order.
%  %viz.export('Task3_final_montecarlo');
%  %%
% %MC = run_montecarlo_final(n=500, L_p=15, L_v=15, Gamma_theta=diag([14.4 14.4 5.76]), Gamma_lambda=[0.075;0.075;0.015], label='montecarlo_final_scaled0p1', export='Task3_montecarlo_scaled0p1');
% %MC = run_montecarlo_final(resume='results\montecarlo_final_20260925_125048.mat');
% MC = run_montecarlo_final(n=500, Gamma_lambda=[0.75;0.75;1.5], label='montecarlo_final_glamz10', parallel=true);