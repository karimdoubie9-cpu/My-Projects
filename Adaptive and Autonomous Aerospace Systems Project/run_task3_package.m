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
%  UPDATE LAWS (driven by the PREDICTION error, not the tracking error - that
%  is what makes it predictor-based, and it decouples adaptation speed from
%  the reference):
%     theta_hat_dot  = Gamma_theta * (v - v_hat)
%     Lambda_hat_dot = Gamma_lambda .* u .* (v - v_hat)/m
%                      - sigma_lambda .* (Lambda_hat - 1)
%
%  The last term is SIGMA-MODIFICATION. Without persistent excitation the
%  estimates can drift on the manifold that fits the data equally well; the
%  leak term pulls them back towards the nominal and guarantees BOUNDEDNESS.
%  It does NOT give parameter convergence - which is exactly what section C
%  below shows happening.
%
%  AUGMENTED INPUT (cancel what has been estimated):
%     u_a = (Lambda_hat^-1 - I)*u_b - Lambda_hat^-1 * theta_hat
%
%  Why only position is augmented (the exam's Remark): the attitude loop is an
%  order of magnitude faster and the package's inertia change is small relative
%  to J, so the attitude error decays before it perturbs the outer loop.
%
%  ------------------------------------------------------------------------
%  SECTIONS
%    A [3.1] baseline, with vs without package - how bad is the problem
%    B [3.3] with package, baseline vs PBMRAC - does the augmentation fix it
%    C [3.3] payload sweep 0.2/0.3/0.4 kg + the Lambda_hat estimate
%  All fly the SAME helix, so any difference is the payload or the controller.
%  ========================================================================

clear; close all; clc;

%ref = traj.helix(500);
viz.defaults();   % Times New Roman, heavier lines (see +viz/defaults.m)

T_TEST = 500;       % or whatever final duration you want

HELIX_RADIUS   = 10;
HELIX_RATE     = 0.2;
HELIX_CLIMB    = 2;

HELIX_CENTER   = [10; 10];
HELIX_PHASE0   = 0;
APPROACH_TIME  = 20;

ref = traj.helix(T_TEST, ...
    radius        = HELIX_RADIUS, ...
    rate          = HELIX_RATE, ...
    climb         = HELIX_CLIMB, ...
    center        = HELIX_CENTER, ...
    phase0        = HELIX_PHASE0, ...
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
%%
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
    phase0        = HELIX_PHASE0, ...
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

    % Report the DEEPEST excursion, not the final value: the sigma-modification
    % term -sigma_lambda*(Lambda_hat - 1) biases the estimate back towards 1
    % whenever the excitation is not persistent, so Lambda_hat(end) merely
    % reports that bias (it is 1.000 for every payload) and says nothing about
    % the payload.  The transient minimum is the quantity that scales with mass.
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
      'Because the helix is a gentle, weakly exciting manoeuvre, the sigma-\n' ...
      'modification dominates once the transient decays and pulls Lambda_hat back\n' ...
      'to 1. That is the expected behaviour of sigma-mod without persistent\n' ...
      'excitation: it guarantees boundedness, not parameter convergence.\n']);
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
    phase0        = HELIX_PHASE0, ...
    start_p       = [0;0;0], ...
    approach_time = APPROACH_TIME);
LAM_FAULT = [0.8; 0.6; 0.75; 0.78; 0.5; 0.9];
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

fprintf('\n=== Lambda_hat unconstrained range ===\n');
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
%    * same 0.25 kg package from P;
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
V_GUST = [10; 0; 0];                   % [m/s] peak inertial wind velocity

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
%% NEW MONTECARLO TEST

MC = run_montecarlo_oldhelix( ...
    n        = 20, ...
    duration = 100, ...
    model    = 'Final');
%%
load('results/montecarlo_20260919_163739.mat','MC');

k = 9;

P9 = harness.config( ...
    package_mass = MC.samples.mass(k), ...
    package_r_c  = MC.samples.r_c(k,:).');

lam9 = disturb.fault( ...
    MC.ref.time, ...
    t0 = 0, ...
    lambda = MC.samples.lambda(k,:).');

R9 = harness.scenario( ...
    MC.ref, ...
    config  = P9, ...
    package = 1, ...
    mrac    = 1, ...
    robust  = 1, ...
    lambda  = lam9, ...
    model   = 'Final', ...
    label   = 'MC run 9 - robust PBMRAC');

viz.errors({R9},{'MC run 9'}, ...
    name='MC run 9 - tracking errors');

viz.actuators({R9},{'MC run 9'},P9.UAV, ...
    name='MC run 9 - actuators');

viz.attitude({R9},{'MC run 9'}, ...
    name='MC run 9 - attitude');

viz.effectiveness({R9},{'MC run 9'}, ...
    lambda=lam9, ...
    t_fault=0, ...
    name='MC run 9 - adaptive effectiveness');

harness.report({R9},{'MC run 9'}, ...
    title='Monte Carlo run 9 diagnostic');
%%
% Original asymmetric run 9
lam_asym = disturb.fault( ...
    MC.ref.time, ...
    t0 = 0, ...
    lambda = MC.samples.lambda(k,:).');

R9_asym = harness.scenario( ...
    MC.ref, ...
    config  = P9, ...
    package = 1, ...
    mrac    = 1, ...
    robust  = 1, ...
    lambda  = lam_asym, ...
    model   = 'Final', ...
    label   = 'Asymmetric');

% Same mean effectiveness, but identical on all six rotors
R9_uniform = harness.scenario( ...
    MC.ref, ...
    config  = P9, ...
    package = 1, ...
    mrac    = 1, ...
    robust  = 1, ...
    lambda  = lam_uniform, ...
    model   = 'Final', ...
    label   = 'Uniform');

runs = {R9_asym, R9_uniform};
tags = {'Asymmetric \lambda_r', 'Uniform mean \lambda_r'};

viz.errors(runs, tags, ...
    name='MC run 9 - asymmetric vs uniform');

viz.actuators(runs, tags, P9.UAV, ...
    name='MC run 9 - actuator comparison');

viz.attitude(runs, tags, ...
    name='MC run 9 - attitude comparison');

harness.report(runs, tags, ...
    title='MC run 9 - asymmetric vs uniform rotor effectiveness');

%% ========================================================================
%  MC RUN 9 - ASYMMETRIC VS UNIFORM ROTOR EFFECTIVENESS
%
%  Purpose:
%    1. Quantify the effective wrench mapping created by the six rotor
%       efficiencies:
%
%           w_actual = F2 * diag(lambda_r) * F2_pinv * w_cmd
%
%    2. Compare asymmetric rotor degradation against a uniform degradation
%       having exactly the same mean effectiveness.
%
%    3. Evaluate the pure-hover parasitic moments caused by asymmetry.
%
%    4. Check whether the degraded vehicle is physically capable of hover if
%       the true rotor effectiveness were known by the allocator.
%
%    5. Re-run the two nonlinear simulations with robust PBMRAC and compare:
%         - tracking
%         - actuator usage
%         - attitude
%         - adaptive Lambda estimates
%
%  IMPORTANT:
%    This script does NOT alter the controller or gains.
%    It is purely a diagnostic.
%  ========================================================================

clearvars;
close all;
clc;

% ------------------------------------------------------------------------
% USER SETTINGS
% -------------------------------------------------------------------------

MC_FILE = 'results/montecarlo_20260919_163739.mat';
MODEL   = 'Final';

k = 9;      % Monte Carlo sample to inspect

RUN_SIMULATIONS = true;


% ========================================================================
%  1. LOAD EXACT MONTE CARLO SAMPLE
%  ========================================================================

load(MC_FILE, 'MC');

mass_pkg   = MC.samples.mass(k);
r_c_pkg    = MC.samples.r_c(k,:).';
lambda_asym = MC.samples.lambda(k,:).';

lambda_mean = mean(lambda_asym);
lambda_uni  = lambda_mean * ones(6,1);

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' MC RUN %d - INPUT PARAMETERS\n', k);
fprintf('=====================================================================\n');

fprintf('Package mass           : %.6f kg\n', mass_pkg);
fprintf('Package CoM            : [%+.6f  %+.6f  %+.6f] m\n', ...
    r_c_pkg(1), r_c_pkg(2), r_c_pkg(3));

fprintf('\nAsymmetric rotor effectiveness:\n');
fprintf('  lambda = [');
fprintf(' %.6f', lambda_asym);
fprintf(' ]\n');

fprintf('Mean effectiveness      : %.6f\n', lambda_mean);

fprintf('Uniform comparison case : [');
fprintf(' %.6f', lambda_uni);
fprintf(' ]\n');


% ========================================================================
%  2. BUILD EXACT VEHICLE CONFIGURATION
%  ========================================================================

P = harness.config( ...
    package_mass = mass_pkg, ...
    package_r_c  = r_c_pkg);

UAV     = P.UAV;
UAV_pkg = P.UAV_pkg;

F2      = UAV.F2;
F2pinv  = UAV.F2_pinv;

g       = P.par.g;
kf      = UAV.k_f;

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' ALLOCATION MATRIX INFORMATION\n');
fprintf('=====================================================================\n');

fprintf('size(F2) = %d x %d\n', size(F2,1), size(F2,2));
fprintf('rank(F2) = %d\n', rank(F2));

disp('F2 =');
disp(F2);

fprintf('||F2*F2_pinv - I||_F = %.3e\n', ...
    norm(F2*F2pinv - eye(size(F2,1)), 'fro'));


% ========================================================================
%  3. EFFECTIVE WRENCH MAPPING
%
%      rotor command:
%          T_r,cmd = F2_pinv * w_cmd
%
%      degraded rotor thrust:
%          T_r,actual = diag(lambda_r) * T_r,cmd
%
%      actual wrench:
%          w_actual = F2 * diag(lambda_r) * F2_pinv * w_cmd
%
%      therefore
%
%          A_lambda = F2 * diag(lambda_r) * F2_pinv
%
%      and
%
%          w_actual = A_lambda * w_cmd
%
%  ========================================================================

A_asym = F2 * diag(lambda_asym) * F2pinv;
A_uni  = F2 * diag(lambda_uni)  * F2pinv;

A_expected_uni = lambda_mean * eye(4);

Delta_asym = A_asym - lambda_mean*eye(4);

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' EFFECTIVE COMMANDED-WRENCH -> DELIVERED-WRENCH MAP\n');
fprintf('=====================================================================\n');

disp('A_asym = F2 * diag(lambda_asym) * F2_pinv');
disp(A_asym);

disp('A_uniform = F2 * diag(lambda_mean) * F2_pinv');
disp(A_uni);

disp('Expected uniform mapping = lambda_mean * I');
disp(A_expected_uni);

fprintf('Uniform-map consistency error ||A_uni-lambda*I||_F = %.3e\n', ...
    norm(A_uni - A_expected_uni, 'fro'));

disp('Delta_asym = A_asym - lambda_mean*I');
disp(Delta_asym);


% ========================================================================
%  4. DIAGONAL EFFECTIVENESS AND OFF-DIAGONAL COUPLING
%  ========================================================================

diag_asym = diag(A_asym);

offdiag = A_asym - diag(diag_asym);

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' DIAGONAL GAINS AND CROSS-COUPLING\n');
fprintf('=====================================================================\n');

fprintf('Diagonal terms of A_asym:\n');
fprintf('  collective -> collective : %.6f\n', A_asym(1,1));
fprintf('  roll       -> roll       : %.6f\n', A_asym(2,2));
fprintf('  pitch      -> pitch      : %.6f\n', A_asym(3,3));
fprintf('  yaw        -> yaw        : %.6f\n', A_asym(4,4));

fprintf('\nOff-diagonal entries:\n');
disp(offdiag);

fprintf('Largest absolute off-diagonal entry = %.6e\n', ...
    max(abs(offdiag(:))));


% ========================================================================
%  5. PURE HOVER TEST
%
%  Project sign convention:
%
%      T_hover = -m*g
%
%  Wrench ordering assumed from F2:
%
%      w = [T; tau_x; tau_y; tau_z]
%
%  This is the most important physical test:
%
%      command pure collective
%
%      [T_hover; 0; 0; 0]
%
%  and inspect whether asymmetric rotor effectiveness creates moments.
%  ========================================================================

m_tot = UAV_pkg.m;

T_hover = -m_tot * g;

w_hover_cmd = [T_hover; 0; 0; 0];

w_hover_asym = A_asym * w_hover_cmd;
w_hover_uni  = A_uni  * w_hover_cmd;

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' PURE HOVER WRENCH TEST\n');
fprintf('=====================================================================\n');

fprintf('Total plant mass = %.6f kg\n', m_tot);
fprintf('Required hover collective = %.6f N\n', T_hover);

fprintf('\nCOMMAND:\n');
fprintf('  T     = %+12.6f N\n',  w_hover_cmd(1));
fprintf('  tau_x = %+12.6f Nm\n', w_hover_cmd(2));
fprintf('  tau_y = %+12.6f Nm\n', w_hover_cmd(3));
fprintf('  tau_z = %+12.6f Nm\n', w_hover_cmd(4));

fprintf('\nDELIVERED - ASYMMETRIC ROTORS:\n');
fprintf('  T     = %+12.6f N\n',  w_hover_asym(1));
fprintf('  tau_x = %+12.6f Nm\n', w_hover_asym(2));
fprintf('  tau_y = %+12.6f Nm\n', w_hover_asym(3));
fprintf('  tau_z = %+12.6f Nm\n', w_hover_asym(4));

fprintf('\nDELIVERED - UNIFORM ROTORS:\n');
fprintf('  T     = %+12.6f N\n',  w_hover_uni(1));
fprintf('  tau_x = %+12.6f Nm\n', w_hover_uni(2));
fprintf('  tau_y = %+12.6f Nm\n', w_hover_uni(3));
fprintf('  tau_z = %+12.6f Nm\n', w_hover_uni(4));

tau_par_asym = w_hover_asym(2:4);
tau_par_norm = norm(tau_par_asym);

fprintf('\nAsymmetric pure-hover parasitic moment:\n');
fprintf('  [tau_x tau_y tau_z] = [%+.6f  %+.6f  %+.6f] Nm\n', ...
    tau_par_asym(1), tau_par_asym(2), tau_par_asym(3));

fprintf('  ||tau_parasitic|| = %.6f Nm\n', tau_par_norm);

fprintf('\nRelative collective effectiveness in hover = %.4f\n', ...
    w_hover_asym(1)/w_hover_cmd(1));


% ========================================================================
%  6. DIRECT ROTOR-LEVEL HOVER CHECK
%
%  First compute what the NORMAL allocator commands for pure hover.
%  Then pass those rotor thrusts through the asymmetric efficiencies.
%  ========================================================================

Tr_hover_cmd = F2pinv * w_hover_cmd;

Tr_hover_actual_asym = diag(lambda_asym) * Tr_hover_cmd;
Tr_hover_actual_uni  = diag(lambda_uni)  * Tr_hover_cmd;

w_hover_asym_direct = F2 * Tr_hover_actual_asym;
w_hover_uni_direct  = F2 * Tr_hover_actual_uni;

Omega_hover_cmd = sqrt(abs(Tr_hover_cmd) / kf);

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' ROTOR-LEVEL HOVER MAPPING\n');
fprintf('=====================================================================\n');

fprintf('\nNominal allocator rotor thrust commands [N]:\n');
for i = 1:6
    fprintf('  rotor %d : %+10.6f N\n', i, Tr_hover_cmd(i));
end

fprintf('\nCorresponding rotor speeds [rad/s]:\n');
for i = 1:6
    fprintf('  rotor %d : %10.3f rad/s\n', i, Omega_hover_cmd(i));
end

fprintf('\nAfter ASYMMETRIC effectiveness, actual rotor thrusts [N]:\n');
for i = 1:6
    fprintf('  rotor %d : %+10.6f N   lambda = %.4f\n', ...
        i, Tr_hover_actual_asym(i), lambda_asym(i));
end

fprintf('\nCross-check direct rotor mapping vs A_asym mapping:\n');
fprintf('  ||w_direct - w_A|| = %.3e\n', ...
    norm(w_hover_asym_direct - w_hover_asym));

fprintf('Uniform cross-check:\n');
fprintf('  ||w_direct - w_A|| = %.3e\n', ...
    norm(w_hover_uni_direct - w_hover_uni));


% ========================================================================
%  7. FAULT-AWARE HOVER FEASIBILITY TEST
%
%  Ask a different question:
%
%    If the allocator KNEW the true six rotor efficiencies, could it find
%    rotor commands that produce the required pure-hover wrench?
%
%  We solve
%
%      F2 * diag(lambda) * T_r_cmd = w_hover
%
%  using the pseudoinverse.
%
%  If this solution is within the rotor-speed limits, then the physical plant
%  still has enough authority for hover. A later closed-loop failure therefore
%  cannot be explained simply by "the degraded plant cannot fly".
%  ========================================================================

G_asym = F2 * diag(lambda_asym);

Tr_faultaware = pinv(G_asym) * w_hover_cmd;

w_faultaware = G_asym * Tr_faultaware;

Omega_faultaware = sqrt(abs(Tr_faultaware) / kf);

same_thrust_direction = all(sign(T_hover) * Tr_faultaware >= -1e-10);

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' FAULT-AWARE HOVER FEASIBILITY TEST\n');
fprintf('=====================================================================\n');

fprintf('Residual ||G*T_r - w_hover|| = %.3e\n', ...
    norm(w_faultaware - w_hover_cmd));

fprintf('\nRequired rotor thrust commands if lambda were known:\n');
for i = 1:6
    fprintf('  rotor %d : %+10.6f N   -> %8.2f rad/s\n', ...
        i, Tr_faultaware(i), Omega_faultaware(i));
end

fprintf('\nRotor speed limits:\n');
fprintf('  Omega_min = %.2f rad/s\n', UAV.Omega_min);
fprintf('  Omega_max = %.2f rad/s\n', UAV.Omega_max);

fprintf('\nFault-aware solution:\n');
fprintf('  min Omega = %.2f rad/s\n', min(Omega_faultaware));
fprintf('  max Omega = %.2f rad/s\n', max(Omega_faultaware));
fprintf('  same thrust direction as normal hover = %d\n', same_thrust_direction);

if max(Omega_faultaware) <= UAV.Omega_max && same_thrust_direction
    fprintf('  VERDICT: static hover is feasible under this rotor degradation.\n');
else
    fprintf('  VERDICT: static hover is NOT feasible within the simple rotor limits.\n');
end


% ========================================================================
%  8. WRENCH REQUIRED IF WE INVERT THE 4x4 EFFECTIVE MAP
%
%  This is another way to see the coupled compensation required:
%
%      A_asym * w_cmd_required = w_hover_desired
%
%  If A were simply lambda*I, only collective would increase.
%  If asymmetric, corrective moments will generally also be required.
%  ========================================================================

if rcond(A_asym) > 1e-10
    w_cmd_comp = A_asym \ w_hover_cmd;

    fprintf('\n');
    fprintf('=====================================================================\n');
    fprintf(' COUPLED WRENCH COMPENSATION REQUIRED FOR PERFECT HOVER\n');
    fprintf('=====================================================================\n');

    fprintf('Required commanded wrench under A_asym:\n');
    fprintf('  T_cmd     = %+12.6f N\n',  w_cmd_comp(1));
    fprintf('  tau_x_cmd = %+12.6f Nm\n', w_cmd_comp(2));
    fprintf('  tau_y_cmd = %+12.6f Nm\n', w_cmd_comp(3));
    fprintf('  tau_z_cmd = %+12.6f Nm\n', w_cmd_comp(4));

    fprintf('\nCompare with simple scalar compensation:\n');
    fprintf('  T_hover/lambda_mean = %+12.6f N\n', ...
        T_hover/lambda_mean);

    fprintf('\nCheck A_asym*w_cmd_comp:\n');
    disp(A_asym*w_cmd_comp);
else
    warning('A_asym is close to singular; inverse-wrench test skipped.');
end


% ========================================================================
%  9. DIMENSIONLESS COUPLING VIEW
%
%  Matrix entries mix force and moment units, so comparing raw off-diagonal
%  numbers directly is not always meaningful.
%
%  Construct characteristic wrench scales based on hover thrust and the
%  geometry encoded in F2.
%  ========================================================================

arm_x = max(abs(F2(2,:)));
arm_y = max(abs(F2(3,:)));
yaw_scale_coeff = max(abs(F2(4,:)));

w_scale = [ ...
    abs(T_hover);
    abs(T_hover)*max(arm_x, eps);
    abs(T_hover)*max(arm_y, eps);
    abs(T_hover)*max(yaw_scale_coeff, eps)];

S = diag(w_scale);

A_asym_norm = S \ A_asym * S;
A_uni_norm  = S \ A_uni  * S;

Delta_norm = A_asym_norm - lambda_mean*eye(4);

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' DIMENSIONLESS NORMALIZED WRENCH MAPPING\n');
fprintf('=====================================================================\n');

fprintf('Characteristic scales [T, tau_x, tau_y, tau_z]:\n');
fprintf('  [%g  %g  %g  %g]\n', w_scale);

disp('Normalized A_asym:');
disp(A_asym_norm);

disp('Normalized Delta = A_asym_norm - lambda_mean*I:');
disp(Delta_norm);

off_norm = Delta_norm - diag(diag(Delta_norm));

fprintf('Largest normalized off-diagonal coupling = %.4f\n', ...
    max(abs(off_norm(:))));

fprintf('Normalized off-diagonal Frobenius norm = %.4f\n', ...
    norm(off_norm, 'fro'));


% ========================================================================
%  10. OPTIONAL NONLINEAR SIMULATION REPLAY
%  ========================================================================

if RUN_SIMULATIONS

    fprintf('\n');
    fprintf('=====================================================================\n');
    fprintf(' RUNNING NONLINEAR ASYMMETRIC / UNIFORM COMPARISON\n');
    fprintf('=====================================================================\n');

    lam_asym_sig = disturb.fault( ...
        MC.ref.time, ...
        t0 = 0, ...
        lambda = lambda_asym);

    lam_uni_sig = disturb.fault( ...
        MC.ref.time, ...
        t0 = 0, ...
        lambda = lambda_uni);

    R_asym = harness.scenario( ...
        MC.ref, ...
        config  = P, ...
        package = 1, ...
        mrac    = 1, ...
        robust  = 1, ...
        lambda  = lam_asym_sig, ...
        model   = MODEL, ...
        label   = 'Asymmetric lambda_r');

    R_uni = harness.scenario( ...
        MC.ref, ...
        config  = P, ...
        package = 1, ...
        mrac    = 1, ...
        robust  = 1, ...
        lambda  = lam_uni_sig, ...
        model   = MODEL, ...
        label   = 'Uniform mean lambda_r');

    runs = {R_asym, R_uni};
    tags = {'Asymmetric \lambda_r', 'Uniform mean \lambda_r'};


    % ---- Tracking report ---------------------------------------------

    harness.report( ...
        runs, tags, ...
        title = sprintf('MC run %d - asymmetric vs uniform effectiveness', k));


    % ---- Actuator usage ----------------------------------------------

    viz.actuators( ...
        runs, tags, UAV, ...
        name = sprintf('MC run %d - actuator comparison', k));


    % ---- Attitude ----------------------------------------------------

    viz.attitude( ...
        runs, tags, ...
        name = sprintf('MC run %d - attitude comparison', k));


    % ---- Tracking errors ---------------------------------------------

    viz.errors( ...
        runs, tags, ...
        name = sprintf('MC run %d - tracking comparison', k));


    % ---- Adaptive effectiveness -------------------------------------

    viz.effectiveness( ...
        {R_asym}, {'Asymmetric'}, ...
        lambda = lam_asym_sig, ...
        t_fault = 0, ...
        name = sprintf('MC run %d - asymmetric adaptive effectiveness', k));

    viz.effectiveness( ...
        {R_uni}, {'Uniform'}, ...
        lambda = lam_uni_sig, ...
        t_fault = 0, ...
        name = sprintf('MC run %d - uniform adaptive effectiveness', k));


    % ====================================================================
    %  11. PRINT LAMBDA-HAT SUMMARY
    %  ====================================================================

    fprintf('\n');
    fprintf('=====================================================================\n');
    fprintf(' ADAPTIVE LAMBDA-HAT SUMMARY\n');
    fprintf('=====================================================================\n');

    printLambdaSummary(R_asym, 'ASYMMETRIC');
    printLambdaSummary(R_uni,  'UNIFORM');


    % ====================================================================
    %  12. PRINT OPTIONAL LOGGED WRENCH / THRUST INFORMATION
    %  ====================================================================

    fprintf('\n');
    fprintf('=====================================================================\n');
    fprintf(' AVAILABLE LOGGED CONTROL SIGNALS\n');
    fprintf('=====================================================================\n');

    inspectSignal(R_asym, 'f_c',   'asymmetric f_c');
    inspectSignal(R_asym, 'T_c',   'asymmetric T_c');
    inspectSignal(R_asym, 'T_out', 'asymmetric T_out');
    inspectSignal(R_asym, 'w_a',   'asymmetric w_a');

    inspectSignal(R_uni, 'f_c',   'uniform f_c');
    inspectSignal(R_uni, 'T_c',   'uniform T_c');
    inspectSignal(R_uni, 'T_out', 'uniform T_out');
    inspectSignal(R_uni, 'w_a',   'uniform w_a');

end


% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function printLambdaSummary(R, label)

fprintf('\n%s\n', label);

if ~isfield(R, 'lambda_hat')
    fprintf('  lambda_hat not logged.\n');
    return
end

L = R.lambda_hat;

for j = 1:size(L,2)
    fprintf('  lambda_hat_%d: min=%8.4f  max=%8.4f  final=%8.4f\n', ...
        j, min(L(:,j)), max(L(:,j)), L(end,j));
end

end


function inspectSignal(R, fieldname, label)

fprintf('\n%s\n', label);

if ~isfield(R, fieldname)
    fprintf('  not logged\n');
    return
end

X = R.(fieldname);

fprintf('  size = %d x %d\n', size(X,1), size(X,2));

for j = 1:size(X,2)
    x = X(:,j);

    fprintf('  component %d: min=%+.6g  max=%+.6g  RMS=%.6g\n', ...
        j, min(x), max(x), sqrt(mean(x.^2)));
end

end
%% FANCY PORCO DIO
% ========================================================================
%  MC RUN 9 - CORRECTED FAILURE CHRONOLOGY
%
%  Purpose
%  -------
%  Reproduce the exact catastrophic asymmetric Monte Carlo run and determine
%  the temporal ordering of:
%
%    1) Lambda_hat departure
%    2) planner collective demand
%    3) requested rotor saturation
%    4) actual actuator saturation
%    5) attitude departure
%    6) body-rate growth
%    7) position-error divergence
%
%  Nothing in the controller is changed.
%
%  IMPORTANT:
%  - Euler angles logged by harness.unpack are ZYX:
%
%        R.euler = [yaw pitch roll]
%
%  - Rotor thrusts use the signed convention:
%
%        Omega = sqrt(abs(T_r)/k_f)
%
%  ========================================================================

clearvars;
close all;
clc;


%% ========================================================================
%  USER SETTINGS
%  ========================================================================

MC_FILE = 'results/montecarlo_20260919_163739.mat';
MODEL   = 'Final';

k = 9;

% Robust PBMRAC Lambda projection
LAMBDA_INNER_LOW = 0.45;
LAMBDA_OUTER_LOW = 0.35;

% Attitude planner collective cap
T_CAP = 22.9308;              % [N]

% Event thresholds
EP1 = 0.10;                   % [m]
EP2 = 0.50;                   % [m]
EP3 = 1.00;                   % [m]
EP4 = 5.00;                   % [m]

ATT1 = 5;                     % [deg]
ATT2 = 10;                    % [deg]
ATT3 = 30;                    % [deg]

OMEGA1 = 1;                   % [rad/s]
OMEGA2 = 3;                   % [rad/s]

RUN_STANDARD_VIZ = true;


%% ========================================================================
%  1. LOAD EXACT MONTE CARLO SAMPLE
%  ========================================================================

load(MC_FILE,'MC');

mass_pkg    = MC.samples.mass(k);
r_c_pkg     = MC.samples.r_c(k,:).';
lambda_asym = MC.samples.lambda(k,:).';

fprintf('\n');
fprintf('============================================================\n');
fprintf(' MC RUN %d - CORRECTED FAILURE CHRONOLOGY\n',k);
fprintf('============================================================\n');

fprintf('Package mass : %.6f kg\n',mass_pkg);
fprintf('Package CoM  : [%+.6f %+.6f %+.6f] m\n', ...
    r_c_pkg(1),r_c_pkg(2),r_c_pkg(3));

fprintf('lambda       : [');
fprintf(' %.6f',lambda_asym);
fprintf(' ]\n');


%% ========================================================================
%  2. BUILD EXACT PLANT CONFIGURATION
%  ========================================================================

P = harness.config( ...
    package_mass = mass_pkg, ...
    package_r_c  = r_c_pkg);

UAV = P.UAV;


%% ========================================================================
%  3. REPLAY EXACT ASYMMETRIC ROBUST PBMRAC RUN
%  ========================================================================

lam_sig = disturb.fault( ...
    MC.ref.time, ...
    t0 = 0, ...
    lambda = lambda_asym);

R = harness.scenario( ...
    MC.ref, ...
    config  = P, ...
    package = 1, ...
    mrac    = 1, ...
    robust  = 1, ...
    lambda  = lam_sig, ...
    model   = MODEL, ...
    label   = 'MC run 9 asymmetric');


%% ========================================================================
%  4. CHECK REQUIRED LOGGED SIGNALS
%  ========================================================================

required = { ...
    'p','p_d','v','v_d', ...
    'lambda_hat','f_c','T_c', ...
    'T_out','w_a','euler','omegab'};

for ii = 1:numel(required)

    if ~isfield(R,required{ii})
        error('Required logged signal "%s" is missing.',required{ii});
    end

end


%% ========================================================================
%  5. POSITION / VELOCITY TRACKING ERRORS
%  ========================================================================

t_ep = R.t;

ep = R.p - R.p_d;
ev = R.v - R.v_d;

ep_norm = vecnorm(ep,2,2);
ev_norm = vecnorm(ev,2,2);


%% ========================================================================
%  6. ADAPTIVE EFFECTIVENESS STATE
%  ========================================================================

t_lambda = R.lambda_hat_t;
Lambda   = R.lambda_hat;

lambda_min = min(Lambda,[],2);
lambda_max = max(Lambda,[],2);


%% ========================================================================
%  7. OUTER-LOOP FORCE COMMAND
%
%  harness.unpack defines f_c as the commanded inertial force from the
%  outer loop.
%  ========================================================================

t_fc = R.f_c_t;
fc   = R.f_c;

fc_norm = vecnorm(fc,2,2);


%% ========================================================================
%  8. PLANNER COLLECTIVE COMMAND
%  ========================================================================

t_Tc = R.T_c_t;
Tc   = R.T_c(:,1);


%% ========================================================================
%  9. EULER ATTITUDE
%
%  harness.unpack stores ZYX Euler angles:
%
%       [yaw pitch roll]
%
%  Therefore:
%
%       yaw   = column 1
%       pitch = column 2
%       roll  = column 3
%  ========================================================================

t_att = R.euler_t;
Edeg  = rad2deg(R.euler);

yaw_deg   = Edeg(:,1);
pitch_deg = Edeg(:,2);
roll_deg  = Edeg(:,3);

tilt_deg = max(abs([roll_deg pitch_deg]),[],2);

fprintf('\nInitial Euler [yaw pitch roll] = ');
fprintf('[%+.3f %+.3f %+.3f] deg\n', ...
    yaw_deg(1),pitch_deg(1),roll_deg(1));


%% ========================================================================
%  10. BODY ANGULAR VELOCITY
%  ========================================================================

t_omega = R.omegab_t;

omega_body = R.omegab;
omega_norm = vecnorm(omega_body,2,2);


%% ========================================================================
%  11. PREDICTOR ERROR, IF AVAILABLE
%  ========================================================================

has_pred = isfield(R,'e_hat_PRED') && ...
           isfield(R,'e_hat_PRED_t') && ...
           size(R.e_hat_PRED,2) >= 6;

if has_pred

    t_pred = R.e_hat_PRED_t;

    ep_pred = R.e_hat_PRED(:,1:3);
    ev_pred = R.e_hat_PRED(:,4:6);

    ep_pred_norm = vecnorm(ep_pred,2,2);
    ev_pred_norm = vecnorm(ev_pred,2,2);

else

    fprintf('\nNOTE: e_hat_PRED is not logged; ');
    fprintf('predictor-error timing unavailable.\n');

    t_pred       = [];
    ep_pred_norm = [];
    ev_pred_norm = [];

end


%% ========================================================================
%  12. RECONSTRUCT REQUESTED ROTOR SPEEDS BEFORE SATURATION
%
%  R.w_a is the six-component control wrench:
%
%       [Fx Fy Fz Mx My Mz]
%
%  The allocation coordinates are:
%
%       u_alloc = [T Mx My Mz]
%
%  In this project's sign convention positive body-z force corresponds to
%  negative collective T:
%
%       T = -Fz
%
%  Hence:
%
%       u_alloc = [-w_a,z  Mx My Mz]
%
%  Nominal allocator:
%
%       T_r,cmd = F2_pinv * u_alloc
%
%  Rotor speed:
%
%       Omega_cmd = sqrt(abs(T_r,cmd)/k_f)
%
%  This is the UNSATURATED rotor demand.
%  ========================================================================

t_wa = R.w_a_t;
wa   = R.w_a;

u_alloc_req = [ ...
    -wa(:,3), ...
     wa(:,4), ...
     wa(:,5), ...
     wa(:,6)];

Tr_req = (UAV.F2_pinv * u_alloc_req.').';

Omega_req = sqrt(abs(Tr_req) / UAV.k_f);

Omega_req_max = max(Omega_req,[],2);
Omega_req_min = min(Omega_req,[],2);

upper_req_sat = any(Omega_req >= UAV.Omega_max,2);
lower_req_sat = any(Omega_req <= UAV.Omega_min,2);
any_req_sat   = upper_req_sat | lower_req_sat;


%% ========================================================================
%  13. ACTUAL ACTUATOR LIMIT FROM T_out
%
%  T_out is the thrust actually produced after the rotor effectiveness.
%
%  Since:
%
%       T_out_i = lambda_i * T_after_limit_i
%
%  recover the equivalent pre-effectiveness, post-saturation thrust:
%
%       T_after_limit_i = T_out_i/lambda_i
%
%  and therefore
%
%       Omega_actual = sqrt(abs(T_after_limit_i)/k_f)
%
%  This quantity can NEVER exceed Omega_max because it is downstream of the
%  physical saturation.
%  ========================================================================

t_out = R.T_out_t;
Tout  = R.T_out;

lambda_row = reshape(lambda_asym,1,[]);

Tr_after_limit = Tout ./ lambda_row;

Omega_actual = sqrt(abs(Tr_after_limit) / UAV.k_f);

Omega_actual_max = max(Omega_actual,[],2);
Omega_actual_min = min(Omega_actual,[],2);

% Small numerical tolerance because logged thrust may differ by roundoff.
upper_actual_sat = any( ...
    Omega_actual >= 0.999*UAV.Omega_max, ...
    2);

lower_actual_sat = any( ...
    Omega_actual <= 1.001*UAV.Omega_min, ...
    2);

any_actual_sat = upper_actual_sat | lower_actual_sat;


%% ========================================================================
%  14. PRINT ACTUATOR RECONSTRUCTION CHECK
%  ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' ACTUATOR RECONSTRUCTION CHECK\n');
fprintf('============================================================\n');

fprintf('Requested rotor speed:\n');
fprintf('  min over run = %.3f rad/s\n',min(Omega_req(:)));
fprintf('  max over run = %.3f rad/s\n',max(Omega_req(:)));

fprintf('\nVehicle limits:\n');
fprintf('  Omega_min    = %.3f rad/s\n',UAV.Omega_min);
fprintf('  Omega_max    = %.3f rad/s\n',UAV.Omega_max);

fprintf('\nRequested-command saturation fraction:\n');
fprintf('  upper        = %.2f %%\n',100*mean(upper_req_sat));
fprintf('  lower        = %.2f %%\n',100*mean(lower_req_sat));
fprintf('  either       = %.2f %%\n',100*mean(any_req_sat));

fprintf('\nPost-limit / actual rotor speed:\n');
fprintf('  min over run = %.3f rad/s\n',min(Omega_actual(:)));
fprintf('  max over run = %.3f rad/s\n',max(Omega_actual(:)));

fprintf('\nActual-limit occupancy:\n');
fprintf('  upper        = %.2f %%\n',100*mean(upper_actual_sat));
fprintf('  lower        = %.2f %%\n',100*mean(lower_actual_sat));
fprintf('  either       = %.2f %%\n',100*mean(any_actual_sat));


%% ========================================================================
%  15. EVENT TABLE
%
%  IMPORTANT:
%     names = string array
%     times = numeric double array
%
%  This fixes the previous "Conversion to cell from double" error.
%  ========================================================================

names = strings(0,1);
times = zeros(0,1);


% ------------------------------------------------------------------------
% Predictor
% -------------------------------------------------------------------------

if has_pred

    [names,times] = addEvt( ...
        names,times, ...
        "predictor |e_v| > 0.05 m/s", ...
        firstAbove(t_pred,ev_pred_norm,0.05));

    [names,times] = addEvt( ...
        names,times, ...
        "predictor |e_v| > 0.10 m/s", ...
        firstAbove(t_pred,ev_pred_norm,0.10));

    [names,times] = addEvt( ...
        names,times, ...
        "predictor |e_v| > 0.50 m/s", ...
        firstAbove(t_pred,ev_pred_norm,0.50));

end


% ------------------------------------------------------------------------
% Lambda_hat
% -------------------------------------------------------------------------

[names,times] = addEvt( ...
    names,times, ...
    "min Lambda_hat < 0.90", ...
    firstBelow(t_lambda,lambda_min,0.90));

[names,times] = addEvt( ...
    names,times, ...
    "min Lambda_hat < 0.80", ...
    firstBelow(t_lambda,lambda_min,0.80));

[names,times] = addEvt( ...
    names,times, ...
    "min Lambda_hat < 0.70", ...
    firstBelow(t_lambda,lambda_min,0.70));

[names,times] = addEvt( ...
    names,times, ...
    "Lambda enters projection region (<0.45)", ...
    firstBelow(t_lambda,lambda_min,LAMBDA_INNER_LOW));

[names,times] = addEvt( ...
    names,times, ...
    "Lambda near lower wall (<0.36)", ...
    firstBelow(t_lambda,lambda_min,0.36));


% ------------------------------------------------------------------------
% Planner collective
% -------------------------------------------------------------------------

[names,times] = addEvt( ...
    names,times, ...
    "planner T_c > 90% cap", ...
    firstAbove(t_Tc,Tc,0.90*T_CAP));

[names,times] = addEvt( ...
    names,times, ...
    "planner T_c > 99% cap", ...
    firstAbove(t_Tc,Tc,0.99*T_CAP));


% ------------------------------------------------------------------------
% Rotor demand / actuator saturation
% -------------------------------------------------------------------------

[names,times] = addEvt( ...
    names,times, ...
    "requested rotor demand reaches Omega_max", ...
    firstTrue(t_wa,upper_req_sat));

[names,times] = addEvt( ...
    names,times, ...
    "requested rotor demand reaches Omega_min", ...
    firstTrue(t_wa,lower_req_sat));

[names,times] = addEvt( ...
    names,times, ...
    "actual rotor sits at upper limit", ...
    firstTrue(t_out,upper_actual_sat));

[names,times] = addEvt( ...
    names,times, ...
    "actual rotor sits at lower limit", ...
    firstTrue(t_out,lower_actual_sat));


% ------------------------------------------------------------------------
% Attitude
% -------------------------------------------------------------------------

[names,times] = addEvt( ...
    names,times, ...
    "roll/pitch excursion > 5 deg", ...
    firstAbove(t_att,tilt_deg,ATT1));

[names,times] = addEvt( ...
    names,times, ...
    "roll/pitch excursion > 10 deg", ...
    firstAbove(t_att,tilt_deg,ATT2));

[names,times] = addEvt( ...
    names,times, ...
    "roll/pitch excursion > 30 deg", ...
    firstAbove(t_att,tilt_deg,ATT3));


% ------------------------------------------------------------------------
% Angular velocity
% -------------------------------------------------------------------------

[names,times] = addEvt( ...
    names,times, ...
    "|omega| > 1 rad/s", ...
    firstAbove(t_omega,omega_norm,OMEGA1));

[names,times] = addEvt( ...
    names,times, ...
    "|omega| > 3 rad/s", ...
    firstAbove(t_omega,omega_norm,OMEGA2));


% ------------------------------------------------------------------------
% Position tracking
% -------------------------------------------------------------------------

[names,times] = addEvt( ...
    names,times, ...
    "|e_p| > 0.10 m", ...
    firstAbove(t_ep,ep_norm,EP1));

[names,times] = addEvt( ...
    names,times, ...
    "|e_p| > 0.50 m", ...
    firstAbove(t_ep,ep_norm,EP2));

[names,times] = addEvt( ...
    names,times, ...
    "|e_p| > 1.00 m", ...
    firstAbove(t_ep,ep_norm,EP3));

[names,times] = addEvt( ...
    names,times, ...
    "|e_p| > 5.00 m", ...
    firstAbove(t_ep,ep_norm,EP4));


% ------------------------------------------------------------------------
% Outer-loop force
% -------------------------------------------------------------------------

[names,times] = addEvt( ...
    names,times, ...
    "||f_c|| > planner thrust cap", ...
    firstAbove(t_fc,fc_norm,T_CAP));

[names,times] = addEvt( ...
    names,times, ...
    "||f_c|| > 50 N", ...
    firstAbove(t_fc,fc_norm,50));

[names,times] = addEvt( ...
    names,times, ...
    "||f_c|| > 100 N", ...
    firstAbove(t_fc,fc_norm,100));


%% ========================================================================
%  16. SORT AND PRINT CHRONOLOGY
%  ========================================================================

Evt = table( ...
    names, ...
    times, ...
    'VariableNames',{'Event','FirstTime_s'});

Evt = sortrows(Evt,'FirstTime_s');

fprintf('\n');
fprintf('============================================================\n');
fprintf(' CORRECTED CHRONOLOGICAL ORDER\n');
fprintf('============================================================\n\n');

disp(Evt);


%% ========================================================================
%  17. PRINT HEADLINE TIMES
%  ========================================================================

t_L90   = firstBelow(t_lambda,lambda_min,0.90);
t_L70   = firstBelow(t_lambda,lambda_min,0.70);
t_L45   = firstBelow(t_lambda,lambda_min,0.45);
t_L36   = firstBelow(t_lambda,lambda_min,0.36);

t_T90   = firstAbove(t_Tc,Tc,0.90*T_CAP);
t_T99   = firstAbove(t_Tc,Tc,0.99*T_CAP);

t_reqUp = firstTrue(t_wa,upper_req_sat);
t_actUp = firstTrue(t_out,upper_actual_sat);

t_att5  = firstAbove(t_att,tilt_deg,5);
t_att10 = firstAbove(t_att,tilt_deg,10);
t_att30 = firstAbove(t_att,tilt_deg,30);

t_w1    = firstAbove(t_omega,omega_norm,1);
t_w3    = firstAbove(t_omega,omega_norm,3);

t_ep01  = firstAbove(t_ep,ep_norm,0.10);
t_ep05  = firstAbove(t_ep,ep_norm,0.50);
t_ep1   = firstAbove(t_ep,ep_norm,1.00);

t_fcCap = firstAbove(t_fc,fc_norm,T_CAP);


fprintf('\n');
fprintf('============================================================\n');
fprintf(' HEADLINE FAILURE TIMES\n');
fprintf('============================================================\n');

printTime('Lambda_hat < 0.90',t_L90);
printTime('Lambda_hat < 0.70',t_L70);
printTime('Lambda enters projection region (<0.45)',t_L45);
printTime('Lambda approaches lower wall (<0.36)',t_L36);

fprintf('\n');

printTime('Planner T_c > 90%% cap',t_T90);
printTime('Planner T_c > 99%% cap',t_T99);

fprintf('\n');

printTime('Requested rotor reaches Omega_max',t_reqUp);
printTime('Actual rotor reaches upper stop',t_actUp);

fprintf('\n');

printTime('Roll/pitch > 5 deg',t_att5);
printTime('Roll/pitch > 10 deg',t_att10);
printTime('Roll/pitch > 30 deg',t_att30);

fprintf('\n');

printTime('|omega| > 1 rad/s',t_w1);
printTime('|omega| > 3 rad/s',t_w3);

fprintf('\n');

printTime('|e_p| > 0.10 m',t_ep01);
printTime('|e_p| > 0.50 m',t_ep05);
printTime('|e_p| > 1.00 m',t_ep1);

fprintf('\n');

printTime('||f_c|| > planner cap',t_fcCap);


%% ========================================================================
%  18. AUTOMATIC FAILURE-ONSET WINDOW
%  ========================================================================

headline_times = [ ...
    t_L90, ...
    t_T90, ...
    t_reqUp, ...
    t_actUp, ...
    t_att5, ...
    t_w1, ...
    t_ep01, ...
    t_L45, ...
    t_ep1];

headline_times = headline_times(isfinite(headline_times));

if isempty(headline_times)

    t_zoom_0 = 0;
    t_zoom_f = min(20,R.t(end));

else

    t_zoom_0 = max(0,min(headline_times)-1);
    t_zoom_f = min(R.t(end),max(headline_times)+5);

end

fprintf('\nSuggested detailed zoom window: %.3f -> %.3f s\n', ...
    t_zoom_0,t_zoom_f);


%% ========================================================================
%  19. FULL CHRONOLOGY FIGURE
%  ========================================================================

figure( ...
    'Name','MC run 9 - corrected failure chronology', ...
    'Color','w', ...
    'Position',[80 30 1200 980]);

tl = tiledlayout(6,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


% ------------------------------------------------------------------------
% Lambda
% -------------------------------------------------------------------------

nexttile;

plot(t_lambda,lambda_min,'LineWidth',1.4);
hold on;

yline(LAMBDA_INNER_LOW,'--','inner lower');
yline(LAMBDA_OUTER_LOW,'--','outer lower');

grid on;

ylabel('min(lambda hat)');

title('Adaptive effectiveness');


% ------------------------------------------------------------------------
% Planner collective
% -------------------------------------------------------------------------

nexttile;

plot(t_Tc,Tc,'LineWidth',1.4);
hold on;

yline(T_CAP,'--','planner cap');

grid on;

ylabel('T_c [N]');

title('Planner collective demand');


% ------------------------------------------------------------------------
% Rotor requested vs actual
% -------------------------------------------------------------------------

nexttile;

plot(t_wa,Omega_req_max,'LineWidth',1.4);
hold on;

plot(t_out,Omega_actual_max,'LineWidth',1.2);

yline(UAV.Omega_max,'--','Omega max');

grid on;

ylabel('Omega [rad/s]');

legend( ...
    'requested max rotor', ...
    'actual max rotor', ...
    'Omega max', ...
    'Location','best');

title('Rotor demand and physical saturation');


% ------------------------------------------------------------------------
% Attitude
% -------------------------------------------------------------------------

nexttile;

plot(t_att,abs(roll_deg),'LineWidth',1.2);
hold on;

plot(t_att,abs(pitch_deg),'LineWidth',1.2);

yline(5,'--');
yline(10,'--');
yline(30,'--');

grid on;

ylabel('angle [deg]');

legend( ...
    '|roll|', ...
    '|pitch|', ...
    'Location','best');

title('Attitude departure');


% ------------------------------------------------------------------------
% Angular velocity
% -------------------------------------------------------------------------

nexttile;

plot(t_omega,omega_norm,'LineWidth',1.4);
hold on;

yline(1,'--');
yline(3,'--');

grid on;

ylabel('|omega| [rad/s]');

title('Body angular-rate growth');


% ------------------------------------------------------------------------
% Tracking
% -------------------------------------------------------------------------

nexttile;

semilogy(t_ep,ep_norm + eps,'LineWidth',1.4);
hold on;

yline(0.10,'--');
yline(0.50,'--');
yline(1.00,'--');

grid on;

ylabel('|e_p| [m]');
xlabel('time [s]');

title('Position tracking divergence');


title(tl,'MC run 9 - corrected temporal failure sequence');


%% ========================================================================
%  20. FAILURE ONSET ZOOM
%  ========================================================================

figure( ...
    'Name','MC run 9 - corrected failure onset zoom', ...
    'Color','w', ...
    'Position',[100 40 1200 950]);

tl2 = tiledlayout(6,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


nexttile;

plot(t_lambda,lambda_min,'LineWidth',1.5);
hold on;
yline(0.90,'--');
yline(0.70,'--');
yline(0.45,'--');
yline(0.35,'--');

xlim([t_zoom_0 t_zoom_f]);
grid on;

ylabel('min lambda');

title('Lambda hat');


nexttile;

plot(t_Tc,Tc,'LineWidth',1.5);
hold on;
yline(T_CAP,'--');

xlim([t_zoom_0 t_zoom_f]);
grid on;

ylabel('T_c [N]');

title('Collective demand');


nexttile;

plot(t_wa,Omega_req_max,'LineWidth',1.5);
hold on;
plot(t_out,Omega_actual_max,'LineWidth',1.2);

yline(UAV.Omega_max,'--');

xlim([t_zoom_0 t_zoom_f]);
grid on;

ylabel('Omega [rad/s]');

legend('requested','actual','limit','Location','best');

title('Rotor speed demand');


nexttile;

plot(t_att,abs(roll_deg),'LineWidth',1.4);
hold on;
plot(t_att,abs(pitch_deg),'LineWidth',1.4);

yline(5,'--');
yline(10,'--');
yline(30,'--');

xlim([t_zoom_0 t_zoom_f]);
grid on;

ylabel('angle [deg]');

legend('|roll|','|pitch|','Location','best');

title('Attitude');


nexttile;

plot(t_omega,omega_norm,'LineWidth',1.5);
hold on;
yline(1,'--');
yline(3,'--');

xlim([t_zoom_0 t_zoom_f]);
grid on;

ylabel('|omega|');

title('Angular velocity');


nexttile;

semilogy(t_ep,ep_norm + eps,'LineWidth',1.5);
hold on;
yline(0.10,'--');
yline(0.50,'--');
yline(1.00,'--');

xlim([t_zoom_0 t_zoom_f]);
grid on;

ylabel('|e_p| [m]');
xlabel('time [s]');

title('Tracking error');


title(tl2,'MC run 9 - failure onset');


%% ========================================================================
%  21. FORCE COMMAND FIGURE
%  ========================================================================

figure( ...
    'Name','MC run 9 - force growth', ...
    'Color','w', ...
    'Position',[150 120 950 600]);

plot(t_fc,fc_norm,'LineWidth',1.5);
hold on;

yline(T_CAP,'--','planner thrust cap');
yline(50,'--','50 N');
yline(100,'--','100 N');

grid on;

xlabel('time [s]');
ylabel('norm(f_c) [N]');

title('Outer-loop force command growth');


%% ========================================================================
%  22. OPTIONAL PREDICTOR-ERROR FIGURE
%  ========================================================================

if has_pred

    figure( ...
        'Name','MC run 9 - predictor error', ...
        'Color','w', ...
        'Position',[170 140 950 600]);

    semilogy(t_pred,ev_pred_norm + eps,'LineWidth',1.5);

    grid on;

    xlabel('time [s]');
    ylabel('predictor velocity error [m/s]');

    title('PBMRAC predictor velocity error');

end


%% ========================================================================
%  23. STANDARD PROJECT DIAGNOSTICS
%  ========================================================================

if RUN_STANDARD_VIZ

    viz.actuators( ...
        {R}, ...
        {'MC run 9'}, ...
        UAV, ...
        name='MC run 9 - actuator cross-check');

    viz.attitude( ...
        {R}, ...
        {'MC run 9'}, ...
        name='MC run 9 - attitude cross-check');

    viz.effectiveness( ...
        {R}, ...
        {'MC run 9'}, ...
        lambda=lam_sig, ...
        t_fault=0, ...
        name='MC run 9 - effectiveness cross-check');

    harness.report( ...
        {R}, ...
        {'MC run 9'}, ...
        title='MC run 9 corrected chronology');

end


%% ========================================================================
%  24. FINAL INFORMATION TO COPY BACK
%  ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' COPY THESE RESULTS BACK\n');
fprintf('============================================================\n');

fprintf(['1. CORRECTED CHRONOLOGICAL ORDER table\n' ...
         '2. HEADLINE FAILURE TIMES\n' ...
         '3. ACTUATOR RECONSTRUCTION CHECK\n' ...
         '4. In particular, compare reconstructed requested max Omega\n' ...
         '   with the existing viz.actuators value (~1489 rad/s).\n']);


% ========================================================================
%  LOCAL FUNCTIONS
%
%  If this section is pasted inside a larger MATLAB script, keep ALL local
%  functions at the very END of that .m file.
%  ========================================================================

function t0 = firstAbove(t,x,threshold)

    t = t(:);
    x = x(:);

    idx = find(x >= threshold,1,'first');

    if isempty(idx)
        t0 = NaN;
    else
        t0 = t(idx);
    end

end


function t0 = firstBelow(t,x,threshold)

    t = t(:);
    x = x(:);

    idx = find(x <= threshold,1,'first');

    if isempty(idx)
        t0 = NaN;
    else
        t0 = t(idx);
    end

end


function t0 = firstTrue(t,x)

    t = t(:);
    x = logical(x(:));

    idx = find(x,1,'first');

    if isempty(idx)
        t0 = NaN;
    else
        t0 = t(idx);
    end

end


function [names,times] = addEvt(names,times,name,t)

    if isfinite(t)

        names(end+1,1) = string(name);
        times(end+1,1) = double(t);

    end

end


function printTime(name,t)

    if isfinite(t)
        fprintf('%-46s : %9.4f s\n',name,t);
    else
        fprintf('%-46s : never\n',name);
    end

end

%% GAIN ANALYSIS
% ========================================================================
%  PBMRAC ROBUSTNESS STUDY:
%  ASYMMETRY CONTINUATION x ADAPTIVE GAIN TUNING
%
%  PURPOSE
%  -------
%  Separate two possible causes of the catastrophic asymmetric-rotor cases:
%
%    A) STRUCTURAL MODEL MISMATCH
%       The PBMRAC assumes a diagonal translational effectiveness Lambda_f,
%       while unequal rotor efficiencies introduce force/moment coupling.
%
%    B) OVER-AGGRESSIVE ADAPTATION / PREDICTOR TUNING
%       Large L, Gamma_theta and Gamma_lambda may react too aggressively to
%       unmodelled dynamics, actuator lag and saturation.
%
%
%  ASYMMETRY CONTINUATION
%  ----------------------
%
%      lambda_i(alpha)
%          = lambda_mean + alpha*(lambda_i_run9 - lambda_mean)
%
%      alpha = 0  -> perfectly uniform rotor degradation
%      alpha = 1  -> exact catastrophic MC run 9
%
%  Mean rotor efficiency is IDENTICAL for every alpha.
%
%
%  IMPORTANT
%  ---------
%  - Projection is unchanged.
%  - e-modification is unchanged.
%  - sigma_e remains 0.02.
%  - NO sigma-modification is introduced.
%  - This is a DIAGNOSTIC study, not final controller tuning.
%
%  ========================================================================

clearvars;
close all;
clc;

MC_FILE = 'results/montecarlo_20260919_163739.mat';
MODEL   = 'Final';

MC_INDEX = 9;

% Diagnostic duration.
%
% Run 9 becomes obviously bad by ~10-16 s, so 40 s is sufficient to detect
% the robustness boundary without running every case for 100 s.
%
% Once a promising tuning is identified, validate it again over the full
% final trajectory duration.
T_DIAG = 40;

% Asymmetry continuation
ALPHA = 0:0.1:1;

% Actual planner collective limit
T_CAP = 22.9308;                          % [N]

% Actual rotor upper limit
OMEGA_TOL = 0.999;


% ========================================================================
%  1. LOAD EXACT BAD MONTE CARLO CASE
%  ========================================================================

load(MC_FILE,'MC');

k = MC_INDEX;

m_pkg       = MC.samples.mass(k);
r_c_pkg     = MC.samples.r_c(k,:).';
lambda_bad  = MC.samples.lambda(k,:).';

lambda_mean = mean(lambda_bad);

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' PBMRAC ASYMMETRY / GAIN ROBUSTNESS STUDY\n');
fprintf('=====================================================================\n');

fprintf('MC sample                  : %d\n',k);
fprintf('Package mass               : %.6f kg\n',m_pkg);

fprintf('Package CoM                : [%+.5f %+.5f %+.5f] m\n', ...
    r_c_pkg(1),r_c_pkg(2),r_c_pkg(3));

fprintf('Run-9 rotor efficiencies   : [');
fprintf(' %.4f',lambda_bad);
fprintf(' ]\n');

fprintf('Mean rotor efficiency      : %.6f\n',lambda_mean);

fprintf('Diagnostic duration        : %.1f s\n',T_DIAG);

fprintf('alpha                       : ');
fprintf(' %.2f',ALPHA);
fprintf('\n');


% ========================================================================
%  2. REFERENCE
%
%  Use exactly the same OLD helix stored in the Monte Carlo campaign.
%  Restrict it to T_DIAG seconds.
%  ========================================================================

idx_ref = MC.ref.time <= T_DIAG;

ref = MC.ref;

ref.time = MC.ref.time(idx_ref);

ref.signals.values = MC.ref.signals.values(idx_ref,:);


% ========================================================================
%  3. CURRENT FRIEND TUNING
%  ========================================================================

L_CURRENT = 150;

GTH_CURRENT = diag([500 500 200]) * 0.8 * 1.2;

% Equivalent:
% diag([480 480 192])

GLAM_CURRENT = [50;
                50;
                10] * 0.1;

% Equivalent:
% [5; 5; 1]

SIGMA_E_CURRENT = 0.02;


fprintf('\nCurrent adaptive tuning:\n');

fprintf('  L_p = L_v = %.3f\n',L_CURRENT);

fprintf('  Gamma_theta diag = [%g %g %g]\n', ...
    diag(GTH_CURRENT));

fprintf('  Gamma_lambda     = [%g %g %g]\n', ...
    GLAM_CURRENT);

fprintf('  sigma_e          = %.4f  (e-modification)\n', ...
    SIGMA_E_CURRENT);


% ========================================================================
%  4. GAIN CASES
%
%  The first cases are ONE-FACTOR-AT-A-TIME tests.
%  This is important: if we change everything together and performance
%  improves, we cannot tell WHICH gain caused the improvement.
%
%  Then we add a few combined lower-bandwidth configurations.
%
%  Columns:
%
%     [ L_scale, Gamma_theta_scale, Gamma_lambda_scale ]
%
%  ========================================================================

cases = struct([]);

n = 0;

% ------------------------------------------------------------------------
% Reference / current tuning
% ------------------------------------------------------------------------

n=n+1;
cases(n).name = 'CURRENT';
cases(n).Ls   = 1.00;
cases(n).Gths = 1.00;
cases(n).Gls  = 1.00;


% ------------------------------------------------------------------------
% Gamma_lambda ONLY
% ------------------------------------------------------------------------

n=n+1;
cases(n).name = 'Glam x0.50';
cases(n).Ls   = 1.00;
cases(n).Gths = 1.00;
cases(n).Gls  = 0.50;

n=n+1;
cases(n).name = 'Glam x0.25';
cases(n).Ls   = 1.00;
cases(n).Gths = 1.00;
cases(n).Gls  = 0.25;

n=n+1;
cases(n).name = 'Glam x0.10';
cases(n).Ls   = 1.00;
cases(n).Gths = 1.00;
cases(n).Gls  = 0.10;


% ------------------------------------------------------------------------
% Gamma_theta ONLY
% ------------------------------------------------------------------------

n=n+1;
cases(n).name = 'Gtheta x0.50';
cases(n).Ls   = 1.00;
cases(n).Gths = 0.50;
cases(n).Gls  = 1.00;

n=n+1;
cases(n).name = 'Gtheta x0.25';
cases(n).Ls   = 1.00;
cases(n).Gths = 0.25;
cases(n).Gls  = 1.00;

n=n+1;
cases(n).name = 'Gtheta x0.10';
cases(n).Ls   = 1.00;
cases(n).Gths = 0.10;
cases(n).Gls  = 1.00;


% ------------------------------------------------------------------------
% Predictor gains L ONLY
% ------------------------------------------------------------------------

n=n+1;
cases(n).name = 'L x0.50';
cases(n).Ls   = 0.50;
cases(n).Gths = 1.00;
cases(n).Gls  = 1.00;

n=n+1;
cases(n).name = 'L x0.25';
cases(n).Ls   = 0.25;
cases(n).Gths = 1.00;
cases(n).Gls  = 1.00;

n=n+1;
cases(n).name = 'L x0.10';
cases(n).Ls   = 0.10;
cases(n).Gths = 1.00;
cases(n).Gls  = 1.00;


% ------------------------------------------------------------------------
% COMBINED lower-bandwidth cases
% ------------------------------------------------------------------------

n=n+1;
cases(n).name = 'ALL x0.50';
cases(n).Ls   = 0.50;
cases(n).Gths = 0.50;
cases(n).Gls  = 0.50;

n=n+1;
cases(n).name = 'ALL x0.25';
cases(n).Ls   = 0.25;
cases(n).Gths = 0.25;
cases(n).Gls  = 0.25;

n=n+1;
cases(n).name = 'ALL x0.10';
cases(n).Ls   = 0.10;
cases(n).Gths = 0.10;
cases(n).Gls  = 0.10;


NCASE  = numel(cases);
NALPHA = numel(ALPHA);


% ========================================================================
%  5. PRINT EXACT GAIN TABLE
%  ========================================================================

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' GAIN CASES\n');
fprintf('=====================================================================\n');

fprintf('%-16s %8s   %-24s %-20s\n', ...
    'case','L','diag(Gamma_theta)','Gamma_lambda');

fprintf('%s\n',repmat('-',1,80));

for c = 1:NCASE

    Lc   = L_CURRENT * cases(c).Ls;
    Gthc = GTH_CURRENT * cases(c).Gths;
    Glc  = GLAM_CURRENT * cases(c).Gls;

    fprintf('%-16s %8.2f   [%6.1f %6.1f %6.1f]   [%5.2f %5.2f %5.2f]\n', ...
        cases(c).name, ...
        Lc, ...
        diag(Gthc), ...
        Glc);

end


% ========================================================================
%  6. PREALLOCATE RESULTS
%  ========================================================================

shape = [NALPHA NCASE];

RMS_P        = nan(shape);
PEAK_P       = nan(shape);
SS_P         = nan(shape);

MAX_TILT     = nan(shape);
MAX_OMEGA_B  = nan(shape);

MIN_LAMBDA   = nan(shape);

TC_CAP_FRAC  = nan(shape);

ACT_UP_FRAC  = nan(shape);

MAX_FC       = nan(shape);

SIM_OK       = false(shape);


% ========================================================================
%  7. BASELINE SWEEP
%
%  Baseline does NOT depend on adaptive gains, so run it only once per alpha.
%  This tells us whether the physical asymmetric plant itself remains
%  controllable.
%  ========================================================================

BASE_RMS       = nan(NALPHA,1);
BASE_PEAK      = nan(NALPHA,1);
BASE_TILT      = nan(NALPHA,1);
BASE_ACTUP     = nan(NALPHA,1);
BASE_OK        = false(NALPHA,1);


fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' BASELINE ASYMMETRY CONTINUATION\n');
fprintf('=====================================================================\n');


for ia = 1:NALPHA

    alpha = ALPHA(ia);

    lambda_i = lambda_mean + ...
               alpha*(lambda_bad - lambda_mean);

    lamSig = disturb.fault( ...
        ref.time, ...
        t0 = 0, ...
        lambda = lambda_i);

    P = harness.config( ...
        package_mass = m_pkg, ...
        package_r_c  = r_c_pkg);

    try

        Rb = harness.scenario( ...
            ref, ...
            config  = P, ...
            package = 1, ...
            mrac    = 0, ...
            lambda  = lamSig, ...
            model   = MODEL, ...
            label   = sprintf('Baseline alpha %.2f',alpha));

        m = harness.metrics(Rb);

        BASE_RMS(ia)  = m.rms_p;
        BASE_PEAK(ia) = m.peak_p;

        if isfield(Rb,'euler')

            E = rad2deg(Rb.euler);

            % ZYX log: [yaw pitch roll]
            BASE_TILT(ia) = max( ...
                max(abs(E(:,2:3)),[],2) );

        end

        BASE_ACTUP(ia) = actuatorUpperFraction( ...
            Rb, ...
            P.UAV, ...
            lambda_i);

        BASE_OK(ia) = true;

        fprintf( ...
            'Baseline alpha=%4.2f | RMS=%8.4f m | peak=%8.4f m | tilt=%7.2f deg | sat=%6.1f%%\n', ...
            alpha, ...
            BASE_RMS(ia), ...
            BASE_PEAK(ia), ...
            BASE_TILT(ia), ...
            100*BASE_ACTUP(ia));

    catch ME

        fprintf( ...
            'Baseline alpha=%4.2f FAILED numerically: %s\n', ...
            alpha, ...
            ME.message);

    end

end


% ========================================================================
%  8. ADAPTIVE GAIN x ASYMMETRY SWEEP
%  ========================================================================

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' PBMRAC GAIN x ASYMMETRY SWEEP\n');
fprintf('=====================================================================\n');


totalRuns = NCASE*NALPHA;
runCount  = 0;


for c = 1:NCASE

    fprintf('\n-------------------------------------------------------------\n');
    fprintf('CASE %d/%d : %s\n',c,NCASE,cases(c).name);
    fprintf('-------------------------------------------------------------\n');

    for ia = 1:NALPHA

        runCount = runCount + 1;

        alpha = ALPHA(ia);

        % ---------------------------------------------------------------
        % Exact asymmetry interpolation
        % ---------------------------------------------------------------

        lambda_i = lambda_mean + ...
                   alpha*(lambda_bad - lambda_mean);

        % Numerical sanity check: mean MUST remain fixed.
        if abs(mean(lambda_i)-lambda_mean) > 1e-12
            error('Mean efficiency changed during alpha continuation.');
        end


        % ---------------------------------------------------------------
        % Build plant
        % ---------------------------------------------------------------

        P = harness.config( ...
            package_mass = m_pkg, ...
            package_r_c  = r_c_pkg);


        % ---------------------------------------------------------------
        % Override ONLY adaptive tuning.
        %
        % Do NOT rebuild design.adaptive here.
        %
        % This preserves:
        %   - current projection settings
        %   - current e-modification implementation
        %   - controller believed mass
        %   - every other model setting
        % ---------------------------------------------------------------

        P.adapt.L_p = L_CURRENT * cases(c).Ls;
        P.adapt.L_v = L_CURRENT * cases(c).Ls;

        P.adapt.Gamma_theta = ...
            GTH_CURRENT * cases(c).Gths;

        P.adapt.Gamma_lambda = ...
            GLAM_CURRENT * cases(c).Gls;


        % ---------------------------------------------------------------
        % Keep e-modification fixed
        % ---------------------------------------------------------------

        if isfield(P.adapt,'sigma_e')

            P.adapt.sigma_e = SIGMA_E_CURRENT;

        else

            warning(['P.adapt.sigma_e does not exist in this copy of ' ...
                     'harness.config. Check that you are using the current ' ...
                     'projection + e-modification version of the project.']);

        end


        % ---------------------------------------------------------------
        % Check that projection fields remain present.
        % ---------------------------------------------------------------

        projection_fields = { ...
            'lambda_proj_center', ...
            'lambda_proj_r0', ...
            'lambda_proj_r1'};

        for jj = 1:numel(projection_fields)

            if ~isfield(P.adapt,projection_fields{jj})

                warning( ...
                    'Missing projection field P.adapt.%s', ...
                    projection_fields{jj});

            end

        end


        % ---------------------------------------------------------------
        % Rotor uncertainty signal
        % ---------------------------------------------------------------

        lamSig = disturb.fault( ...
            ref.time, ...
            t0 = 0, ...
            lambda = lambda_i);


        % ---------------------------------------------------------------
        % Run
        % ---------------------------------------------------------------

        try

            R = harness.scenario( ...
                ref, ...
                config  = P, ...
                package = 1, ...
                mrac    = 1, ...
                robust  = 1, ...
                lambda  = lamSig, ...
                model   = MODEL, ...
                label   = sprintf('%s alpha %.2f', ...
                                  cases(c).name,alpha));


            % -----------------------------------------------------------
            % Tracking metrics
            % -----------------------------------------------------------

            m = harness.metrics(R);

            RMS_P(ia,c)  = m.rms_p;
            PEAK_P(ia,c) = m.peak_p;
            SS_P(ia,c)   = m.ss_p;


            % -----------------------------------------------------------
            % Attitude
            % -----------------------------------------------------------

            if isfield(R,'euler')

                E = rad2deg(R.euler);

                % ZYX = [yaw pitch roll]
                tilt = max(abs(E(:,2:3)),[],2);

                MAX_TILT(ia,c) = max(tilt);

            end


            % -----------------------------------------------------------
            % Body angular velocity
            % -----------------------------------------------------------

            if isfield(R,'omegab')

                MAX_OMEGA_B(ia,c) = ...
                    max(vecnorm(R.omegab,2,2));

            end


            % -----------------------------------------------------------
            % Adaptive Lambda state
            % -----------------------------------------------------------

            if isfield(R,'lambda_hat')

                MIN_LAMBDA(ia,c) = ...
                    min(R.lambda_hat(:));

            end


            % -----------------------------------------------------------
            % Planner collective saturation
            % -----------------------------------------------------------

            if isfield(R,'T_c')

                TC_CAP_FRAC(ia,c) = ...
                    mean(R.T_c(:,1) >= 0.99*T_CAP);

            end


            % -----------------------------------------------------------
            % Physical upper actuator-limit occupancy
            % -----------------------------------------------------------

            ACT_UP_FRAC(ia,c) = ...
                actuatorUpperFraction( ...
                    R, ...
                    P.UAV, ...
                    lambda_i);


            % -----------------------------------------------------------
            % Outer-loop force magnitude
            % -----------------------------------------------------------

            if isfield(R,'f_c')

                MAX_FC(ia,c) = ...
                    max(vecnorm(R.f_c,2,2));

            end


            SIM_OK(ia,c) = ...
                all(isfinite([ ...
                    RMS_P(ia,c), ...
                    PEAK_P(ia,c), ...
                    MAX_TILT(ia,c)]));


            fprintf( ...
                '[%3d/%3d] %-15s alpha=%4.2f | RMS=%8.4f | peak=%8.3f | tilt=%7.2f | minLam=%6.3f | TcCap=%6.1f%% | rotorCap=%6.1f%%\n', ...
                runCount,totalRuns, ...
                cases(c).name, ...
                alpha, ...
                RMS_P(ia,c), ...
                PEAK_P(ia,c), ...
                MAX_TILT(ia,c), ...
                MIN_LAMBDA(ia,c), ...
                100*TC_CAP_FRAC(ia,c), ...
                100*ACT_UP_FRAC(ia,c));


        catch ME

            fprintf( ...
                '[%3d/%3d] %-15s alpha=%4.2f | SIMULATION ERROR: %s\n', ...
                runCount,totalRuns, ...
                cases(c).name, ...
                alpha, ...
                ME.message);

        end

    end
end


% ========================================================================
%  9. SUMMARY TABLE
%  ========================================================================

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' PERFORMANCE AT FULL RUN-9 ASYMMETRY, alpha = 1\n');
fprintf('=====================================================================\n');

fprintf('%-16s %10s %10s %10s %10s %10s %10s\n', ...
    'case', ...
    'RMS[m]', ...
    'peak[m]', ...
    'tilt[deg]', ...
    'minLam', ...
    'TcCap[%]', ...
    'rotCap[%]');

fprintf('%s\n',repmat('-',1,95));

ia1 = NALPHA;

for c = 1:NCASE

    fprintf('%-16s %10.4f %10.3f %10.2f %10.3f %10.1f %10.1f\n', ...
        cases(c).name, ...
        RMS_P(ia1,c), ...
        PEAK_P(ia1,c), ...
        MAX_TILT(ia1,c), ...
        MIN_LAMBDA(ia1,c), ...
        100*TC_CAP_FRAC(ia1,c), ...
        100*ACT_UP_FRAC(ia1,c));

end


% ========================================================================
%  10. ESTIMATE ROBUSTNESS BOUNDARIES
%
%  We DO NOT define one arbitrary "failure" criterion.
%
%  Instead print several physically interpretable thresholds separately.
%  ========================================================================

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' ASYMMETRY THRESHOLD DIAGNOSTICS\n');
fprintf('=====================================================================\n');

fprintf('%-16s %12s %12s %12s %12s\n', ...
    'case', ...
    'RMS>1m', ...
    'peak>1m', ...
    'tilt>30deg', ...
    'Lam<0.45');

fprintf('%s\n',repmat('-',1,70));


for c = 1:NCASE

    a_rms  = firstAlpha(ALPHA,RMS_P(:,c)      > 1.0);
    a_peak = firstAlpha(ALPHA,PEAK_P(:,c)     > 1.0);
    a_tilt = firstAlpha(ALPHA,MAX_TILT(:,c)   > 30);
    a_lam  = firstAlpha(ALPHA,MIN_LAMBDA(:,c) < 0.45);

    fprintf('%-16s %12s %12s %12s %12s\n', ...
        cases(c).name, ...
        fmtAlpha(a_rms), ...
        fmtAlpha(a_peak), ...
        fmtAlpha(a_tilt), ...
        fmtAlpha(a_lam));

end


% ========================================================================
%  11. FIGURE - RMS POSITION ERROR vs ASYMMETRY
%  ========================================================================

figure( ...
    'Name','Asymmetry continuation - RMS position error', ...
    'Color','w', ...
    'Position',[80 80 1100 650]);

hold on;
grid on;

plot(ALPHA,BASE_RMS, ...
    'k--', ...
    'LineWidth',2.0, ...
    'DisplayName','Baseline');

for c = 1:NCASE

    plot(ALPHA,RMS_P(:,c), ...
        '-o', ...
        'LineWidth',1.3, ...
        'MarkerSize',4, ...
        'DisplayName',cases(c).name);

end

set(gca,'YScale','log');

xlabel('asymmetry factor alpha');
ylabel('RMS position error [m]');

title('Tracking sensitivity to rotor-to-rotor asymmetry');

legend('Location','eastoutside');


% ========================================================================
%  12. FIGURE - MAX ATTITUDE EXCURSION
%  ========================================================================

figure( ...
    'Name','Asymmetry continuation - attitude excursion', ...
    'Color','w', ...
    'Position',[100 100 1100 650]);

hold on;
grid on;

plot(ALPHA,BASE_TILT, ...
    'k--', ...
    'LineWidth',2.0, ...
    'DisplayName','Baseline');

for c = 1:NCASE

    plot(ALPHA,MAX_TILT(:,c), ...
        '-o', ...
        'LineWidth',1.3, ...
        'MarkerSize',4, ...
        'DisplayName',cases(c).name);

end

yline(30,'k:','30 deg diagnostic');

xlabel('asymmetry factor alpha');
ylabel('max roll/pitch excursion [deg]');

title('Attitude robustness vs rotor asymmetry');

legend('Location','eastoutside');


% ========================================================================
%  13. FIGURE - LAMBDA PROJECTION ACTIVITY
%  ========================================================================

figure( ...
    'Name','Asymmetry continuation - minimum Lambda hat', ...
    'Color','w', ...
    'Position',[120 120 1100 650]);

hold on;
grid on;

for c = 1:NCASE

    plot(ALPHA,MIN_LAMBDA(:,c), ...
        '-o', ...
        'LineWidth',1.3, ...
        'MarkerSize',4, ...
        'DisplayName',cases(c).name);

end

yline(0.45,'k--','projection transition');
yline(0.35,'k:','outer projection wall');

xlabel('asymmetry factor alpha');
ylabel('minimum Lambda hat');

title('Adaptive-state response to increasing structural mismatch');

legend('Location','eastoutside');


% ========================================================================
%  14. FIGURE - ACTUATOR-LIMIT OCCUPANCY
%  ========================================================================

figure( ...
    'Name','Asymmetry continuation - actuator saturation', ...
    'Color','w', ...
    'Position',[140 140 1100 650]);

hold on;
grid on;

plot(ALPHA,100*BASE_ACTUP, ...
    'k--', ...
    'LineWidth',2.0, ...
    'DisplayName','Baseline');

for c = 1:NCASE

    plot(ALPHA,100*ACT_UP_FRAC(:,c), ...
        '-o', ...
        'LineWidth',1.3, ...
        'MarkerSize',4, ...
        'DisplayName',cases(c).name);

end

xlabel('asymmetry factor alpha');
ylabel('upper actuator-limit occupancy [%]');

title('Physical actuator saturation vs rotor asymmetry');

legend('Location','eastoutside');


% ========================================================================
%  15. FIGURE - PLANNER COLLECTIVE CAP
%  ========================================================================

figure( ...
    'Name','Asymmetry continuation - collective cap', ...
    'Color','w', ...
    'Position',[160 160 1100 650]);

hold on;
grid on;

for c = 1:NCASE

    plot(ALPHA,100*TC_CAP_FRAC(:,c), ...
        '-o', ...
        'LineWidth',1.3, ...
        'MarkerSize',4, ...
        'DisplayName',cases(c).name);

end

xlabel('asymmetry factor alpha');
ylabel('time at >=99% collective cap [%]');

title('Adaptive command aggressiveness vs rotor asymmetry');

legend('Location','eastoutside');


% ========================================================================
%  16. SAVE RESULTS
%  ========================================================================

Study = struct();

Study.alpha = ALPHA;
Study.cases = cases;

Study.package_mass = m_pkg;
Study.package_r_c  = r_c_pkg;

Study.lambda_bad  = lambda_bad;
Study.lambda_mean = lambda_mean;

Study.current.L_p          = L_CURRENT;
Study.current.L_v          = L_CURRENT;
Study.current.Gamma_theta  = GTH_CURRENT;
Study.current.Gamma_lambda = GLAM_CURRENT;
Study.current.sigma_e      = SIGMA_E_CURRENT;

Study.rms_p       = RMS_P;
Study.peak_p      = PEAK_P;
Study.ss_p        = SS_P;

Study.max_tilt    = MAX_TILT;
Study.max_omega_b = MAX_OMEGA_B;

Study.min_lambda  = MIN_LAMBDA;

Study.tc_cap_frac = TC_CAP_FRAC;
Study.act_up_frac = ACT_UP_FRAC;

Study.max_fc      = MAX_FC;
Study.sim_ok      = SIM_OK;

Study.baseline.rms_p       = BASE_RMS;
Study.baseline.peak_p      = BASE_PEAK;
Study.baseline.max_tilt    = BASE_TILT;
Study.baseline.act_up_frac = BASE_ACTUP;
Study.baseline.sim_ok      = BASE_OK;

outfile = fullfile( ...
    'results', ...
    ['pbmrac_asymmetry_gainstudy_' ...
     datestr(now,'yyyymmdd_HHMMSS') '.mat']);

save(outfile,'Study');

fprintf('\nSaved study to:\n%s\n',outfile);


% ========================================================================
%  17. WHAT TO SEND BACK
%  ========================================================================

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' SEND ME THESE RESULTS\n');
fprintf('=====================================================================\n');

fprintf(['1) PERFORMANCE AT FULL RUN-9 ASYMMETRY table\n' ...
         '2) ASYMMETRY THRESHOLD DIAGNOSTICS table\n' ...
         '3) Any obvious change in where RMS / tilt / saturation blow up\n' ...
         '4) If possible, screenshots of the first four figures\n']);


% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function frac = actuatorUpperFraction(R,UAV,lambda)

    if ~isfield(R,'T_out')

        frac = NaN;
        return

    end

    lambda = reshape(lambda,1,[]);

    % T_out is downstream of rotor effectiveness.
    %
    % Undo lambda to recover post-actuator-limit rotor thrust.
    Tpre = R.T_out ./ lambda;

    Omega = sqrt(abs(Tpre)/UAV.k_f);

    atUpper = any( ...
        Omega >= 0.999*UAV.Omega_max, ...
        2);

    frac = mean(atUpper);

end


function a = firstAlpha(alpha,mask)

    idx = find(mask,1,'first');

    if isempty(idx)
        a = NaN;
    else
        a = alpha(idx);
    end

end


function s = fmtAlpha(a)

    if isnan(a)
        s = 'none';
    else
        s = sprintf('%.2f',a);
    end

end

%%
load('results/montecarlo_20260919_163739.mat','MC');

k = 9;

%% Exact old catastrophic physical sample
P9 = harness.config( ...
    package_mass = MC.samples.mass(k), ...
    package_r_c  = MC.samples.r_c(k,:).');

%% Force the final candidate gains explicitly
P9.adapt.L_p = 150;
P9.adapt.L_v = 150;

P9.adapt.Gamma_theta = diag([480 480 192]) * 0.30;
P9.adapt.Gamma_lambda = [5; 5; 1] * 0.15;

P9.adapt.sigma_e = 0.02;   % e-modification on theta

%% Exact old MC rotor degradation
lam9 = disturb.fault( ...
    MC.ref.time, ...
    t0 = 0, ...
    lambda = MC.samples.lambda(k,:).');

%% Replay on OLD Monte Carlo helix
R9 = harness.scenario( ...
    MC.ref, ...
    config  = P9, ...
    package = 1, ...
    mrac    = 1, ...
    robust  = 1, ...
    lambda  = lam9, ...
    model   = 'Final', ...
    label   = 'MC run 9 - robust PBMRAC');

viz.errors({R9},{'MC run 9'}, ...
    name='MC run 9 - tracking errors');

viz.actuators({R9},{'MC run 9'},P9.UAV, ...
    name='MC run 9 - actuators');

viz.attitude({R9},{'MC run 9'}, ...
    name='MC run 9 - attitude');

viz.effectiveness({R9},{'MC run 9'}, ...
    lambda=lam9, ...
    t_fault=0, ...
    name='MC run 9 - adaptive effectiveness');

harness.report({R9},{'MC run 9'}, ...
    title='Monte Carlo run 9 diagnostic');
%% NEW MONTECARLO
MC = run_montecarlo_final(n=500);
%% ---- Save the figures ------------------------------------------------
% Written to Figures/Task3_package/ , named from each figure's Name property and
% numbered so the folder sorts in slide order.
viz.export('Task3_dummy_montecarlo_check');
