%% Spacecraft Guidance and Navigation (2025/2026)
% Assignment 1 – Exercise 3: Impulsive Transfer to Distant Prograde Orbit
% Author : Karim Ahmed Mohamed Doubie
%
% Objective:
%   Point 1 – Propagate the spacecraft state from t0 to tA using an
%             N-body model centred on Earth (J2000 frame) and compare
%             with the DPO ephemeris at arrival.
%   Point 2 – Formal problem statement (see report).
%   Point 3 – Solve the multiple-shooting impulsive transfer optimisation
%             via fmincon. The decision vector contains five impulses,
%             three intermediate position nodes, the final node, and the
%             free final time. The objective minimises total Delta-V squared.
%   Point 3 (Alternative) – Solve the same transfer problem with an alternative
%             multiple-shooting formulation adopting a full-state node decision
%             vector. The objective minimises total Delta-V squared.
%% ====================== Initialisation ==================================
% Step #0.1: Environment setup and figure LaTeX rendering options
close all; clearvars; clc;
format long;


% ── Global plot defaults ──────────────────────────────────────────────────
set(groot, 'DefaultAxesFontSize',                 16);
set(groot, 'DefaultAxesTitleFontSizeMultiplier',   1);
set(groot, 'DefaultTextFontSize',                 12);
set(groot, 'DefaultLineLineWidth',               1.5);
set(0,     'DefaultTextInterpreter',         'latex');
set(0,     'DefaultLegendInterpreter',       'latex');
set(0,     'DefaultAxesTickLabelInterpreter','latex');
% Step #0.2: Load necessary SPICE kernels
% --- SPICE kernels ---
cspice_kclear();
cspice_furnsh('assignment01.tm');

% Step #0.3: Define mission epochs and state vectors at t0 and list of bodies
% --- Mission epochs ---
t0_str = '2032-11-14T05:16:23.639';
tA_str = '2033-02-14T07:45:05.552';
et0    = cspice_str2et(t0_str);
etf    = cspice_str2et(tA_str);      

% --- Initial state (Earth-centred J2000, km / km/s) ---
r0 = [4882.128777349525;
      4030.951188974080;
      1444.063996364637];
v0 = [-5.761404259200368;
       9.117460162775972;
      -2.399524966790392];
x0 = [r0; v0];
% --- Third-body perturber list (Earth excluded) ---
bodies = {'Mars Barycenter'; 'Jupiter Barycenter'; 'Saturn Barycenter'; ...
          'Sun'; 'Mercury'; 'Venus'; 'Moon'; ...
          'Uranus Barycenter'; 'Neptune Barycenter'; 'Pluto Barycenter'};

% --- ODE integration options ---
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
                 'MaxStep', (etf - et0) / (3*500));
             
%% ============================ Point 1 ==================================
fprintf('\n==================== Point 1 ====================\n');

% Step #1.1: Propagate the spacecraft state from t0 to tA with N-body dynamics
frame     = 'J2000';
earth_id  = 399;
target_id = -1032046;
GM_Earth = cspice_bodvrd('EARTH', 'GM', 1);

[tt, xx] = ode113(@(t,x) nbody_shift_rhs(t, x, bodies, frame, earth_id, GM_Earth), ...
                  [et0, etf], x0, options);

% Step #1.2: Retrieve the nominal target DPO state at the arrival epoch
% Target DPO state at arrival
rv_dpo    = cspice_spkez(target_id, etf, frame, 'NONE', earth_id);
%target ephemeris x={r,v}
r_dpo_p1 = rv_dpo(1:3);
  %miss distance at tA
  dr_p1 = xx(end,1:3).' - r_dpo_p1;
v_dpo_p1 = rv_dpo(4:6);
  % miss velocity at tA
  dv_p1 = xx(end,4:6).' - v_dpo_p1;

lbl_i = 'State Vector Error Analysis';
met_i = ["Pos Error x [km]"; "Pos Error y [km]"; "Pos Error z [km]"; ...
         "Vel Error x [km/s]"; "Vel Error y [km/s]"; "Vel Error z [km/s]"; ...
         "||r(tA) - r_DPO|| [km]"; "||v(tA) - v_DPO|| [m/s]"];
val_i = [dr_p1(1); dr_p1(2); dr_p1(3); ...
         dv_p1(1); dv_p1(2); dv_p1(3); ...
         norm(dr_p1); norm(dv_p1)*1e3];

% Step #1.4: Output performance tables and visualization
fprintf('\n=== %s ===\n', lbl_i);
disp(table(met_i, val_i, 'VariableNames', {'Metric', 'Value'}));

% Figure 1: propagated arc vs DPO position at tA
    n_pts       = 500;
    t_grid      = linspace(et0, etf, n_pts);
    rr_moon     = zeros(3, n_pts);
    rr_dpo_arc  = nan(3, n_pts);
    
    for k = 1:n_pts
        rr_moon(:,k) = cspice_spkpos('Moon', t_grid(k), frame, 'NONE', num2str(earth_id));
    end
    
    figure; hold on; grid on; axis equal;
    
    plot3(rr_moon(1,:)*1e-5, rr_moon(2,:)*1e-5, rr_moon(3,:)*1e-5, ...
          '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
          'DisplayName', '$\mathrm{Moon\ Orbit}$');
    
    plot3(xx(:,1)*1e-5, xx(:,2)*1e-5, xx(:,3)*1e-5, ...
          'Color', [0.85 0.33 0.1], 'LineWidth', 1.5, ...
          'DisplayName', '$\mathrm{Spacecraft\ Trajectory}$');
    
    plot3(0, 0, 0, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8, ...
          'DisplayName', '$\mathrm{Earth}$');
    
    plot3(r_dpo_p1(1)*1e-5, r_dpo_p1(2)*1e-5, r_dpo_p1(3)*1e-5, ...
          's', 'Color', [0.0 0.6 0.6], 'MarkerFaceColor', [0.0 0.6 0.6], 'MarkerSize', 7, ...
          'DisplayName', '$\mathrm{DPO\ at\ }t_A$');
    
    xlabel('$x\ [\times10^5\ \mathrm{km}]$');
    ylabel('$y\ [\times10^5\ \mathrm{km}]$');
    zlabel('$z\ [\times10^5\ \mathrm{km}]$');
    title('$\mathrm{N\mbox{-}body\ Dynamics\ Propagation\ in\ J2000\ Frame}$');
    legend('Location', 'best');
    view(3);

%% ============================ Point 3 Primary Solution ==================================
fprintf('\n==================== Point 3 ====================\n');

% Step #3.1: Define multiple shooting segments and boundary conditions
% --- Multiple-shooting parameters ---
N   = 5;                              
tA_ub_str = '2033-02-14T07:48:08.552';
etf_ub    = cspice_str2et(tA_ub_str);
dT_nom    = etf    - et0;             % nominal ToF [s]
dT_ub     = etf_ub - et0;            % upper-bound ToF [s]

% Step #3.2: Formulate node state warm-starts by sequential unforced propagation
% % --- Build uniform initial time grid and propagate unforced arcs ---
ti   = et0 + (0:N-1) * dT_nom/(N-1);   % N points, same spacing as original
xi   = zeros(6, N-1);
curr = x0;
for k = 1:N-1
    [~, Xtmp] = ode113(@(t,x) nbody_shift_rhs(t,x,bodies,frame,earth_id,GM_Earth), ...
                       [ti(k) ti(k+1)], curr, options);
    curr    = Xtmp(end,:).';
    xi(:,k) = curr;
end
% --- Initial impulse guesses (zero — unforced arcs are continuous) ---
dv_guess = zeros(3, N);    % [3 x N]:  dv at nodes 0..4

% Step #3.3: Flatten the optimization variables into the decision vector
% --- Flatten optimisation decision vector ---
dv0g = dv_guess(:,1);  dv1g = dv_guess(:,2);
dv2g = dv_guess(:,3);  dv3g = dv_guess(:,4);  dv4g = dv_guess(:,5);

r1g = xi(1:3,1);  r2g = xi(1:3,2);
r3g = xi(1:3,3);  rfg = xi(1:3,4);


guess_0 = [dv0g; r1g; dv1g; r2g; dv2g; r3g; dv3g; rfg; dv4g; dT_nom];


% Step #3.4: Set up objective function and optimization bounds
% --- Objective: minimise sum of squared impulse norms (total DeltaV^2) ---
deltaV_total = @(z) sum(vecnorm([z(1:3), z(7:9), z(13:15), z(19:21), z(25:27)], 2, 1).^2);

% --- Nonlinear constraints (defect + terminal matching) ---
nonlconst= @(z)MS_constraints(z, r0, v0, et0, frame, earth_id, bodies, target_id, GM_Earth,N);
% --- Bounds on decision vector ---
lb = -inf(28,1);  ub = inf(28,1);
lb(28) = 0.7 * dT_nom;
ub(28) = 1.2 * dT_ub;

% Step #3.5: Execute transfer optimization via fmincon
% --- fmincon problem struct ---
problem.objective   = deltaV_total;
problem.x0          = guess_0;
problem.Aineq       = [];   problem.bineq = [];
problem.Aeq         = [];   problem.beq   = [];
problem.lb          = lb;   problem.ub    = ub;
problem.nonlcon     = nonlconst;
problem.solver      = 'fmincon';
problem.options     = optimoptions('fmincon', ...
    'Display',                'iter', ...
    'ConstraintTolerance',    1e-8,   ...
    'StepTolerance',          1e-10,  ...
    'MaxFunctionEvaluations', 5e4);

fprintf('\n--- Running fmincon ---\n');
tic;
[x, ~, exitflag, output] = fmincon(problem);
toc;

% Step #3.6: Verify and unpack optimal solution parameters
% --- Check constraint residuals ---
[~, ceq_sol] = nonlconst(x);
fprintf('\nexitflag = %d   iterations = %d   max|ceq| = %.3e\n', ...
        exitflag, output.iterations, max(abs(ceq_sol)));

% --- Unpack solution vector ---
x_sol = x;                        % full 28-element solution 
dv0 = x_sol(1:3);   r1  = x_sol(4:6);
dv1 = x_sol(7:9);   r2  = x_sol(10:12);
dv2 = x_sol(13:15); r3  = x_sol(16:18);
dv3 = x_sol(19:21); rf  = x_sol(22:24);
dv4 = x_sol(25:27); dT  = x_sol(28);

% impulses: [3 x N] matrix of Delta-V vectors at each node
impulses = [dv0, dv1, dv2, dv3, dv4];

etf_opt   = et0 + dT;
total_dv  = sum(vecnorm(impulses, 2, 1)) * 1e3;   % [m/s]

% Updated time grid with optimal ToF
ti_opt    = et0 + (0:N-1) * dT/(N-1);

% Step #3.7: Format and print optimized solution summaries
tformat   = 'YYYY-MM-DDTHR:MN:SC.###::UTC';
t0c  = cspice_timout(ti_opt(1), tformat);
t1c  = cspice_timout(ti_opt(2), tformat);
t2c  = cspice_timout(ti_opt(3), tformat);
t3c  = cspice_timout(ti_opt(4), tformat);
tfc  = cspice_timout(etf_opt,   tformat);

v0_opt = v0 + dv0;
r_err  = ceq_sol(13:15) * 1e3;    % [m]
v_err  = ceq_sol(16:18);          % [km/s]

%  Optimised Transfer Summary
lbl_1 = 'Optimised Transfer Summary';
met_1 = ["tf (UTC)"; "ToF [days]"; "Total Delta-V [m/s]"; "||r(tf)-r_DPO|| [km]"; "||v(tf)-v_DPO|| [m/s]"];
val_1 = [string(tfc); string(dT/86400); string(total_dv); string(norm(r_err)*1e-3); string(norm(v_err)*1e3)];
fprintf('\n=== %s ===\n', lbl_1);
disp(table(met_1, val_1, 'VariableNames', {'Metric', 'Value'}));

%  Post-insertion velocity
lbl_2 = 'Post-insertion velocity [km/s]';
met_2 = ["v0+ x"; "v0+ y"; "v0+ z"];
val_2 = [v0_opt(1); v0_opt(2); v0_opt(3)];
fprintf('\n=== %s ===\n', lbl_2);
disp(table(met_2, val_2, 'VariableNames', {'Metric', 'Value'}));

%  Optimised impulses [m/s]
tformat_imp  = 'YYYY-MM-DDTHR:MN:SC.###::UTC';
imp_labels   = {'dv0', 'dv1', 'dv2', 'dv3', 'dv4'};
imp_times    = [ti_opt, etf_opt]; 

fprintf('\n=== Optimised Impulses [m/s] ===\n');
for k = 1:size(impulses, 2)
    tc  = cspice_timout(imp_times(k), tformat_imp);
    dvi = impulses(:, k) * 1e3;
    fprintf('%s  %s:  [%15.10f  %15.10f  %15.10f]  ||dv|| = %12.10f m/s\n', ...
            tc, imp_labels{k}, dvi(1), dvi(2), dvi(3), norm(dvi));
end


% Step #3.8: Rebuild optimised trajectory
options_altb = odeset('RelTol',1e-9,'AbsTol',1e-9);

r_nodes  = [r0,  r1,  r2,  r3];
dv_nodes = [dv0, dv1, dv2, dv3];
t_start  = ti_opt(1:end-1);
t_end    = [ti_opt(2:end), etf_opt];

traj = [];
Xs   = [];
for k = 1 : N-1
    if k == 1
        v_init = v0 + dv_nodes(:,k);
    else
        v_init = Xs(end,4:6).' + dv_nodes(:,k);
    end
    [~, Xs] = ode113(@(t,x) nbody_shift_rhs(t,x,bodies,frame,earth_id,GM_Earth), ...
                     [t_start(k), t_end(k)], [r_nodes(:,k); v_init], options_altb);
    traj = [traj; Xs];
end

rv_dpo    = cspice_spkez(target_id, etf_opt, frame, 'NONE', earth_id);
r_dpo_opt = rv_dpo(1:3);

% Figure 2: optimised 3-D transfer trajectory WITH ballistic overlay + Moon + DPO
    n_pts2      = 500;
    t_grid2     = linspace(et0, etf_opt, n_pts2);
    rr_moon2    = zeros(3, n_pts2);
    rr_dpo_arc2 = nan(3, n_pts2);
    
    for k = 1:n_pts2
        rr_moon2(:,k) = cspice_spkpos('Moon', t_grid2(k), frame, 'NONE', num2str(earth_id));
    end
    
    figure; hold on; grid on; axis equal;
    
    plot3(rr_moon2(1,:), rr_moon2(2,:), rr_moon2(3,:), ...
          '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
          'DisplayName', '$\mathrm{Moon\ Orbit}$');
    plot3(xx(:,1), xx(:,2), xx(:,3), '--', 'Color', [0 0.5 0], 'LineWidth', 1.5, ...
          'DisplayName', '$\mathrm{Ballistic\ (Point\ 1)}$');
    
    plot3(traj(:,1), traj(:,2), traj(:,3), 'Color', [0.8 0.4 0], 'LineWidth', 1.5, ...
          'DisplayName', '$\mathrm{Optimized\ Trajectory}$');
    
    plot3(r0(1), r0(2), r0(3), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6, ...
          'DisplayName', '$r_0$');
    plot3(r1(1), r1(2), r1(3), 'cd', 'MarkerFaceColor', 'c', 'MarkerSize', 6, ...
          'DisplayName', '$r_i$');
    plot3(r2(1), r2(2), r2(3), 'cd', 'MarkerFaceColor', 'c', 'MarkerSize', 6, 'HandleVisibility', 'off');
    plot3(r3(1), r3(2), r3(3), 'cd', 'MarkerFaceColor', 'c', 'MarkerSize', 6, 'HandleVisibility', 'off');
    
    plot3(r_dpo_opt(1), r_dpo_opt(2), r_dpo_opt(3), 'ko', 'MarkerFaceColor', [0.0 0.6 0.6], ...
          'MarkerSize', 6, 'DisplayName', '$\mathrm{DPO\ at\ }t_f$');
    plot3(xx(end,1), xx(end,2), xx(end,3), 's', 'Color', [0 0.5 0], ...
          'MarkerFaceColor', [0 0.5 0], 'MarkerSize', 8, 'DisplayName', '$\mathrm{Ballistic\ arrival}$');
    
    xlabel('$x\ [\mathrm{km}]$'); ylabel('$y\ [\mathrm{km}]$'); zlabel('$z\ [\mathrm{km}]$');
    title('$\mathrm{Transfer\ Comparison:\ Optimized\ vs.\ Ballistic}$');
    legend('Location', 'best');
    view(3);


%% ============================ Point 3 – Alternative  ====
fprintf('\n==================== Point 3 – Alternative ====================\n');
% Decision vector layout (28 elements):
%   z(1:3)   = v0+  (post-impulse velocity at t0)
%   z(4:9)   = [r1;v1], z(10:15) = [r2;v2], z(16:21) = [r3;v3], z(22:27) = [r4;v4]
%   z(28)    = dT [s]

t_nodes_alt = linspace(et0, etf, N+1);
z0_alt      = v0;
x_node      = [r0; v0];
for k = 2:N
    [~, Xtmp] = ode113(@(t,x) nbody_shift_rhs(t,x,bodies,frame,earth_id,GM_Earth), ...
                       [t_nodes_alt(k-1) t_nodes_alt(k)], x_node, options);
    x_node = Xtmp(end,:).';
    z0_alt = [z0_alt; x_node];
end
z0_alt = [z0_alt; dT_nom];

% --- Objective: true total Delta-V (sum of velocity jump norms) ---
obj_alt = @(z)MS_alt_cost(z, r0, v0, et0, N, frame, earth_id, bodies, target_id, GM_Earth);

% --- Constraints: position defects + terminal position match ---
nlc_alt = @(z) MS_alt_constraints(z, r0, et0, N, frame, earth_id, bodies, target_id, GM_Earth);

% --- Bounds ---
lb_alt     = -inf(28,1);
ub_alt     =  inf(28,1);
lb_alt(28) = 0.7 * dT_nom;
ub_alt(28) = 1.2 * dT_ub;

% --- fmincon options ---
opts_alt  = optimoptions('fmincon', ...
    'Algorithm', ...
    'interior-point', ...
    'Display',                'iter', ...
    'ConstraintTolerance',    1e-8,   ...
     'StepTolerance',          1e-12,  ...  
    'MaxFunctionEvaluations', 5e4);

% --- fmincon problem2 struct ---
problem2.objective   = obj_alt;
problem2.x0          = z0_alt;
problem2.Aineq       = [];   problem2.bineq = [];
problem2.Aeq         = [];   problem2.beq   = [];
problem2.lb          = lb_alt;   problem2.ub    = ub_alt;
problem2.nonlcon     = nlc_alt;
problem2.solver      = 'fmincon';
problem2.options     = opts_alt;

fprintf('\n--- Running fmincon (alternative) ---\n');
tic;
[z_alt, fval_alt, exitflag_alt, output_alt] = fmincon(problem2);
toc;

% --- Check constraint residuals ---
[~, ceq_sol_alt] = nlc_alt(z_alt);
fprintf('\nexitflag = %d   iterations = %d   max|ceq| = %.3e\n', ...
        exitflag_alt, output_alt.iterations, max(abs(ceq_sol_alt)));

% --- Unpack ---
v0_plus_alt = z_alt(1:3);
dT_alt      = z_alt(28);
etf_alt     = et0 + dT_alt;
dv0_alt     = v0_plus_alt - v0;

% Intermediate nodes
r1_alt = z_alt(4:6);   v1_alt = z_alt(7:9);
r2_alt = z_alt(10:12); v2_alt = z_alt(13:15);
r3_alt = z_alt(16:18); v3_alt = z_alt(19:21);
r4_alt = z_alt(22:24); v4_alt = z_alt(25:27);
h_alt  = dT_alt / N;
ti_a   = et0 + (0:N)*h_alt;
options_alt = odeset('RelTol',1e-10,'AbsTol',1e-12);

nodes_r = [r0,          r1_alt, r2_alt, r3_alt, r4_alt];
nodes_v = [v0_plus_alt, v1_alt, v2_alt, v3_alt, v4_alt];

impulses_alt = [dv0_alt, zeros(3,N)];
traj_alt = [];  Xa = [];

for k = 1:N
    [~, Xa] = ode113(@(t,x) nbody_shift_rhs(t,x,bodies,frame,earth_id,GM_Earth), ...
                     [ti_a(k) ti_a(k+1)], [nodes_r(:,k); nodes_v(:,k)], options_alt);
    traj_alt = [traj_alt; Xa];
    if k < N
        impulses_alt(:,k+1) = nodes_v(:,k+1) - Xa(end,4:6).';
    end
end

rv_dpo_alt       = cspice_spkez(target_id, etf_alt, frame, 'NONE', earth_id);
impulses_alt(:,N+1) = rv_dpo_alt(4:6) - Xa(end,4:6).';
total_dv_alt     = sum(vecnorm(impulses_alt,2,1)) * 1e3;

r_err_alt = Xa(end,1:3).' - rv_dpo_alt(1:3);
v_err_alt = Xa(end,4:6).' + impulses_alt(:,N+1) - rv_dpo_alt(4:6);

fprintf('\nTotal Delta-V (alt) = %.4f m/s\n', total_dv_alt);
fprintf('||r(tf)-r_DPO|| (alt) = %.4e km\n',  norm(r_err_alt));
fprintf('||v(tf)-v_DPO|| (alt) = %.4e m/s\n', norm(v_err_alt)*1e3);

tformat_alt    = 'YYYY-MM-DDTHR:MN:SC.###::UTC';
impulse_labels = {'dv0','dv1','dv2','dv3','dv4','dv5'};
impulse_times_alt = [ti_a, etf_alt];
fprintf('\n=== Optimised Impulses – Alternative [m/s] ===\n');
for k = 1:N+1
    tc  = cspice_timout(impulse_times_alt(k), tformat_alt);
    dvi = impulses_alt(:,k) * 1e3;
    fprintf('%s  %s:  [%10.4f  %10.4f  %10.4f]  ||dv|| = %10.4f m/s\n', ...
            tc, impulse_labels{k}, dvi(1), dvi(2), dvi(3), norm(dvi));
end

% --- Figure: Alternative optimised trajectory ---
n_pts3      = 500;
t_grid3     = linspace(et0, etf_alt, n_pts3);
rr_moon3    = zeros(3, n_pts3);
rr_dpo_arc3 = nan(3, n_pts3);
for k = 1:n_pts3
    rr_moon3(:,k) = cspice_spkpos('Moon', t_grid3(k), frame, 'NONE', num2str(earth_id));
    try
        rv_tmp           = cspice_spkez(target_id, t_grid3(k), frame, 'NONE', earth_id);
        rr_dpo_arc3(:,k) = rv_tmp(1:3);
    catch
    end
end

figure; hold on; grid on; axis equal;

plot3(rr_moon3(1,:), rr_moon3(2,:), rr_moon3(3,:), ...
      '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
      'DisplayName', '$\mathrm{Moon\ Orbit}$');

plot3(xx(:,1), xx(:,2), xx(:,3), '--', 'Color', [0 0.5 0], 'LineWidth', 1.5, ...
      'DisplayName', '$\mathrm{Ballistic\ (Point\ 1)}$');

plot3(traj_alt(:,1), traj_alt(:,2), traj_alt(:,3), 'Color', [0.8 0.4 0], 'LineWidth', 1.5, ...
      'DisplayName', '$\mathrm{Optimized\ Trajectory\ (Alt)}$');

plot3(r0(1), r0(2), r0(3), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6, ...
      'DisplayName', '$r_0$');
plot3(r1_alt(1), r1_alt(2), r1_alt(3), 'cd', 'MarkerFaceColor', 'c', 'MarkerSize', 6, ...
      'DisplayName', '$r_i$');
plot3(r2_alt(1), r2_alt(2), r2_alt(3), 'cd', 'MarkerFaceColor', 'c', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot3(r3_alt(1), r3_alt(2), r3_alt(3), 'cd', 'MarkerFaceColor', 'c', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot3(r4_alt(1), r4_alt(2), r4_alt(3), 'cd', 'MarkerFaceColor', 'c', 'MarkerSize', 6, 'HandleVisibility', 'off');

plot3(rv_dpo_alt(1), rv_dpo_alt(2), rv_dpo_alt(3), 'ko', 'MarkerFaceColor', [0.0 0.6 0.6], ...
      'MarkerSize', 6, 'DisplayName', '$\mathrm{DPO\ at\ }t_f$');
plot3(xx(end,1), xx(end,2), xx(end,3), 's', 'Color', [0 0.5 0], ...
      'MarkerFaceColor', [0 0.5 0], 'MarkerSize', 8, 'DisplayName', '$\mathrm{Ballistic\ arrival}$');

xlabel('$x\ [\mathrm{km}]$'); ylabel('$y\ [\mathrm{km}]$'); zlabel('$z\ [\mathrm{km}]$');
title('$\mathrm{Transfer\ Comparison:\ Optimized\ (Alternative)\ vs.\ Ballistic}$');
legend('Location', 'best');
view(3);
%% ========================== Local Functions ============================

function J = MS_alt_cost(z, r0, v0, et0, N, frame, earth_id, bodies, target_id, GM_Earth)
% AIM: Computes the total Delta-V cost for the alternative multiple-shooting
%      formulation using full-state node parameterization. The objective is
%      the true sum of velocity-jump norms (not squared).
% INPUTS:
%   z         [28x1] Decision vector [v0+{3}; r1,v1{6}; ...; r4,v4{6}; dT{1}]
%   r0        [3x1]  Initial position [km]
%   v0        [3x1]  Initial velocity [km/s]
%   et0       [1x1]  Initial epoch [s] (SPICE ET)
%   etf       [1x1]  Nominal final epoch [s] (used for MaxStep reference)
%   N         [1x1]  Number of segments
%   frame     [char] Reference frame (e.g. 'J2000')
%   earth_id  [1x1]  SPICE ID for Earth (399)
%   bodies    [cell] Perturber body names
%   target_id [1x1]  SPICE ID for DPO target
%   GM_Earth  [1x1]  Earth gravitational parameter [km^3/s^2]
% OUTPUT:
%   J         [1x1]  Total Delta-V cost [km/s]

    % Step #B.1: Unpack decision vector into node states and time grid
    v0p   = z(1:3);
    dT    = z(end);
    h     = dT / N;
    ti    = et0 + (0:N) * h;
    nodes = reshape(z(4:end-1), 6, N-1);   % [r1;v1, r2;v2, ...]

    r_nodes = [r0,  nodes(1:3,:)];          % [3 x N] arc start positions
    v_nodes = [v0p, nodes(4:6,:)];          % [3 x N] arc start velocities

    % Step #B.2: Set ODE integration options
    opts_c = odeset('RelTol',1e-10,'AbsTol',1e-12);

% Step #B.3: Accumulate initial impulse cost (dv0)
    J  = norm(v0p - v0)^2;
    xf = zeros(6,1);

    % Step #B.4: Propagate each arc and accumulate interior impulse costs
    for k = 1:N
        [~, Xk] = ode113(@(t,x) nbody_shift_rhs(t,x,bodies,frame,earth_id,GM_Earth), ...
                         [ti(k), ti(k+1)], [r_nodes(:,k); v_nodes(:,k)], opts_c);
        xf = Xk(end,:).';
        if k < N
            % Interior impulse: velocity jump at node k+1
            J = J + norm(v_nodes(:,k+1) - xf(4:6))^2;
        end
    end

    % Step #B.5: Add terminal acquisition impulse (dv5) to match DPO velocity
    rv_dpo = cspice_spkez(target_id, et0+dT, frame, 'NONE', earth_id);
    J = J + norm(rv_dpo(4:6) - xf(4:6))^2;
end


function [c, ceq] = MS_alt_constraints(z, r0, et0, N, frame, earth_id, bodies, target_id, GM_Earth)
% AIM: Formulates equality constraints for the alternative multiple-shooting
%      formulation. Only position defects are enforced as hard constraints;
%      velocity continuity at intermediate nodes is NOT required.
% INPUTS:
%   z         [28x1] Decision vector [v0+{3}; r1,v1{6}; ...; r4,v4{6}; dT{1}]
%   r0        [3x1]  Initial position [km]
%   et0       [1x1]  Initial epoch [s] (SPICE ET)
%   N         [1x1]  Number of segments
%   frame     [char] Reference frame (e.g. 'J2000')
%   earth_id  [1x1]  SPICE ID for Earth (399)
%   bodies    [cell] Perturber body names
%   target_id [1x1]  SPICE ID for DPO target
%   GM_Earth  [1x1]  Earth gravitational parameter [km^3/s^2]
% OUTPUTS:
%   c         []     Inequality constraints (none)
%   ceq       [Nx1]  Position defects at N-1 intermediate nodes + terminal match

    % Step #A.1: Unpack time variables and scaling
    dT               = z(end);
    h                = dT / N;
    position_scaling = 1000;            % [km] dimensional scaling
    opts_c           = odeset('RelTol',1e-10,'AbsTol',1e-12);

    % Step #A.2: Unpack node states from decision vector
    nodes   = reshape(z(4:end-1), 6, N-1);  % [r1;v1, r2;v2, ...]
    r_nodes = [r0,     nodes(1:3,:)];        % [3 x N] arc start positions
    v_nodes = [z(1:3), nodes(4:6,:)];        % [3 x N] arc start velocities

    % Step #A.3: Initialise constraint arrays
    c   = [];
    ceq = [];

    % Step #A.4: Propagate each arc and assemble position defect constraints
    for k = 1:N
        t0k = et0 + (k-1)*h;
        [~, Yk] = ode113(@(t,x) nbody_shift_rhs(t,x,bodies,frame,earth_id,GM_Earth), ...
                         [t0k, t0k+h], [r_nodes(:,k); v_nodes(:,k)], opts_c);
        rf_k = Yk(end,1:3).';

        if k < N
            % Interior defect: propagated position must match next node position
            r_next = nodes(1:3, k);
        else
            % Step #A.5: Terminal defect — propagated position must match DPO
            rv_dpo = cspice_spkez(target_id, et0+dT, frame, 'NONE', earth_id);
            r_next = rv_dpo(1:3);
        end

        ceq = [ceq; (rf_k - r_next) / position_scaling];
    end
end

function dxdt = nbody_shift_rhs(t, x, bodies, frame, earth_id, GM_Earth)
% AIM: Computes Earth-centred N-body acceleration using the Encke q-formulation.
%      This formulation minimizes numerical error for third-body perturbations.
% INPUTS:
%   t        [1x1]   SPICE ephemeris time [s].
%   x        [6x1]   State vector [r; v] in Earth-centred J2000 [km, km/s].
%   bodies   [cell]  Cell array of strings containing perturber names.
%   frame    [char]  Reference frame (e.g., 'J2000').
%   earth_id [1x1]   SPICE ID for Earth.
%   GM_Earth     [1x1]   Earth gravitational parameter [km^3/s^2].
% OUTPUT:
%   dxdt     [6x1]   State time derivative [vel; acc] [km/s, km/s^2].

    % Step #D.1.1: Calculate two-body acceleration (Earth gravity)
    rr   = x(1:3);  vv = x(4:6);
    dxdt = zeros(6,1);
    dxdt(1:3) = vv;
    dxdt(4:6) = -GM_Earth * rr / norm(rr)^3;    % Earth two-body term

    % Step #D.1.2: Retrieve Earth position relative to Solar System Barycentre (SSB)
    rvE  = cspice_spkez(earth_id, t, frame, 'NONE', 0);
    rr_E = rvE(1:3);
    rr_sc_ssb = rr + rr_E;                   % spacecraft wrt SSB
    % Step #D.1.3: Apply third-body perturbations via Encke q-formulation
    for j = 1:numel(bodies)
       % Compute third-body vector wrt Earth
        rr_b  = cspice_spkpos(bodies{j}, t, frame, 'NONE', 'SSB');
        rho_j = rr_b - rr_E;                 % body wrt Earth
        if norm(rho_j) < 1e-10, continue; end

        dd_j  = rr_sc_ssb - rr_b;            % spacecraft wrt body

        % Step #D.1.4: Compute numerical correction factors (q_i, f_i)
        q_i = rr.' * (rr - 2*rho_j) / (rho_j.' * rho_j);
        f_i = q_i * (3 + 3*q_i + q_i^2) / (1 + (1 + q_i)^1.5);
        % Step #D.1.5: Add perturbation acceleration
        GM_b = cspice_bodvcd(cspice_bodn2c(bodies{j}), 'GM', 1);
        dxdt(4:6) = dxdt(4:6) - GM_b(1) * (rr + rho_j*f_i) / norm(dd_j)^3;
    end
end

function [c, ceq] = MS_constraints(z, r0, v0, et0, frame, earth_id, bodies, target_id, GM_Earth,N)
% AIM: Formulates equality constraints for multiple-shooting trajectory optimization.
% INPUTS:
%   z        [28x1] Decision vector 
    %   z(1:3)=dv0, z(4:6)=r1,  z(7:9)=dv1,  z(10:12)=r2, z(13:15)=dv2,
    %   z(16:18)=r3, z(19:21)=dv3, z(22:24)=rf, z(25:27)=dv4, z(28)=dT
%   r0, v0   [3x1]  Initial spacecraft state [km, km/s].
%   et0      [1x1]  Initial mission epoch [s].
%   ... (Remaining are integration environment parameters).
% OUTPUT:
%   c        [0x0]  Inequality constraints (none).
%   ceq      [21x1] Equality constraints (segment defects + terminal matching).

    dv0 = z(1:3);   r1  = z(4:6);
    dv1 = z(7:9);   r2  = z(10:12);
    dv2 = z(13:15); r3  = z(16:18);
    dv3 = z(19:21); rf  = z(22:24);
    dv4 = z(25:27); dT  = z(28);
      
    etf_free         = et0 + dT;
    h                = dT / (N-1);
    ti               = et0 + (0:N-2)*h;        % N-1 node times
    %position_scaling = norm(r0);
   position_scaling = 1000;

    opts_c = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',h/200);
        
    % Collect waypoints and impulsive ΔVs into arrays (3 × N-1)
    r_nodes  = [r0,  r1,  r2,  r3 ];
    dv_nodes = [dv0, dv1, dv2, dv3];
    t_end    = [ti(2:end), etf_free];   % end time of each arc
    
    % Propagate the N-1 = 4 arcs
    xf = zeros(6, N-1);
    for k = 1 : N-1
        if k == 1
            v_init = v0 + dv_nodes(:,k);           % first arc: use given v0
        else
            v_init = xf(4:6, k-1) + dv_nodes(:,k); % subsequent: inherit velocity
        end
        [~, Xf] = ode113(@(t,x) nbody_shift_rhs(t,x,bodies,frame,earth_id,GM_Earth), ...
                         [ti(k), t_end(k)], [r_nodes(:,k); v_init], opts_c);
        xf(:,k) = Xf(end,:).';
    end
    
    % Target positions at the end of each arc
    rf_nodes = [r1, r2, r3, rf];               % 3 × N-1
    
    rv_dpo = cspice_spkez(target_id, etf_free, frame, 'NONE', earth_id);    
    % Build constraint vector
    c   = [];
    ceq = [];
    for k = 1 : N-1
        ceq = [ceq; (xf(1:3,k) - rf_nodes(:,k)) / position_scaling];
    end
    ceq = [ceq;
           (rf          - rv_dpo(1:3)) / position_scaling;   % position match at tf
           xf(4:6,end) + dv4 - rv_dpo(4:6)];                % velocity match at tf
end


