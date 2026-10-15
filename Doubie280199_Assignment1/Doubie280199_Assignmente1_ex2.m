%% Spacecraft Guidance and Navigation (2025/2026)
% Assignment 1 – Exercise 2: Low-Thrust Optimal Transfer around Eros
% Author : Karim Ahmed Mohamed Doubie
%
% Objective:
%   Point 1 – Debris density profile q(rho) and dimensional boundary states.
%   Point 2 – Non-dimensional scaling (DU, TU, VU, MU).
%   Point 3 – PMP formulation and zero-finding problem.
%   Point 4 – Planar (Delta_i = 0) optimal transfer via PMP single-shooting.
%   Point 5 – Inclined target (Delta_i = 3.5 deg) via Rx-rotation continuation.
%   Point 6 – N-body (SPICE) vs Kepler comparison over one month.

%% Step #0: Environment Initialisation and Configuration
clc; clearvars; close all;
format long;

% ── Global plot defaults ──────────────────────────────────────────────────
set(groot, 'DefaultAxesFontSize',                 16);
set(groot, 'DefaultAxesTitleFontSizeMultiplier',   1);
set(groot, 'DefaultTextFontSize',                 12);
set(groot, 'DefaultLineLineWidth',               1.5);
set(0,     'DefaultTextInterpreter',         'latex');
set(0,     'DefaultLegendInterpreter',       'latex');
set(0,     'DefaultAxesTickLabelInterpreter','latex');

% ── SPICE kernels ─────────────────────────────────────────────────────────
cspice_kclear();
cspice_furnsh('kernels/2000433.bsp');
cspice_furnsh('kernels/naif0012.tls');
cspice_furnsh('kernels/pck00010.tpc');
cspice_furnsh('kernels/de432s.bsp');
cspice_furnsh('kernels/gm_de432.tpc');

%% =======================================================================
%  Step #1 – Debris Density Profile and Boundary States
%% =======================================================================
fprintf('\n==================== Point 1 ====================\n');

% Step #1.1: Physical parameters (Table 4)
k1 = 7.393750e-3;             % [DU^-1]
k2 = 7.500000e-3;             % [DU^2]
k3 = 3.696875e-4;             % [DU^-1]
k4 = 6.250000e-4;             % [DU^2]
T_max = 21.59e-6;             % [N]
Isp   = 382.82;               % [s]
m0    = 25.0;                 % [kg]
g0    = 9.80665;              % [m/s^2]
R_a        = 17.00;           % [km]  Eros mean radius
mu_phys    = 3.5e-4;          % [km^3/s^2]
h_i        = 52.25;           % [km]  departure orbit altitude
h_f        = 36.15;           % [km]  arrival  orbit altitude
rho_dust_peak_A = 23.314 + R_a;    % [km]  first  dust-peak radius
rho_dust_peak_B = 42.170 + R_a;    % [km]  second dust-peak radius

% Step #1.2: Reference units and density profile  (DU = rho_i = h_i + R_a,  MU = m0)
MU = m0;                      % [kg]
DU = h_i + R_a;               % [km]  = 69.25 km ;               
q_debris = @(r_dist) k1 ./ (k2 + ((r_dist - rho_dust_peak_A)./DU).^2) + ...
                   k3 ./ (k4 + ((r_dist - rho_dust_peak_B)./DU).^2);

% Step #1.4: Locate density peaks numerically
rho_scan_A = linspace(20, 50, 10^4*7);
rho_scan_B = linspace(50, 70, 10^4*7);

[~, idx_A]  = max(q_debris(rho_scan_A));
[~, idx_B]  = max(q_debris(rho_scan_B));

rho_A_peak  = rho_scan_A(idx_A);
rho_B_peak  = rho_scan_B(idx_B);


fprintf('\n=== Peak Analysis ===\n');
met_i = ["Peak A [km]"; "Peak A [DU]"; "Peak B [km]"; "Peak B [DU]"];
val_i = [rho_A_peak; rho_A_peak/DU; rho_B_peak; rho_B_peak/DU];
disp(table(met_i, val_i, 'VariableNames', {'Metric','Value'}));

% Step #1.5: Visualise debris density profile
rho_values = linspace(20, 70, 3000);
figure();
hold on; grid on;

plot(rho_values, q_debris(rho_values), 'b-', 'LineWidth', 1.5, 'DisplayName', '$q(\rho)$');
xline(rho_A_peak, 'r--', 'LineWidth', 1.2, 'DisplayName', '$\rho_A$');
xline(rho_B_peak, 'k--', 'LineWidth', 1.2, 'DisplayName', '$\rho_B$');

% Annotate peaks with text arrows instead of numbers in legend
text(rho_A_peak + 0.5, max(q_debris(rho_values))*0.85, '$\rho_A$', ...
     'Interpreter','latex', 'FontSize',14, 'Color','r');
text(rho_B_peak + 0.5, max(q_debris(rho_values))*0.85, '$\rho_B$', ...
     'Interpreter','latex', 'FontSize',14, 'Color','k');

xlabel('$\rho\,[\mathrm{km}]$');
ylabel('$q(\rho)\,[\mathrm{DU}^{-1}]$');
title('Debris density profile $q(\rho)$');
legend('Location','best');

% Step #1.6: Dimensional boundary states
theta_i = pi/4;    % 45 deg
theta_f = 0;       %  0 deg

ri = R_a + h_i;   vi = sqrt(mu_phys/ri);
rf = R_a + h_f;   vf = sqrt(mu_phys/rf);

x_end = [ rf*cos(theta_f);  rf*sin(theta_f); 0; ...
          -vf*sin(theta_f);  vf*cos(theta_f); 0];
x_start = [ ri*cos(theta_i);  ri*sin(theta_i); 0; ...
          -vi*sin(theta_i);  vi*cos(theta_i); 0];


fprintf('\n=== State Vector Comparison ===\n');
met_i = ["x_start [km, km/s]"; "x_end [km, km/s]"];
val_i = [string(mat2str(x_start.', 12)); string(mat2str(x_end.', 12))];
disp(table(met_i, val_i, 'VariableNames', {'Metric','Value'}));


%% =======================================================================
%  Step #2 – Adimensionalisation
%% =======================================================================
fprintf('\n==================== Point 2 ====================\n');

% Step #2.1: Time and velocity units
TU = sqrt(DU^3 / mu_phys);            % [s]
VU = DU / TU;                          % [km/s]

% Step #2.2: Non-dimensional control parameters
mdot_max = T_max / (Isp * g0);        % [kg/s]
beta     = (mdot_max / MU) * TU;      % [-]

% T* = T_max * TU^2 / (MU * DU[m])
T_adim = T_max * TU^2 / (MU * (DU * 1e3));   % [-]

% Step #2.3: Non-dimensional peak radj
rho_A = rho_dust_peak_A / DU;
rho_B = rho_dust_peak_B / DU;

% Step #2.4: Non-dimensional states
x0_adim = [x_start(1:3)/DU; x_start(4:6)/VU; 1.0];   % [r*; v*; m*]
xf_adim = [x_end(1:3)/DU; x_end(4:6)/VU];          % [r*; v*]

fprintf('\n=== Normalisation Parameters ===\n');
met_i = ["DU [km]"; "TU [s]"; "VU [km/s]"; "T*"; "beta"; "rho_A*"; "rho_B*"];
val_i = [DU; TU; VU; T_adim; beta; rho_A; rho_B];
disp(table(met_i, val_i, 'VariableNames', {'Metric','Value'}));

% Step #2.5: Pack all parameters into the problem structure
prm.mu_adim   = 1.0;
prm.T_aim    = T_adim;
prm.beta      = beta;
prm.k1 = k1;  prm.k2 = k2;
prm.k3 = k3;  prm.k4 = k4;
prm.start_state  = x0_adim;
prm.end_state = xf_adim;
prm.rho_A     = rho_A;
prm.rho_B     = rho_B;


%% =======================================================================
%  Step #4 – Planar Optimal Transfer (Delta_i = 0)
%% =======================================================================
fprintf('\n==================== Point 4 ====================\n');

% Step #4.1: Solver configuration (Levenberg-Marquardt, tight tolerances)
fsolve_opt = optimoptions('fsolve', ...
    'Display',               'iter', ...
    'Algorithm',             'levenberg-marquardt', ...
    'FunctionTolerance',     1e-11, ...
    'StepTolerance',         1e-11, ...
    'MaxIterations',         500,   ...
    'MaxFunctionEvaluations', 10000);

% Steps #4.2–4.6: Initial guess, solve, propagate, evaluate
reiterate = false;
preconverged_guess = [-0.20445356732;% p_r0 x
                      -0.93475583137;% p_r0 y
                      0.00000000000; % p_r0 z
                      -0.01060250418;% p_v0 x
                      -0.64919708175;% p_v0 y
                      0.00000000000; % p_v0 z
                       2.25113531411;% p_m0
                       14.87498885648];% tf* [-]

[sol_planar, F_planar, solver_flag_planar, ...
 t_plan, Y_plan, metrics_plan] = ...
    solve_and_evaluate_arc(reiterate, preconverged_guess, prm, fsolve_opt, ...
                           TU, MU, VU,DU, m0, Isp, g0);

% Unpack for downstream use (Steps 5–6)
p_r0     = sol_planar(1:3);
p_v0     = sol_planar(4:6);
p_m0     = sol_planar(7);
tf_plan  = sol_planar(8);
err_r_km = metrics_plan.err_r_km;
err_v_ms = metrics_plan.err_v_ms;
mf_kg    = metrics_plan.final_mass_kg;
mprop_g  = (m0 - mf_kg) * 1e3;
dv_eq    = Isp * g0 * log(m0 / max(mf_kg, 1e-10));

% Step #4.7: Optimal planar trajectory plot
    % Custom Color Palette
    Sage       = [0.47, 0.53, 0.42];
    Terracotta = [0.80, 0.36, 0.27];
    AmberGold  = [0.93, 0.69, 0.13]; 
    InitialCol = [0.20, 0.20, 0.65]; 
    TargetCol  = [0.70, 0.15, 0.45]; 
    
    figure('Name', 'Planar Transfer 3D Analysis');
    hold on; grid on; axis equal;
    
    % 1. Render Eros as a "Gray Grid" body 
    [sx, sy, sz] = sphere(25);
    R_eros_plot = R_a / DU; 
    hEros = surf(sx*R_eros_plot, sy*R_eros_plot, sz*R_eros_plot);
    set(hEros, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', [0 0 0], ...
               'FaceAlpha', 1, 'LineWidth', 0.5, 'DisplayName', 'Eros Sphere');
    
    % 2. Plot the Initial and Final Parking Orbits (Thickened for evidence)
    % In the planar case, both orbits lie in the Z=0 plane
    th_ring = linspace(0, 2*pi, 1000);
    ri_star = norm(x0_adim(1:3));
    rf_star = norm(xf_adim(1:3));
    
    plot3(ri_star*cos(th_ring), ri_star*sin(th_ring), zeros(size(th_ring)), ...
        'Color', InitialCol, 'LineWidth', 1.5, 'DisplayName', 'Initial Orbit');
    
    plot3(rf_star*cos(th_ring), rf_star*sin(th_ring), zeros(size(th_ring)), ...
        'Color', TargetCol, 'LineWidth', 1.5, 'DisplayName', 'Target Orbit');
    
    % 3. Plot the Planar Trajectory in 3D 
    plot3(Y_plan(:,1), Y_plan(:,2), Y_plan(:,3), ...
        'Color', AmberGold, 'LineWidth', 3.0, 'DisplayName', 'Optimal Path');
    
    % 4. Plot Debris Peak Boundaries
    plot3((rho_A_peak/DU)*cos(th_ring), (rho_A_peak/DU)*sin(th_ring), zeros(size(th_ring)), ...
        'k--', 'LineWidth', 1.2, 'DisplayName', '$\rho_A$ Peak');
    plot3((rho_B_peak/DU)*cos(th_ring), (rho_B_peak/DU)*sin(th_ring), zeros(size(th_ring)), ...
        'k:',  'LineWidth', 1.2, 'DisplayName', '$\rho_B$ Peak');
    
    % 5. Plot Boundary State Markers
    plot3(x0_adim(1), x0_adim(2), x0_adim(3), 'o', ...
        'MarkerEdgeColor', Sage, 'MarkerFaceColor', Sage, 'MarkerSize', 10, 'DisplayName', 'Departure');
    plot3(xf_adim(1), xf_adim(2), xf_adim(3), 'p', ...
        'MarkerEdgeColor', Terracotta, 'MarkerFaceColor', Terracotta, 'MarkerSize', 13, ...
        'MarkerFaceColor', Terracotta, 'DisplayName', 'Arrival');
    
    % Aesthetics and Labels
    view(0, 90); %top-down for planar case
    xlabel('$x^*$ [DU]'); ylabel('$y^*$ [DU]'); zlabel('$z^*$ [DU]');
    title('3D View of the Planar Minimal-Exposure Transfer ($\Delta i = 0^\circ$)');
    legend('Location', 'northeastoutside', 'Interpreter', 'latex');

% Step #4.8: Diagnostic Analysis of the Planar Transfer Arc (Delta_i = 0)
    [S_plan, rho_plan, q_plan, H_plan] = arc_analysis(t_plan, Y_plan, prm);
    H_tf_plan = H_plan(end);
    t_days_plan = t_plan * TU / 86400;
    Sage       = [0.47, 0.53, 0.42];
    Terracotta = [0.80, 0.36, 0.27];
    Olive      = [0.50, 0.50, 0.00];
    figure('Name', 'Planar Optimal Arc Diagnostics');  
     % Group 1: Physical Parameters (Distance and Debris)
    subplot(1,2,1); hold on; grid on; 
    plot(t_days_plan, rho_plan * DU, 'LineWidth', 1.5);
    yline(rho_A_peak, '--', 'Color', [0.3 0.3 0.3]);
    yline(rho_B_peak, ':',  'Color', [0.3 0.3 0.3]);
    xlabel('$t$ [days]'); ylabel('$\rho$ [km]');
    title('Radial Distance Profile $\rho(t)$');

    subplot(1,2,2); hold on; grid on; 
    plot(t_days_plan, q_plan, 'Color', Olive, 'LineWidth', 1.5);
    xlabel('$t$ [days]'); ylabel('$q$ [$\mathrm{DU}^{-1}$]');
    title('Instantaneous Debris Exposure $q(\rho)$');
    % Group 2: Optimization Analysis (Hamiltonian and Switching)
    figure;
    subplot(1,2,1); hold on; grid on; 
    plot(t_days_plan, H_plan - H_tf_plan, 'Color', Terracotta, 'LineWidth', 1.5);
    xlabel('$t$ [days]'); ylabel('$H - H(t_f)$ [-]');
    title('Hamiltonian Stability $\Delta H$');

    subplot(1,2,2); hold on; grid on; 
    plot(t_days_plan, S_plan, 'Color', Sage, 'LineWidth', 1.5);
    yline(0, '-k', 'HandleVisibility', 'off');
    xlabel('$t$ [days]'); ylabel('$S$ [-]');
    title('Pontryagin Switching Logic $S(t)$');
    
    

%% =======================================================================
%  Step #5 – Inclined Target Orbit (Delta_i = 3.5 deg) via Continuation
%% =======================================================================
fprintf('\n==================== Point 5 ====================\n');

% Step #5.1: Pre-allocate result struct
inc_target_deg       = 3.5;
n_continuation_steps = 10;
cont_data = repmat(struct( ...
    ... % --- Convergence diagnostics ---
    'hamiltonian_dev',    [], ...
    'terminal_error',     [], ...
    ... % --- Performance metrics ---
    'delta_v_ms',         [], ...
    'propellant_mass_kg', [], ...
    'final_mass_kg',      [], ...
    ... % --- Optimisation solution ---
    'flight_time_days',   [], ...
    'costate_init',       [], ...
    ... % --- Continuation parameter ---
    'inclination_deg',    [], ...
    'transfer_angle_rad', [], ...
    ... % --- Full trajectory ---
    'state_history',      [], ...
    'time_history',       [] ), n_continuation_steps, 1);

costate_guess = sol_planar;    % warm-start from the planar solution
% Step #5.2: Continuation parameters
last_step_idx        = n_continuation_steps;
inc_grid_rad         = linspace(0, deg2rad(inc_target_deg), n_continuation_steps);

% Rotation matrix around the inertial x-axis
Rotation_x = @(theta_x) [1  0                 0;
                         0  cos(theta_x)  -sin(theta_x);
                         0  sin(theta_x)  cos(theta_x)];

xf_ref = xf_adim;    % store planar target for re-use at each step

% Step #5.3: Continuation solver options (looser display, same algorithm)
cont_opt = optimoptions('fsolve', ...
    'Display',               'off', ...
    'Algorithm',             'levenberg-marquardt', ...
    'FunctionTolerance',     1e-10, ...
    'StepTolerance',         1e-10, ...
    'MaxIterations',         400,   ...
    'MaxFunctionEvaluations', 10000);

% Step #5.4: Homotopy continuation loop over inclination
fprintf('Continuation 0 → %.1f deg in %d steps ...\n', inc_target_deg, n_continuation_steps);
wb = waitbar(0, 'Continuation: step 0 ...');
t_cont_start = tic;

for j = 1:n_continuation_steps
    inclination = inc_grid_rad(j);
    DCM_x       = Rotation_x(inclination);

    % Rotate both position and velocity of the target
    prm.end_state = [DCM_x * xf_ref(1:3); DCM_x * xf_ref(4:6)];

    [costate_opt, residual_bvp, solver_flag] = ...
        fsolve(@(guess) bvp_residual(guess, prm), costate_guess, cont_opt);

    if solver_flag > 0
        costate_guess = costate_opt;
    else
        warning('Step %d (%.3f deg) did not converge (exitflag=%d).', ...
                j, rad2deg(inclination), solver_flag);
    end

    % Step #5.5: Propagate and compute metrics for this step
    odeOpt  = odeset('RelTol',1e-11, 'AbsTol',1e-13);
    Y0_c = [prm.start_state; costate_opt(1:3); costate_opt(4:6); costate_opt(7)];
    [t_c, Y_c] = ode45(@(t,Y) eom_augmented(t,Y,prm), [0 costate_opt(8)], Y0_c, odeOpt);

    Yf_c = Y_c(end,:).';
    rf_c = Yf_c(1:3);  vf_c = Yf_c(4:6);  mf_c = Yf_c(7);

    h_orb           = cross(rf_c, vf_c);
    inc_c           = acos(h_orb(3)/norm(h_orb)) * 180/pi;
    tf_days_c       = costate_opt(8) * TU / 86400;
    mf_kg_c         = mf_c * MU;
    mprop_kg_c      = m0 - mf_kg_c;
    dv_ms_c         = Isp * g0 * log(m0 / max(mf_kg_c, 1e-10));
    [~, Hdev_c]     = check_H_constancy(t_c, Y_c, prm);
    err_c           = [rf_c - prm.end_state(1:3); vf_c - prm.end_state(4:6)];

    cont_data(j).hamiltonian_dev    = Hdev_c;
    cont_data(j).terminal_error     = err_c;
    cont_data(j).delta_v_ms         = dv_ms_c;
    cont_data(j).propellant_mass_kg = mprop_kg_c;
    cont_data(j).final_mass_kg      = mf_kg_c;
    cont_data(j).flight_time_days   = tf_days_c;
    cont_data(j).costate_init       = costate_opt;
    cont_data(j).inclination_deg    = inc_c;
    cont_data(j).transfer_angle_rad = inclination;
    cont_data(j).state_history      = Y_c;
    cont_data(j).time_history       = t_c;

    lbl_i = sprintf('Continuation Step %d/%d', j, n_continuation_steps);
    met_i = ["flag"; "||F||"; "inc [deg]"; "tf [days]"; "mf [kg]"; "dv_eq [m/s]"];
    val_i = [solver_flag; norm(residual_bvp); inc_c; tf_days_c; mf_kg_c; dv_ms_c];

    fprintf('\n=== %s ===\n', lbl_i);
    disp(table(met_i, val_i, 'VariableNames', {'Metric', 'Value'}));

    waitbar(j/n_continuation_steps, wb, ...
            sprintf('Continuation: step %d/%d (%.2f deg)', ...
                    j, n_continuation_steps, rad2deg(inclination)));
end
close(wb);
fprintf('Continuation CPU time: %.2f s\n', toc(t_cont_start));

% Step #5.6: Tabulate continuation results
fprintf('\n=== Continuation Summary ===\n');
val_matrix = zeros(n_continuation_steps, 5);
for j = 1:n_continuation_steps
    val_matrix(j,:) = [cont_data(j).inclination_deg,  cont_data(j).flight_time_days, ...
                       cont_data(j).delta_v_ms,        cont_data(j).final_mass_kg, ...
                       cont_data(j).hamiltonian_dev];
end
disp(array2table(val_matrix, ...
     'VariableNames', {'inc_deg','tf_days','dv_eq_ms','mf_kg','Hdev'}));

% Step #5.7: Extract final inclined solution
DCM_x_final   = Rotation_x(cont_data(last_step_idx).transfer_angle_rad);
prm.end_state = [DCM_x_final * xf_ref(1:3); DCM_x_final * xf_ref(4:6)];

Y_last      = cont_data(last_step_idx).state_history;
t_last      = cont_data(last_step_idx).time_history;
t_last_days = t_last * TU / 86400;
sol_inc     = cont_data(last_step_idx).costate_init;

p_r0_inc = sol_inc(1:3);
p_v0_inc = sol_inc(4:6);
p_m0_inc = sol_inc(7);
tf_inc   = sol_inc(8);

err_r_km_inc = norm(cont_data(last_step_idx).terminal_error(1:3)) * DU;
err_v_ms_inc = norm(cont_data(last_step_idx).terminal_error(4:6)) * VU * 1e3;


fprintf('\n=== Inclined Solution (Delta_i = %.1f deg) ===\n', inc_target_deg);
metrics = { ...
    'p_r0',             mat2str(p_r0_inc, 10);
    'p_v0',             mat2str(p_v0_inc, 10);
    'p_m0',             sprintf('%.10f', p_m0_inc);
    'tf* [-]',          sprintf('%.10f', tf_inc);
    'tf* [min]',        sprintf('%.10f', tf_inc * TU/60);
    '||dr_f|| [km]',    sprintf('%.10f', err_r_km_inc);
    '||dv_f|| [m/s]',   sprintf('%.10f', err_v_ms_inc);
    'mf [kg]',          sprintf('%.10f', cont_data(last_step_idx).final_mass_kg);
    'dv_eq [m/s]',      sprintf('%.10f', cont_data(last_step_idx).delta_v_ms);
    'max|H-Hf|',        sprintf('%.10e', cont_data(last_step_idx).hamiltonian_dev)};

disp(cell2table(metrics, 'VariableNames', {'Metric', 'Value'}));


% Step #5.8: 3D Visualization of the Optimal Inclined Transfer
    Sage       = [0.47, 0.53, 0.42];
    Terracotta = [0.80, 0.36, 0.27];
    VividTeal  = [0.00, 0.75, 0.75]; 
    InitialCol = [0.20, 0.20, 0.65]; 
    TargetCol  = [0.70, 0.15, 0.45]; 
    figure('Name', '3D Optimal Arc & Orbit Comparison');
    hold on; grid on; axis equal;
    % 1. Render Eros as a "Gray Grid" body 
    [sx, sy, sz] = sphere(25);
    R_eros_plot = R_a / DU; 
    hEros = surf(sx*R_eros_plot, sy*R_eros_plot, sz*R_eros_plot);
    set(hEros, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', [0 0 0], ...
               'FaceAlpha', 1, 'LineWidth', 0.5, 'DisplayName', 'Eros');
    % 2. Plot the Initial and Final Parking Orbits 
    % Initial Orbit (Reference plane)
    ri_star = norm(x0_adim(1:3));
    plot3(ri_star*cos(th_ring), ri_star*sin(th_ring), zeros(size(th_ring)), ...
        'Color', InitialCol, 'LineWidth', 1.5, 'DisplayName', 'Initial Orbit');
    % Final Orbit (Inclined by 3.5 deg)
    rf_star = norm(prm.end_state(1:3));
    DCM_inc = [1 0 0; 0 cos(deg2rad(3.5)) -sin(deg2rad(3.5)); 0 sin(deg2rad(3.5)) cos(deg2rad(3.5))];
    ring_f = DCM_inc * [rf_star*cos(th_ring); rf_star*sin(th_ring); zeros(size(th_ring))];
    plot3(ring_f(1,:), ring_f(2,:), ring_f(3,:), ...
        'Color', TargetCol, 'LineWidth', 1.5, 'DisplayName', 'Target Orbit');
    % 3. Plot the 3D Trajectory 
    plot3(Y_last(:,1), Y_last(:,2), Y_last(:,3), ...
        'Color', VividTeal, 'LineWidth', 3.0, 'DisplayName', 'Optimal Path');
    % 4. Plot Debris Peak Boundaries
    plot3((rho_A_peak/DU)*cos(th_ring), (rho_A_peak/DU)*sin(th_ring), zeros(size(th_ring)), ...
        'k--', 'LineWidth', 1.2, 'DisplayName', '$\rho_A$ Peak');
    plot3((rho_B_peak/DU)*cos(th_ring), (rho_B_peak/DU)*sin(th_ring), zeros(size(th_ring)), ...
        'k:',  'LineWidth', 1.2, 'DisplayName', '$\rho_B$ Peak');
    % 5. Plot Boundary State Markers 
    plot3(x0_adim(1), x0_adim(2), x0_adim(3), 'o', ...
        'MarkerEdgeColor', Sage, 'MarkerFaceColor', Sage, 'MarkerSize', 10, 'DisplayName', 'Departure');
    plot3(prm.end_state(1), prm.end_state(2), prm.end_state(3), 'p', ...
        'MarkerEdgeColor', Terracotta, 'MarkerFaceColor', Terracotta, 'MarkerSize', 13, ...
        'MarkerFaceColor', Terracotta, 'DisplayName', 'Arrival');
    % Aesthetics and Labels
    view(50, 25); 
    xlabel('$x^*$ [DU]'); ylabel('$y^*$ [DU]'); zlabel('$z^*$ [DU]');
    title('3D Low-Thrust Transfer: Initial vs. Target Plane ($\Delta i = 3.5^\circ$)');
    legend('Location', 'northeastoutside', 'Interpreter', 'latex');
% Step #5.9: Yz Plane projection Visualization of the Optimal Inclined Transfer
    figure('Name', 'YZ Plane Projection');
    hold on; grid on; axis equal;
    % 1. Eros cross-section on YZ plane
    th_ring_yz = linspace(0, 2*pi, 300);
    fill(R_eros_plot*cos(th_ring_yz), R_eros_plot*sin(th_ring_yz), ...
         [0.4 0.4 0.4], 'EdgeColor','k', 'LineWidth',0.5, 'DisplayName','Eros');
    % 2. Initial orbit projection onto YZ
    plot(ri_star*sin(th_ring), zeros(size(th_ring)), ...
         'Color', InitialCol, 'LineWidth', 1.5, 'DisplayName', 'Initial Orbit');
    % 3. Target orbit projection onto YZ (inclined)
    plot(ring_f(2,:), ring_f(3,:), ...
         'Color', TargetCol, 'LineWidth', 1.5, 'DisplayName', 'Target Orbit');
    % 4. Optimal trajectory projected onto YZ
    plot(Y_last(:,2), Y_last(:,3), ...
         'Color', VividTeal, 'LineWidth', 3.0, 'DisplayName', 'Optimal Path');
    % 5. Debris boundaries on YZ
    plot((rho_A_peak/DU)*sin(th_ring), zeros(size(th_ring)), ...
         'k--', 'LineWidth', 1.2, 'DisplayName', '$\rho_A$ Peak');
    plot((rho_B_peak/DU)*sin(th_ring), zeros(size(th_ring)), ...
         'k:',  'LineWidth', 1.2, 'DisplayName', '$\rho_B$ Peak');
    % 6. Departure and Arrival markers
    plot(x0_adim(2), x0_adim(3), 'o', ...
         'MarkerEdgeColor', Sage, 'MarkerFaceColor', Sage, ...
         'MarkerSize', 10, 'DisplayName', 'Departure');
    plot(prm.end_state(2), prm.end_state(3), 'p', ...
         'MarkerEdgeColor', Terracotta, 'MarkerFaceColor', Terracotta, ...
         'MarkerSize', 13, 'DisplayName', 'Arrival');
    xlabel('$y^*$ [DU]'); ylabel('$z^*$ [DU]');
    title('YZ Plane Projection: Plane Change $\Delta i = 3.5^\circ$');
    legend('Location','northeastoutside','Interpreter','latex');

            
% Step #5.9: Planar vs inclined comparison
    [~, Hdev_plan] = check_H_constancy(t_plan, Y_plan, prm);
    fprintf('\n=== Comparison: Planar vs Inclined ===\n');
    met_i = ["Planar: tf*"; "Planar: max|H-Hf|"; "Planar: ||err||"; ...
             "Inclined: tf*"; "Inclined: max|H-Hf|"; "Inclined: ||err||"; ...
             "Inclined: inc [deg]"; "Inclined: dv_eq [m/s]"];
    val_i = [tf_plan; Hdev_plan; norm([err_r_km/DU; err_v_ms/(VU*1e3)]); ...
             tf_inc; cont_data(last_step_idx).hamiltonian_dev; ...
             norm(cont_data(last_step_idx).terminal_error); ...
             cont_data(last_step_idx).inclination_deg; ...
             cont_data(last_step_idx).delta_v_ms];
    disp(table(met_i, val_i, 'VariableNames', {'Metric','Value'}));


% Step #5.11: Inclined arc diagnostic plots
    [S_inc, rho_inc, q_inc, H_inc] = arc_analysis(t_last, Y_last, prm);
    H_tf_inc = H_inc(end);
    Sage = [0.47, 0.53, 0.42]; Terracotta = [0.80, 0.36, 0.27];
    Slate = [0.27, 0.51, 0.71]; Olive = [0.50, 0.50, 0.00];
    
    % Group 1: Physical Parameters (Distance and Debris)
    figure('Name', 'Orbit Evolution (Inclined)');
    subplot(1,2,1); hold on; grid on;
    plot(t_last_days, rho_inc*DU, 'Color', Slate, 'LineWidth', 1.5);
    yline(rho_A_peak, '--', 'Color', [0.3 0.3 0.3]); 
    yline(rho_B_peak, ':',  'Color', [0.3 0.3 0.3]);
    xlabel('$t$ [days]'); ylabel('$\rho$ [km]');
    title('Radial Distance Profile $\rho(t)$ ($\Delta i = 3.5^\circ$)');
    
    subplot(1,2,2); hold on; grid on;
    plot(t_last_days, q_inc, 'Color', Olive, 'LineWidth', 1.5);
    xlabel('$t$ [days]'); ylabel('$q$ [$\mathrm{DU}^{-1}$]');
    title('Instantaneous Debris Exposure $q(\rho)$');
    
    % Group 2: Optimization Diagnostics (Hamiltonian and Switching)
    figure('Name', 'Control Diagnostics (Inclined)');
    subplot(1,2,1); hold on; grid on;
    plot(t_last_days, H_inc - H_tf_inc, 'Color', Terracotta, 'LineWidth', 1.5);
    xlabel('$t$ [days]'); ylabel('$H - H(t_f)$ [-]');
    title('Hamiltonian Stability $\Delta H$ ($\Delta i = 3.5^\circ$)');
    
    subplot(1,2,2); hold on; grid on;
    plot(t_last_days, S_inc, 'Color', Sage, 'LineWidth', 1.5);
    yline(0, '-k', 'HandleVisibility', 'off');
    xlabel('$t$ [days]'); ylabel('$S$ [-]');
    title('Pontryagin Switching Logic $S(t)$');


%% =======================================================================
%  Step #6 – N-body vs Kepler Propagation (one-month check)
%% =======================================================================
fprintf('\n==================== Point 6 ====================\n');

% Step #6.1: Preliminary Manipulation SPICE
eros_id = cspice_bodn2c('EROS');
mu_eros = cspice_bodvcd(eros_id,'GM',1);
mu_eros = mu_eros(1);
et0 = cspice_str2et('2012-01-15T00:00:00');
etf = cspice_str2et('2012-02-15T00:00:00');

% Initial condition: dimensional final state from the inclined arc
Yf_inc = Y_last(end,:).';
           %r0_nb          %v0_nb
x0_nb  = [Yf_inc(1:3)*DU; Yf_inc(4:6)*VU];   % [km; km/s]

% Step #6.2: Third-body perturber data
body_labels = {'SUN', ...
               'EARTH','MOON','MARS BARYCENTER','VENUS', ...
               'URANUS BARYCENTER','NEPTUNE BARYCENTER','MERCURY', ...
               'PLUTO BARYCENTER','JUPITER BARYCENTER','SATURN BARYCENTER'};
bodies = cell(numel(body_labels),1);
for j = 1:numel(body_labels)
    id = cspice_bodn2c(body_labels{j});
    gm = cspice_bodvcd(id,'GM',1);
    bodies{j}.GM   = gm(1);
    bodies{j}.name = body_labels{j};
end

% Step #6.3–6.5: Run comparison via helper
frame      = 'ECLIPJ2000';
integ_opts = odeset('RelTol',1e-10, ...
                    'AbsTol',[1e-6*ones(3,1); 1e-9*ones(3,1)], ...
                    'MaxStep',(etf-et0)/1500);

[t_grid, traj_nbody, traj_kepler, r_error_norm, dR, dT, dN] = ...
    propagate_and_compare(x0_nb, et0, etf, mu_eros, bodies, frame, eros_id, integ_opts);

% Step #6.6: Error statistics
days = (t_grid - et0) / 86400;
fprintf('\n=== Position Error Analysis ===\n');
met_i = ["Final ||r_nb-r_kep|| [km]"; "Max ||r_nb-r_kep|| [km]"; "RMS ||r_nb-r_kep|| [km]"];
val_i = [r_error_norm(end); max(r_error_norm); sqrt(mean(r_error_norm.^2))];
disp(table(met_i, val_i, 'VariableNames', {'Metric','Value'}));

fprintf('\n=== RTN Error Statistics [km] ===\n');
met_i = ["Max R";"Max T";"Max N";"RMS R";"RMS T";"RMS N"];
val_i = [max(abs(dR)); max(abs(dT)); max(abs(dN)); ...
         sqrt(mean(dR.^2)); sqrt(mean(dT.^2)); sqrt(mean(dN.^2))];
disp(table(met_i, val_i, 'VariableNames', {'Metric','Value'}));


% Step #6.7: Plots
% --- Figure 1: 3D Trajectory Comparison (Separated) ---
    figure('Name', 'Validation: N-Body vs Keplerian Dynamics'); 
    hold on; grid on; axis equal;
    % 1. Render Eros as a "Gray Grid" physical body 
    [sx, sy, sz] = sphere(25);
    hEros = surf(sx*R_a, sy*R_a, sz*R_a);
    set(hEros, ...
        'FaceColor', [0.4 0.4 0.4], ... 
        'EdgeColor', [0 0 0],       ... 
        'FaceAlpha', 1,             ... 
        'LineWidth', 0.5);
    plot3(0, 0, 0, 'go', 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % 2. Plot Trajectories (N-body in Red, Kepler in Black dashed)
    plot3(traj_nbody(:,1), traj_nbody(:,2), traj_nbody(:,3), 'r', 'LineWidth', 1.5, 'DisplayName', 'Perturbed (N-Body)');
    plot3(traj_kepler(:,1), traj_kepler(:,2), traj_kepler(:,3), 'k--', 'LineWidth', 1.2, 'DisplayName', 'Reference (Kepler)');
    % Aesthetics and Labels
    xlabel('$x_{\mathrm{Eros}}$ [km]');
    ylabel('$y_{\mathrm{Eros}}$ [km]');
    zlabel('$z_{\mathrm{Eros}}$ [km]');
    title('3D Comparison: Perturbed N-Body Trajectory vs. Ideal Keplerian');
    % Update legend to include the Body
    legend('Eros Body', 'Perturbed (N-Body)', 'Reference (Kepler)', 'Location', 'best');
    view(40, 25);

% --- Figure 2: Comprehensive Error Analysis (Unique Subplots) ---
    figure('Name', 'Detailed Displacement Analysis');
    % Subplot 1: Magnitude of the position error
    subplot(2,1,1); hold on; grid on;
    plot(days, r_error_norm, 'm', 'LineWidth', 1.5); 
    xlim([0 31]); 
    xlabel('$t-t_0$ [days]');
    ylabel('$\|\Delta \mathbf{r}\|$ [km]');
    title('Absolute Divergence over Time');
    
    % Subplot 2: RTN decomposition
    subplot(2,1,2); hold on; grid on;
    plot(days, dR, 'c', 'LineWidth', 1.2); 
    plot(days, dT, 'g', 'LineWidth', 1.2); 
    plot(days, dN, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.2); 
    xlim([0 31]); 
    xlabel('$t-t_0$ [days]');
    ylabel('RTN Components [km]');
    legend('Radial ($dR$)', 'Transverse ($dT$)', 'Normal ($dN$)', 'Location', 'best');
    title('Relative Displacement Projected in the Local RTN Frame');


%% ====================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================
function dY = eom_augmented(~, Y, prm)
% AIM: Compute the time derivative of the augmented state-costate vector
%      based on Pontryagin's Minimum Principle (PMP) for optimal trajectory control.
% INPUTS:
%   ~    : [1x1] Time (unused for autonomous system).
%   Y    : [14x1] Augmented state[r(3); v(3); m; p_r(3); p_v(3); p_m].
%   prm  : [struct] Problem parameters (thrust, gravity, running cost weights).
% OUTPUT:
%   dY   : [14x1] Derivative of the augmented state[r_dot; v_dot; m_dot; p_r_dot; p_v_dot; p_m_dot].

    % Step #1: Unpack augmented state and costate vectors
    r   = Y(1:3);  v   = Y(4:6);  m   = Y(7);
    p_r = Y(8:10); p_v = Y(11:13); p_m = Y(14);

    % Step #2: Pre-compute radial norms and scalar components
    rn   = norm(r);
    r3   = rn^3;
    r5   = rn^5;
    pv_n = norm(p_v);

    % Step #3: Define optimal control (Thrust direction and Bang-bang throttle)
    if pv_n > 0
        u_hat = -p_v / pv_n;
    else
        u_hat = zeros(3,1);
    end

    % Switching function and bang-bang throttle
    S_sw = -(prm.T_aim / max(m,1e-10)) * pv_n - p_m * prm.beta;
    u    = 0.5*(1 - sign(S_sw));

    % Step #4: Compute state equations (Kinematics and Dynamics)
    grav_acc = -(prm.mu_adim / r3) * r;
    thrust_acc = (u * prm.T_aim / max(m, 1e-10)) * u_hat;
    vdot = grav_acc + thrust_acc;
    rdot = v;
    mdot = -u * prm.beta;

    % Step #5: Compute running-cost gradient via chain rule (dq/dr)
    slope_q = -2*prm.k1*(rn - prm.rho_A) / (prm.k2 + (rn - prm.rho_A)^2)^2 ...
           -2*prm.k3*(rn - prm.rho_B) / (prm.k4 + (rn - prm.rho_B)^2)^2;
    grad_q_r = slope_q * (r / max(rn, 1e-10));

    % Step #6: Compute costate equations (Adjoint Dynamics)
    p_r_dot = (prm.mu_adim/r3)*p_v - (3*prm.mu_adim/r5)*dot(r,p_v)*r - grad_q_r;

    p_v_dot = -p_r;
    p_m_dot = -u * (prm.T_aim / max(m,1e-10)^2) * pv_n;
    % Step #7: Assemble final augmented derivative vector
    dY = [rdot; vdot; mdot; p_r_dot; p_v_dot; p_m_dot];
end

% -------------------------------------------------------------------------
function H = eval_H(Y, prm)
% AIM: Evaluates the scalar Hamiltonian H of the system at a given augmented state.
% INPUTS:
%   Y    : [14x1] Augmented state-costate vector[pos(3); vel(3); m; p_r(3); p_v(3); p_m].
%   prm  :[struct] System parameters (mu, thrust constants, weights).
% OUTPUT:
%   H    : [1x1] Scalar Hamiltonian value.

    % Step #1: Unpack state (r, v, m) and costates (p_r, p_v, p_m)
    p_r = Y(8:10); p_v = Y(11:13); p_m = Y(14);
    r   = Y(1:3);  v   = Y(4:6);  m   = Y(7);
    % Step #2: Compute geometric norms and optimal thrust direction
    rn   = norm(r);
    r3   = rn^3;
    pv_n = norm(p_v);

    if pv_n > 0
        u_hat = -p_v / pv_n;
    else
        u_hat = zeros(3,1);
    end
    % Step #3: Calculate switching function and throttle (bang-bang control)
    S_sw = - p_m * prm.beta -(prm.T_aim / max(m,1e-10)) * pv_n ;
    u    = 0.5*(1 - sign(S_sw));
   % Step #4: Assemble the Hamiltonian scalar
    H = prm.k1/(prm.k2 + (rn - prm.rho_A)^2) + prm.k3/(prm.k4 + (rn - prm.rho_B)^2) ...
      + p_r(1)*v(1) + p_r(2)*v(2) + p_r(3)*v(3) ...
      - (prm.mu_adim / r3) * dot(p_v, r) ...
      + u * (prm.T_aim / max(m,1e-10)) * dot(p_v, u_hat) ...
      - u * prm.beta * p_m;
end

% -------------------------------------------------------------------------
function [S_vec, rho_vec, q_vec, H_vec] = arc_analysis(t_vec, Y_mat, prm)
% AIM: Extract and compute diagnostic time-series along the optimal arc 
%      to verify control switching, distance profiles, and Hamiltonian integrity.
% INPUTS:
%   t_vec  : [Nx1] Time vector [adimensional].
%   Y_mat  : [Nx14] Augmented state-costate history matrix.
%   prm    : [struct] System parameters (thrust, gravity, weights).
% OUTPUTS:
%   S_vec   : [Nx1] Switching function history [-].
%   rho_vec : [Nx1] Radial distance history [adimensional].
%   q_vec   : [Nx1] Running cost function history [-].
%   H_vec   : [Nx1] Hamiltonian history [-].

    % Step #1: Extract state and costate components from the augmented history
    pv_mat  = Y_mat(:,11:13);                  % [Nx3]
    pm_vec  = Y_mat(:,14);                     % [Nx1]
    r_mat   = Y_mat(:,1:3);                    % [Nx3]
    m_vec   = Y_mat(:,7);                      % [Nx1]

    % Step #2: Compute radial distance and velocity costate norms
    rho_vec = vecnorm(r_mat, 2, 2);            % ||r|| at each step
    pv_n_vec = vecnorm(pv_mat, 2, 2);          % ||p_v|| at each step

    % Step #3: Calculate the switching function S(t) (vectorised)
    S_vec = -pm_vec.*prm.beta - (prm.T_aim .* pv_n_vec) ./ max(m_vec, 1e-10);

    % Step #4: Calculate the running cost q(rho(t)) (vectorised)
    q_vec  = prm.k1 ./ (prm.k2 + (rho_vec - prm.rho_A).^2) + prm.k3 ./ (prm.k4 + (rho_vec - prm.rho_B).^2);

    % Step #5: Compute Hamiltonian trajectory (scalar evaluation loop)
    n     = length(t_vec);
    H_vec = zeros(n,1);
    for step = 1:n
        H_vec(step) = eval_H(Y_mat(step,:).', prm);
    end
end

% -------------------------------------------------------------------------
function [final_ham, max_ham_drift] = check_H_constancy(t, Y, prm)
% AIM: Validate PMP optimality by computing Hamiltonian constancy over the arc.
% INPUTS:
%   t             : [Nx1] Time vector[adimensional].
%   Y             : [Nx14] State and costate history matrix.
%   prm           : [struct] Problem parameters (thrust, gravity, weights).
% OUTPUTS:
%   final_ham     : [1x1] Terminal Hamiltonian value (should be 0 for free tf).
%   max_ham_drift : [1x1] Maximum absolute deviation of H(t) from H(tf).

  % Step #1: Memory pre-allocation for efficiency
    n = numel(t);
    H = zeros(n,1);
  % Step #2: Compute Hamiltonian at each discrete time step
    for idx = 1:n
        % Pass individual row as column vector for internal eval_H logic
        H(idx) = eval_H(Y(idx,:).', prm);
    end
  % Step #3: Extract boundary value and compute absolute divergence
    final_ham  = H(end);
    max_ham_drift = max(abs(H - final_ham));
end

% -------------------------------------------------------------------------
function sdot = nbody_rhs_eros(t, x, mu_eros, bodies, frame, eros_id)
% AIM: Compute spacecraft EOM relative to Eros including N-body perturbations.
% INPUTS:
%   t       : [1x1] Current time (ET) [s].
%   x       : [6x1] State vector [pos; vel] in Eros-centered frame [km, km/s].
%   mu_eros : [1x1] Eros gravitational parameter [km^3/s^2].
%   bodies  : [cell] Struct array of third-body perturbers (.name, .GM).
%   frame   : [char] Reference inertial frame string.
%   eros_id : [1x1] SPICE ID for the central body (Eros).
% OUTPUTS:
%   sdot    : [6x1] State derivative[vel; acc] [km/s, km/s^2].

    % Step #1: Initialize velocity component of the derivative
    sdot = zeros(6,1);
    % Step #2: Retrieve central body (Eros) inertial position
    sdot(1:3) = x(4:6);

    % Step #3: Calculate central two-body gravitational acceleration
    rE    = cspice_spkez(eros_id, t, frame, 'NONE', 0);
    rE    = rE(1:3);
    % Two-body central acceleration (Eros)
    a_acc = -(mu_eros /  norm(x(1:3))^3) * x(1:3);

    % Step #4: Accumulate third-body perturbation accelerations
    for i = 1:numel(bodies)
        rB_ssb = cspice_spkezr(bodies{i}.name, t, frame, 'NONE', 'SSB');
        a_acc = a_acc + bodies{i}.GM * ...
                (  (rB_ssb(1:3) - (rE + x(1:3))) /norm(rB_ssb(1:3) - (rE + x(1:3)))^3 - (rB_ssb(1:3) - rE)/norm(rB_ssb(1:3) - rE)^3);
    end

    % Step #5: Map total acceleration to state derivative vector
    sdot(4:6) = a_acc;
end



% -------------------------------------------------------------------------
function [t_grid, traj_nb, traj_kep, r_err, dR, dT, dN] = ...
         propagate_and_compare(x0, t0, tf, mu_cb, bodies, frame, cb_id, opts)
% AIM: Performs comparative propagation between N-body and Keplerian dynamics, 
%      then decomposes the resulting position error into the Radial-Transverse-Normal (RTN) frame.
% INPUTS:
%   x0        :[6x1]   Initial state [pos (km); vel (km/s)].
%   t0, tf    : [1x1]   Initial and final SPICE ephemeris times [s].
%   mu_cb     : [1x1]   Central body gravitational parameter [km^3/s^2].
%   bodies    : [cell]  Array of structs containing third-body perturber data (.name, .GM).
%   frame     : [char]  Reference inertial frame string.
%   cb_id     :[1x1]   SPICE ID of the central body.
%   opts      : [struct] ODE solver options.
% OUTPUTS:
%   t_grid    : [Nx1]   Time vector [s].
%   traj_nb   : [Nx6]   N-body state history [km, km/s].
%   traj_kep  : [Nx6]   Keplerian state history [km, km/s].
%   r_err     : [Nx1]   Magnitude of position error [km].
%   dR, dT, dN: [Nx1]   RTN components of position error [km].


    % Step #1: Perform N-body propagation
    [t_grid, traj_nb] = ode113( ...
        @(t,x) nbody_rhs_eros(t,x,mu_cb,bodies,frame,cb_id), ...
        [t0 tf], x0, opts);

    % Step #2: Perform Keplerian propagation on the same time grid
    rhs_kep = @(t,x) [x(4:6); -mu_cb * x(1:3) / norm(x(1:3))^3];
    [~, traj_kep] = ode113(rhs_kep, t_grid, x0, opts);

    % Step #3: Calculate magnitude of the position error vector
    delta_r = traj_nb(:,1:3) - traj_kep(:,1:3);
    r_err   = vecnorm(delta_r, 2, 2);

    % Step #4: Perform RTN frame decomposition
    n  = numel(t_grid);
    dR = zeros(n,1);  
    dT = dR; 
    dN = dR;
    for idx = 1:n
        r_ref = traj_kep(idx,1:3).';
        v_ref = traj_kep(idx,4:6).';
        uR = r_ref / norm(r_ref);                          % Radial
        uN = cross(r_ref,v_ref) / norm(cross(r_ref,v_ref));% Normal (angular momentum)
        uT = cross(uN, uR);                                % Transverse
        % Project error vector onto basis
        dR(idx) = dot(delta_r(idx,:).', uR);
        dT(idx) = dot(delta_r(idx,:).', uT);
        dN(idx) = dot(delta_r(idx,:).', uN);
    end
end
% -------------------------------------------------------------------------

function[sol, F, flag, t_arc, Y_arc, metrics] = ...
         solve_and_evaluate_arc(reiterate, preconverged_guess, prm, fsolve_opt, ...
                                TU, MU, VU, DU, m0, Isp, g0)
% AIM: Coordinates the initialization, shooting-based optimization, 
%      propagation, and performance evaluation of an optimal transfer arc.
% INPUTS:
%   reiterate          : [logical] True to perform random coarse search for guess.
%   preconverged_guess : [8x1]     Initial guess vector[lambda_r0(3); lambda_v0(3); lambda_m0; tf*].
%   prm                : [struct]  Problem parameters and constants.
%   fsolve_opt         : [struct]  Optimization algorithm settings.
%   TU, MU, VU, DU     : [1x1]     Dimensional normalization constants.
%   m0, Isp, g0        :[1x1]     Propulsion and physical constants.
% OUTPUTS:
%   sol   : [8x1] Converged decision vector.
%   F     : [8x1] Residual vector at solution.
%   flag  : [1x1] Solver exitflag.
%   t_arc : [Nx1] Time vector [adimensional].
%   Y_arc : [Nx14] State-costate history.
%   metrics: [struct] Computed terminal performance indicators.

    % Step #1: Initialization of optimization guess
    if reiterate
        sol0 = initialise_costate_guess_(prm);
    else
        sol0 = preconverged_guess;
    end

    % Step #2: Solve the TPBVP via single-shooting method
    [sol, F, flag] = fsolve(@(guess) bvp_residual(guess, prm), sol0, fsolve_opt);

    if flag <= 0
        warning('Shooting did not converge (exitflag = %d).', flag);
    end

    % Step #3: Unpack and report converged costates
    metrics.p_r0    = sol(1:3);
    metrics.p_v0    = sol(4:6);
    metrics.p_m0    = sol(7);
    metrics.tf_adim = sol(8);
    metrics.tf_min  = sol(8) * TU / 60;

    fprintf('\n=== Shooting Solution ===\n');
    met_i =["p_r0"; "p_v0"; "p_m0"; "tf* [-]"; "tf* [min]"; "||F||"];
    val_i =[string(mat2str(metrics.p_r0.',10)); ...
             string(mat2str(metrics.p_v0.',10)); ...
             metrics.p_m0; metrics.tf_adim; metrics.tf_min; norm(F)];
    disp(table(met_i, val_i, 'VariableNames', {'Metric','Value'}));
    


    % Step #4: Propagate the converged optimal arc
    ode_opt = odeset('RelTol', 1e-11, 'AbsTol', 1e-13);
    Y0_arc  =[prm.start_state; metrics.p_r0; metrics.p_v0; metrics.p_m0];
    [t_arc, Y_arc] = ode45(@(t,Y) eom_augmented(t,Y,prm), ...
                           [0, metrics.tf_adim], Y0_arc, ode_opt);

    % Step #5: Compute terminal and physical performance metrics
    Yf = Y_arc(end,:).';
    metrics.err_r_km        = norm(Yf(1:3) - prm.end_state(1:3)) * DU;
    metrics.err_v_ms        = norm(Yf(4:6) - prm.end_state(4:6)) * VU * 1e3;
    metrics.final_mass_kg   = Yf(7) * MU;
    metrics.propellant_mass_g = (m0 - metrics.final_mass_kg) * 1e3;
    metrics.delta_v_ms      = Isp * g0 * log(m0 / max(metrics.final_mass_kg, 1e-10));

    fprintf('\n=== Terminal and Physical Metrics ===\n');
    met_i = ["||dr_f|| [km]"; "||dv_f|| [m/s]"; ...
             "mf [kg]"; "m_prop [g]"; "dv_eq [m/s]"];
    val_i =[metrics.err_r_km;  metrics.err_v_ms; metrics.final_mass_kg; ...
             metrics.propellant_mass_g; metrics.delta_v_ms];
    disp(table(met_i, val_i, 'VariableNames', {'Metric','Value'}));
end

% -----------------------------------------------------------------------
function sol0 = initialise_costate_guess_(prm)
% AIM: Generate near-optimal initial costates and flight time using a 
%      Monte Carlo (random shooting) approach to warm-start the TPBVP solver.
% INPUTS:
%   prm  : [struct] Problem parameters (dynamics, start/end states).
% OUTPUTS:
%   sol0 : [8x1] Vector[lambda_r0(3); lambda_v0(3); lambda_m0; tf*].

    % Step #1: Configuration and workspace setup
    % rng(42);
    num_samples    = 5*10000;
    r_i_adim       = norm(prm.start_state(1:3));
    T_orb_adim     = 2*pi*sqrt(r_i_adim^3);
    min_residual   = 1e20;
    sol0           = zeros(8,1);
    opts_cs        = odeset('RelTol',1e-8,'AbsTol',1e-10);

    fprintf('Coarse random search (%d candidates, seed 42) ...\n', num_samples);
    wb_cs = waitbar(0,'Coarse search: candidate 0 ...');

    % Step #2: Iterative random shooting loop
    for j = 1:num_samples
       
        % Step #2.1: Sample random costate space and flight time
        p_r_cand = -100 + 200*rand(3,1);
        p_v_cand = -100 + 200*rand(3,1);
        p_m_cand = -100 + 100*rand;

        tf_cand  = (1 + 1.5*rand) * T_orb_adim;
        Y0_cand  =[prm.start_state; p_r_cand; p_v_cand; p_m_cand];

        % Step #2.2: Integrate augmented dynamics using current candidates
        [~, Y] = ode45(@(t,Y) eom_augmented(t,Y,prm), [0 tf_cand], Y0_cand, opts_cs);
        Yf = Y(end,:).';
        
        % Step #2.3: Evaluate candidate performance and update best result
        if ~any(isnan(Yf)) && all(isfinite(Yf))
            e_cand = norm([Yf(1:3)-prm.end_state(1:3); ...
                           Yf(4:6)-prm.end_state(4:6)]);
            if e_cand < min_residual
                min_residual = e_cand;
                sol0 =[p_r_cand; p_v_cand; p_m_cand; tf_cand];
            end
        end

        % Step #2.4: Update user progress bar
        waitbar(j/num_samples, wb_cs, ...
            sprintf('Coarse search: candidate %d/%d  ||e||=%.2e', ...
                    j, num_samples, min_residual));
    end
      
    % Step #3: Cleanup waitbar and finalize output
    close(wb_cs);
    fprintf('Best candidate ||e|| = %.6e\n', min_residual);
end
% -----------------------------------------------------------------------

function F = bvp_residual(guess, prm)
% AIM: Computes the 8-component shooting residual for the PMP boundary value problem (BVP).
% INPUTS:
%   guess : [8x1] Unknowns vector[lambda_r0(3); lambda_v0(3); lambda_m0; tf*].
%   prm   : [struct] Problem parameters (initial/target states, dynamics constants).
% OUTPUT:
%   F     : [8x1] Residual vector. Convergence to F=0 enforces terminal 
%           conditions (pos/vel) and transversality conditions (mass/Hamiltonian).

   % Step #1: Unpack decision variables from the guess vector
    p_r0 = guess(1:3);
    p_v0 = guess(4:6);
    p_m0 = guess(7);
    tf   = guess(8);

   % Step #2: Validate flight time and handle non-physical regions
    if tf <= 0
        F = 1e8 * ones(8,1);
        return;
    end
    % Step #3: Integrate augmented system (states and costates) forward to tf*
    Y0   = [prm.start_state; p_r0; p_v0; p_m0];
    opts = odeset('RelTol',1e-10, 'AbsTol',1e-13);
    flag = true;
    [~, Y] = ode45(@(t,YY) eom_augmented(t,YY,prm), [0 tf], Y0, opts);
    if any(~isfinite(Y(end,:))) || any(isnan(Y(end,:)))
        flag = false;
    end

    if ~flag
        F = 1e8 * ones(8,1);
        return;
    end
    % Step #4: Extract terminal state/costates and evaluate target conditions
    Yf   = Y(end,:).';
    rf   = Yf(1:3);  vf  = Yf(4:6);  mf  = Yf(7);
    p_mf = Yf(14);

    % Step #5: Compute terminal Hamiltonian (Transversality condition) and residuals
    final_ham = eval_H(Yf, prm);
    pos_resid=rf  - prm.end_state(1:3);
    vel_resid=vf  - prm.end_state(4:6);

    F = [pos_resid;
         vel_resid;
         p_mf;
         final_ham];
end
