%%  Spacecraft Guidance and Navigation 2025/2026
%  Assignment 1 – Exercise 1: Lambert Problem in the CR3BP
%  Author : Karim Ahmed Mohamed Doubie
%
%  Objective:
%   Point 1 – Solve the Lambert problem for Δt = 0.9 via STM-based
%            differential correction. Report the corrected initial velocity
%             and the terminal arrival velocity.
%            Plot the converged arc .
%   Point 2 – Apply numerical continuation from Δt = 0.9 to Δt = 0.3 to generate 
%            a one-parameter family of connecting arcs. Report the initial  and 
%            arrival velocities at Δt = 0.3. Plot the trajectory family .
%   Point 3 – Find Δt_min that minimises the injection manoeuvre Δv = ||v_target – vf(Δt)||.
%% STEP #0: INITIAL SETTING AND INITIAL DEFINITION
clc; clearvars; close all;

set(groot, 'DefaultAxesFontSize',                  16);
set(groot, 'DefaultAxesTitleFontSizeMultiplier',    1);
set(groot, 'DefaultTextFontSize',                  12);
set(0,     'DefaultTextInterpreter',         'latex');
set(0,     'DefaultLegendInterpreter',       'latex');
set(0,     'DefaultAxesTickLabelInterpreter','latex');
cspice_kclear();
cspice_furnsh('assignment01.tm');
format long g;

% Step #0.1: Dimensional conversion factors
DU          = 384400;                          % [km]    Distance Unit (Earth-Moon)
GM_Eearth   =  cspice_bodvrd('EARTH','GM',1);            
GM_Moon        = cspice_bodvrd('MOON' ,'GM',1);       
TU          = sqrt(DU^3 / (GM_Eearth + GM_Moon));      % [s]     Time Unit
TU_days     = TU / 86400;                      % [days]
VU          = DU / TU;                         % [km/s]  Velocity Unit

fprintf('Scale factors (DE432):  DU = %.0f km   TU = %.8f days   VU = %.8f km/s\n\n', ...
        DU, TU_days, VU);



%% STEP #1: LAMBERT PROBLEM SOLUTION (NOMINAL CASE Δt = 0.9)
fprintf('=============================================================\n');
fprintf('  PART 1  –  Lambert Problem   (Δt = 0.9)\n');
fprintf('=============================================================\n\n');

% Step #1.1: Define boundary conditions and initial guess
mu  = 0.012150;            % CR3BP mass ratio
r1  = [-0.3;  0;    0  ];  % Initial position  [adim]
r2  = [ 0.008;-0.01; 0  ]; % Target  position  [adim]
tof = 0.9;                 % Time of flight    [adim]
x0  = [r1; 0; 1; 0];       % Initial state guess [r1; v_guess]

fprintf('r1 = [%.4f  %.4f  %.4f]\n', r1);
fprintf('r2 = [%.4f  %.4f  %.4f]\n', r2);
fprintf('tof = %.4f\n\n', tof);


% Step #0.2: Differential corrector settings
Nmax    = 70;       % Maximum Newton-Raphson iterations
tolerance_position = 1e-11;    % Position convergence tolerance

% Step #1.2: Newton-Raphson differential corrector loop
i = 0;  err_r = Inf;
while i < Nmax
    i = i + 1;
    [xf, PHIf] = integrate_cr3bp_stm([r1; x0(4:6)], tof, mu);
    err_r = r2 - xf(1:3);
    if norm(err_r) <= tolerance_position, break; end
    % Core shooting update: dv0 = PHIf_rv \ err_r
    dv0       = PHIf(1:3,4:6) \ err_r;
    x0(4:6)   = x0(4:6) + dv0;
    fprintf('  Iter %2d   |F| = %.3e\n', i, norm(err_r));
end

% Step #1.3: Convergence check and result extraction
if norm(err_r) <= tolerance_position
    fprintf('Converged in %d iterations.\n', i);
else
    warning('Part 1: did not converge.');
end

% Store Part 1 solution
v_departure = x0(4:6);         % Converged initial velocity
v_arrival = xf(4:6);         % Arrival velocity at tf

% integrate_cr3bp_stm for plotting
[~, ~, xx] = integrate_cr3bp_stm([r1; v_departure], tof, mu);

% Result Printing – Part 1
fprintf('\n  Corrected initial velocity  v1  (Δt = 0.9)  [-]\n');
fprintf('  %-22s %-22s %-22s\n', 'vx [-]', 'vy [-]', 'vz [-]');
fprintf('  %+.10f      %+.10f      %+.10f\n\n', v_departure(1), v_departure(2), v_departure(3));

fprintf('Table 2 – Terminal velocity  vf(tf)  (Δt = 0.9)  [-]\n');
fprintf('  %-22s %-22s %-22s\n', 'vx(tf) [-]', 'vy(tf) [-]', 'vz(tf) [-]');
fprintf('  %+.10f      %+.10f      %+.10f\n\n', v_arrival(1), v_arrival(2), v_arrival(3));

% Visualization 
figure('Color','w', 'Name','Lambert solution');
plot(xx(:,1), xx(:,2), 'b-', 'LineWidth',1.5, 'DisplayName','Trajectory');
hold on; grid on; axis equal;
plot(r1(1), r1(2), 'ko', 'MarkerFaceColor','g', 'MarkerSize',8, 'DisplayName','$r_1$');
plot(r2(1), r2(2), 'ko', 'MarkerFaceColor','m', 'MarkerSize',8, 'DisplayName','$r_2$');
xlabel('x [-]'); ylabel('y [-]');
title(sprintf('Lambert arc in the CR3BP synodic frame, $\\Delta t = %.1f$ [-]', tof));
legend('Location','best');

%% STEP #2: NUMERICAL CONTINUATION (Δt : 0.9 to 0.3)
fprintf('=============================================================\n');
fprintf('  PART 2  –  Numerical Continuation   (0.9 to 0.3)\n');
fprintf('=============================================================\n\n');

% Step #2.1: Initialise tof grid and storage matrices
tof_grid1      = linspace(0.9, 0.3, 60);   % tof range [adim]
n_tof  = length(tof_grid1);

V0     = zeros(3, n_tof);           % [3×n_tof] initial velocity at each tof
VF     = zeros(3, n_tof);           % [3×n_tof] arrival velocity at each tof
traj_p2 = cell(n_tof, 1);

V0(:,1) = v_departure;                    % seed from Part 1
colors  = parula(n_tof);            

% Step #2.2: Continuation loop (warm-started Newton-Raphson)
for s = 1:n_tof
    tof_f = tof_grid1(s);
    v1c   = V0(:, max(s-1,1));
    res = Inf;

    for i = 1:Nmax
        [xf, PHIf] = integrate_cr3bp_stm([r1; v1c], tof_f, mu);
        res = r2 - xf(1:3);
        if norm(res) <= tolerance_position, break; end
        v1c = v1c + PHIf(1:3,4:6) \ res;
    end

    if norm(res) <= tolerance_position
        V0(:,s) = v1c;
        VF(:,s) = xf(4:6);
        [~,~,xx]     = integrate_cr3bp_stm([r1; v1c], tof_f, mu);
        traj_p2{s} = xx;
    else
        fprintf('Part 2: corrector failed at tof = %.4f (step %d)\n', tof_f, s);
        break;
    end
end

% Step #2.3: Report terminal solution at Δt = 0.3
if all(V0(:,end) == 0), error('Part 2: last step did not converge.'); end

fprintf(' Corrected initial velocity  v1  (Δt = 0.3)  [-]\n');
fprintf('  %-22s %-22s %-22s\n', 'vx [-]', 'vy [-]', 'vz [-]');
fprintf('  %+.10f      %+.10f      %+.10f\n\n', V0(1,end), V0(2,end), V0(3,end));

fprintf(' Terminal velocity  vf(tf)  (Δt = 0.3)  [-]\n');
fprintf('  %-22s %-22s %-22s\n', 'vx(tf) [-]', 'vy(tf) [-]', 'vz(tf) [-]');
fprintf('  %+.10f      %+.10f      %+.10f\n\n', VF(1,end), VF(2,end), VF(3,end));
% Terminal position at Δt = 0.3 (should match r2 to tolerance)
[xf_end, ~] = integrate_cr3bp_stm([r1; V0(:,end)], tof_grid1(end), mu);
rf_end = xf_end(1:3);
fprintf('Terminal position  rf(tf)  (Δt = 0.3)  [-]\n');
fprintf('  %-22s %-22s %-22s\n', 'rx [-]', 'ry [-]', 'rz [-]');
fprintf('  %+.10f      %+.10f      %+.10f\n\n', rf_end(1), rf_end(2), rf_end(3));
fprintf('  ||rf - r2|| = %.3e  [-]\n\n', norm(rf_end - r2));

% Visualization 
figure('Color','w','Name','Trajectory family');
hold on; grid on; axis equal;
for s = 1:n_tof
    if ~isempty(traj_p2{s})
        plot(traj_p2{s}(:,1), traj_p2{s}(:,2), '-', ...
             'Color',colors(s,:), 'LineWidth',1.2, 'HandleVisibility','off');
    end
end
plot(NaN,NaN,'-','Color',[0.45 0.45 0.45],'LineWidth',1.5, ...
     'DisplayName','Trajectory family ($\Delta t \downarrow$)');
plot(r1(1),r1(2),'ko','MarkerFaceColor','g','MarkerSize',8,'DisplayName','$r_1$');
plot(r2(1),r2(2),'ko','MarkerFaceColor','m','MarkerSize',8,'DisplayName','$r_2$');
xlabel('x [-]'); ylabel('y [-]');
title('CR3BP Lambert family: continuation in time of flight');
legend('Location','best');

%% STEP #3: MINIMUM INJECTION MANOEUVRE
fprintf('=== PART 3 – Minimum Injection Manoeuvre ===\n\n');

v_LEO = [-2.56; -5.16; 3.31];

% Step #3.1: Target LEO velocity and refinement configuration
% --- Search configuration ---
n_coarse       = 750;        % coarse sweep resolution
n_refine       = 250;        % refinement resolution
shrink_factor  = 0.20;       % bracket = 20 % of current interval around minimum
max_zoom       = 5;          % maximum refinement passes
tol_bracket    = 1e-10;      % stop when interval width < this

% --- Storage ---
rec_tof = [];
rec_dv  = [];

% --- Initialise search bounds and warm seed ---
alpha_lo = 0.001;
alpha_hi = 0.300;
v1_warm  = V0(:, end);       % seed from Part 2 terminal solution

zoom_pass = 0;
proceed   = true;

while proceed
    % Step #3.2: Build search grid for this level
    if zoom_pass == 0
        %Step#3.2.1 Coarse sweep: descending tof (better continuation direction)
        tof_grid = linspace(alpha_hi, alpha_lo, n_coarse);
        fprintf('Coarse sweep  –  %d pts  [%.4f → %.4f]\n', n_coarse, alpha_hi, alpha_lo);
        n_pts = n_coarse;
    else
        tof_grid = linspace(alpha_lo, alpha_hi, n_refine);
        fprintf('Zoom pass %d   –  %d pts  [%.10f → %.10f]\n', zoom_pass, n_refine, alpha_lo, alpha_hi);
        n_pts = n_refine;
    end
    success_grid   = false(1, n_pts); v1_cursor = v1_warm;
    dv_grid   = NaN(1, n_pts); v1_grid   = NaN(3, n_pts);
   
    % Step#3.3: Sweep with warm-started corrector 
    for s = 1:n_pts
        tof_s = tof_grid(s);
        v1c   = v1_cursor;
        res   = Inf;

        for i = 1:Nmax
            [xf, PHIf] = integrate_cr3bp_stm([r1; v1c], tof_s, mu);
            res = r2 - xf(1:3);
            if norm(res) <= tolerance_position, break; end
            v1c = v1c + PHIf(1:3,4:6) \ res;
        end

        if norm(res) <= tolerance_position
            dv_grid(s)   = norm(v_LEO - xf(4:6));
            v1_grid(:,s) = v1c;
           success_grid(s)   = true;
            v1_cursor    = v1c;
        end
    end
    rec_tof = [rec_tof, tof_grid]; rec_dv  = [rec_dv,  dv_grid];

    % Step#3.4:Find minimum among converged points
    conv_idx =success_grid & ~isnan(dv_grid);
    if ~any(conv_idx)
        fprintf('  No converged points – terminating.\n');
        break;
    end

    dv_conv_idx_local = dv_grid;
    dv_conv_idx_local(~conv_idx) = Inf;
    [~, idx_best_local] = min(dv_conv_idx_local);

    tof_star = tof_grid(idx_best_local);

    % Step#3.5: Geometric bracket: shrink_factor * interval width around minimum 
    half_w   = shrink_factor * (tof_grid(end) - tof_grid(1)) / 2;
    % handle descending grid
    if tof_grid(1) > tof_grid(end)
        half_w = shrink_factor * (tof_grid(1) - tof_grid(end)) / 2;
    end
    alpha_lo = max(0.001, tof_star - abs(half_w));
    alpha_hi = min(0.300, tof_star + abs(half_w));

    %Step#3.6: Warm-start next pass from LEFT boundary solution 
    idx_left = idx_best_local;   % find nearest converged point to alpha_lo
    [~, idx_left] = min(abs(tof_grid(conv_idx) - alpha_lo));
    idx_left_global = find(conv_idx);
    v1_warm = v1_grid(:, idx_left_global(idx_left));

    fprintf('  Best tof = %.10f   dv = %.10f\n', tof_star, dv_grid(idx_best_local));

    zoom_pass = zoom_pass + 1;
    if zoom_pass > max_zoom || abs(alpha_hi - alpha_lo) < tol_bracket
        proceed = false;
    end
end

% --- Extract global minimum ---
[tof_all, sort_idx] = sort(rec_tof);
dv_all   = rec_dv(sort_idx);
clean_mask     = ~isnan(dv_all);
tof_all  = tof_all(clean_mask);
dv_all   = dv_all(clean_mask);

[dv_min, idx_global_min]  = min(dv_all);
tof_min          = tof_all(idx_global_min);
tof_min_dimensional     = tof_min * TU_days;
dv_min_dimensional       = dv_min  * VU;

fprintf('\n--- Results ---\n');
fprintf('  %-22s %-22s %-22s %-22s\n','dt_min [-]','dv_min [-]','dt_min [days]','dv_min [km/s]');
fprintf('  %-22.10f %-22.10f %-22.10f %-22.10f\n\n', tof_min, dv_min, tof_min_dimensional, dv_min_dimensional);

% --- Plot (full scale) ---
tof_dimensional_all = tof_all * TU_days;
dv_dimensional_all   = dv_all  * VU;

figure('Color','w');
hold on; grid on;
plot(tof_dimensional_all, dv_dimensional_all, 'Color',[1 0.5 0], 'LineWidth', 1.4, 'DisplayName','$\Delta v(\Delta t)$');
plot(tof_min_dimensional, dv_min_dimensional, 'bx', 'MarkerSize', 9, 'LineWidth', 2, 'DisplayName','$\Delta t_{\min}$');
xlabel('$\Delta t$ [days]'); ylabel('$\Delta v$ [km/s]');
title('Injection manoeuvre vs.\ Time of Flight');
legend('Location','best');

% --- Plot (zoomed) ---
figure('Color','w');
hold on; grid on;
plot(tof_dimensional_all, dv_dimensional_all, 'Color',[1 0.5 0], 'LineWidth', 1.4, 'DisplayName','$\Delta v(\Delta t)$');
plot(tof_min_dimensional, dv_min_dimensional, 'bx', 'MarkerSize', 9, 'LineWidth', 2, 'DisplayName','$\Delta t_{\min}$');
xlabel('$\Delta t$ [days]'); ylabel('$\Delta v$ [km/s]');
title('Injection manoeuvre vs.\ Time of Flight (zoom)');
ylim([dv_min_dimensional*0.90, min(prctile(dv_dimensional_all,95), 2.5*dv_min_dimensional)]);
xlim([min(tof_dimensional_all), prctile(tof_dimensional_all,85)]);
legend('Location','best');
%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function [xf, PHIf, xx] = integrate_cr3bp_stm(x0, tof, mu)
%integrate_cr3bp_stm  Integrates the 42-DOF augmented state (state + STM) in the CR3BP.
%   Uses ode113 (high-order Runge-Kutta) for efficient, accurate propagation.
%
%   Inputs:
%       x0  : [6x1]  Initial state [r; v]
%       tof : [1x1]  Time of flight [adim]
%       mu  : [1x1]  CR3BP mass ratio
%   Outputs:
%       xf   : [6x1]   Final state
%       PHIf : [6x6]   Final State Transition Matrix
%       xx   : [Nx6]   (Optional) full state history

    Phi0   = eye(6);
    x0_aug = [x0; Phi0(:)];

    options = odeset('RelTol',1e-13,'AbsTol',1e-14);
    [~, xx_aug] = ode113(@(t,x) xyzCR3BP_STM(t, x, mu), [0 tof], x0_aug, options);

    xf_aug = xx_aug(end,:).';
    xf     = xf_aug(1:6);
    PHIf   = reshape(xf_aug(7:end), 6, 6);

    if nargout > 2
        xx = xx_aug(:,1:6);
    end
end

function dxdt = xyzCR3BP_STM(~, x_aug, mu)
%XYZCR3BP_STM  RHS of the CR3BP + variational equations for the STM.
%   Computes the 6-DOF equations of motion and A(t)·Φ analytically.
%
%   Inputs:
%       x_aug : [42x1]  Augmented state [r; v; Phi(:)]
%       mu    : [1x1]   CR3BP mass ratio
%   Output:
%       dxdt  : [42x1]  Augmented state derivative

    % --- Unpack state and STM ---
    r   = x_aug(1:3);
    v   = x_aug(4:6);
    Phi = reshape(x_aug(7:end), 6, 6);

    x = r(1);  y = r(2);  z = r(3);

    % --- Distances to primaries ---
    d1 = sqrt((x+mu)^2    + y^2 + z^2);
    d2 = sqrt((x-(1-mu))^2+ y^2 + z^2);

    % --- Equations of motion (state derivatives) ---
    ax = 2*v(2) + x - (1-mu)*(x+mu)/d1^3 - mu*(x-(1-mu))/d2^3;
    ay = -2*v(1)+ y - (1-mu)*y/d1^3      - mu*y/d2^3;
    az =        - (1-mu)*z/d1^3           - mu*z/d2^3;

    dxdt_state = [v; ax; ay; az];

    % --- Effective potential second derivatives ---
    dUdx = 1 - (1-mu)/d1^3 + 3*(1-mu)*(x+mu)^2/d1^5 ...
             - mu/d2^3      + 3*mu*(x-(1-mu))^2/d2^5;
    dUdy = 1 - (1-mu)/d1^3 + 3*(1-mu)*y^2/d1^5 ...
             - mu/d2^3      + 3*mu*y^2/d2^5;
    dUdz =   - (1-mu)/d1^3 + 3*(1-mu)*z^2/d1^5 ...
             - mu/d2^3      + 3*mu*z^2/d2^5;

    Uxy = 3*(1-mu)*(x+mu)*y/d1^5 + 3*mu*(x-(1-mu))*y/d2^5;
    Uxz = 3*(1-mu)*(x+mu)*z/d1^5 + 3*mu*(x-(1-mu))*z/d2^5;
    Uyz = 3*(1-mu)*y*z/d1^5       + 3*mu*y*z/d2^5;

    % --- Jacobian A(t) = df/dx  [6×6] ---
    dfdx = zeros(6,6);
    dfdx(1:3, 4:6) = eye(3);
    dfdx(4:6, 1:3) = [dUdx, Uxy,  Uxz;
                       Uxy,  dUdy, Uyz;
                       Uxz,  Uyz,  dUdz];
    dfdx(4,5) =  2;   % Coriolis da_x/dv_y
    dfdx(5,4) = -2;   % Coriolis da_y/dv_x

    % --- STM derivative: Φ̇ = A(t)·Φ ---
    Phidot = dfdx * Phi;

    dxdt = [dxdt_state; Phidot(:)];
end
