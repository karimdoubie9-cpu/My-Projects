%% Spacecraft Guidance and Navigation (2025/2026)
% Assignment 2 – Exercise 3: Sequential Filtering for FUTURE CubeSat
% Author: Karim Ahmed Mohamed Doubie
%
% Objective:
%   Point 1 – Propagate J2-perturbed trajectory and identify AI-IP acquisition
%              windows from the provided CSV measurement epochs.
%   Point 2 – Simulate GNSS measurements by injecting 20 m (1-sigma) Gaussian
%              noise on the SPK truth and verify 3-sigma consistency.
%   Point 3 – UKF with GNSS: process continuous Cartesian position updates
%              via Unscented Transform under J2-perturbed dynamics.
%   Point 4 – UKF with AI-IP: process sparse geodetic measurements (Range,
%              Lon, Lat) with updates restricted to acquisition windows and
%              circular statistics for longitude periodicity.
%   Point 5 – Compare both filters against SPK truth in terms of RMS/final
%              errors and 3-sigma covariance consistency.
%% STEP #0: INITIAL SETTING AND INITIAL DEFINITION
clc; clearvars; close all;

% Step #0.1: Set up paths for external functions and SPICE kernels
addpath('kernels/');
addpath('sgp4/');

% Step #0.2: Configure global LaTeX interpreter for professional plotting
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');
% Set the default font size for all axes (includes tick labels and titles)
set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1);
% Optional: Set the font size for all other text objects (legends, text annotations)
set(groot, 'DefaultTextFontSize', 12);
% Step #0.3: Load the mission meta-kernel
cspice_furnsh('assignment02.tm');
rng(40); % Seed for reproducibility
%Given data:
    %noise for GNSS sensor
    std_pos_m  = 20;               
    std_pos_km = std_pos_m * 1e-3; 
    COV_GNSS = (std_pos_km^2) * eye(3);
    %noise for AI-IP geodetic sensor
    sigma_rng = 10.12;              
    sigma_ang = deg2rad(0.07);      
    COV_AIP = diag([sigma_rng^2, sigma_ang^2, sigma_ang^2]);

%% STEP #1: NOMINAL TRAJECTORY & ACQUISITION WINDOW ANALYSIS (POINT 1)
% Step #1.1: Define simulation timeframe and initial state from Table 6
t0_utc = '2027-05-02T01:00:00';t0 = cspice_str2et(t0_utc);
tf_utc = '2027-05-02T12:00:00';tf = cspice_str2et(tf_utc);
% Initial state in ECI J2000 [km; km/s]
x0 = [-436.489; 872.269; -6904.868; 7.126589; 2.598426; -0.131473];

% Step #1.2: Initialize time grid and gravitational constants
time_step = 5.0;
tspan = (t0 : time_step : tf)';   
N = length(tspan);

lbl_i = 'Propagation Summary';
met_i = ["Duration [h]"; "N"; "t0 [UTC]"; "tf [UTC]"];
val_i = [string((tf-t0)/3600); string(N); string(t0_utc); string(tf_utc)];

fprintf('\n=== %s ===\n', lbl_i);
disp(table(met_i, val_i, 'VariableNames', {'Metric', 'Value'}));
% Step #1.3: Propagate nominal J2 trajectory
fprintf('Propagating nominal J2 trajectory...\n');
[~, ~, trajectory_eci] = keplerian_propagator(x0, tspan, 'EARTH', 'J2');

% Step #1.4: Compute geodetic ground track (Lon/Lat)
[lon_deg, lat_deg] = compute_ground_track(tspan, trajectory_eci);

% Step #1.5: Load AI-IP measurement data and map to simulation time grid
fprintf('Loading ai_ip_measurements.csv...\n');
ai_data       = readtable('ai_ip_measurements.csv');
ai_epochs_et  = cspice_str2et(cellstr(string(ai_data.timetag)));



lon_obs = wrapToPi(deg2rad(ai_data.lon(:)).');%rad
lat_obs = deg2rad(ai_data.lat(:)).';%rad
range_obs  = ai_data.range(:).';%km

y_aip = nan(3, N);   % [Range; Lon; Lat]
ai_visible = false(1, N);
% time grid mapping 
idx_obs = dsearchn(tspan(:), ai_epochs_et(:));
for j = 1:numel(idx_obs)
    k = idx_obs(j);
    y_aip(:,k) = [range_obs(j); lon_obs(j); lat_obs(j)];
    ai_visible(k) = true;
end
% Step #1.6: Identify discrete acquisition window intervals
t_start = []; t_end = [];
if any(ai_visible)
    d       = diff([0 ai_visible 0]);
    t_start = find(d ==  1);
    t_end   = find(d == -1) - 1;
end

% --- Result Printing: Point 1 ---
fprintf('\n--- Point 1: AI-IP Acquisition Data ---\n');
T_summary = table([numel(ai_epochs_et); sum(ai_visible); numel(t_start)], ...
    'VariableNames', {'Value'}, ...
    'RowNames', {'Rows_in_CSV', 'Mapped_Updates_on_Grid', 'Acquisition_Windows'});
disp(T_summary);

num_windows = numel(t_start);
Window       = (1:num_windows)';
Start_UTC    = cell(num_windows, 1);
End_UTC      = cell(num_windows, 1);
Duration_min = zeros(num_windows, 1);

for idx = 1:num_windows
    t_s = tspan(t_start(idx));
    t_e = tspan(t_end(idx));
    Start_UTC{idx}    = cspice_et2utc(t_s, 'ISOC', 0);
    End_UTC{idx}      = cspice_et2utc(t_e, 'ISOC', 0);
    Duration_min(idx) = (t_e - t_s) / 60;
end

T_intervals = table(Window, Start_UTC, End_UTC, Duration_min);
fprintf('Detailed Acquisition Intervals:\n');
disp(T_intervals);

% --- Plotting: Point 1 ---
lon_plot = lon_deg;  lat_plot = lat_deg;
jumps = find(abs(diff(lon_deg)) > 90);
lon_plot(jumps+1) = NaN;
lat_plot(jumps+1) = NaN;

figure('Name',' Ground track');
hold on; grid on;
if isfile('Earth_IMAGE.jpg')
    image([-180 180],[90 -90],imread('Earth_IMAGE.jpg'));
    set(gca,'YDir','normal');
end
h_tr  = plot(lon_plot(~ai_visible), lat_plot(~ai_visible), '.', ...
             'Color',[0.6 0.6 0.6],'MarkerSize',4,'DisplayName','Trajectory');
h_aq  = plot(lon_plot(ai_visible),  lat_plot(ai_visible),  '.', ...
             'Color',[0.85 0.33 0.10],'MarkerSize',8,'DisplayName','Acquisition windows');
h_s   = plot(lon_plot(1),   lat_plot(1),   '*','Color','g','MarkerSize',14,'DisplayName','Start');
h_e   = plot(lon_plot(end), lat_plot(end), '*','Color','r','MarkerSize',14,'DisplayName','End');
xlim([-180 180]); ylim([-90 90]);
xlabel('Longitude $[^\circ]$'); ylabel('Latitude $[^\circ]$');
legend([h_tr,h_aq,h_s,h_e],'Location','southout','Orientation','horizontal');

fprintf('--- Point 1 Complete ---\n');

%% STEP #2: GNSS MEASUREMENT SIMULATION & STOCHASTIC MODELING (POINT 2)
fprintf('\n--- POINT 2: GNSS ---\n');

% Step #2.1: Retrieve SPK reference "Truth" trajectory for error analysis
target_id = '-10003001';
state_reference = cspice_spkezr(target_id, tspan', 'J2000','NONE','EARTH');  
r_reference = state_reference(1:3,:);   
v_reference = state_reference(4:6,:);   

% Step #2.2: Define measurement noise model and generate noisy observations
r_GNSS = r_reference + std_pos_km * randn(3, N);   
% Step #2.3: Validate noise properties and compute error statistics
error_GNSS = (r_GNSS - r_reference)* 1e3;mean_error    = mean(error_GNSS, 2);
rms_error    = std(error_GNSS, 0, 2);valid_ratio = mean(abs(error_GNSS) <= 3*std_pos_m, 2);

% --- Result Printing: Point 2 ---
fprintf('\n--- Point 2: GNSS Noise Statistics Summary ---\n');
fprintf('Target sigma_p = %.1f m\n\n', std_pos_m);
T_gnss_stats = table(mean_error(:), rms_error(:), valid_ratio(:), ...
    'VariableNames', {'Mean_m', 'Std_m', 'Frac_within_3sigma'}, ...
    'RowNames', {'X_axis', 'Y_axis', 'Z_axis'});
disp(T_gnss_stats);

% --- Plotting: Point 2 ---
figure('Name',' GNSS noise check');
labs = {'$|\mathbf{r}_{\mathrm{true},x} - \mathbf{r}_{\mathrm{GNSS},x}|$', ...
        '$|\mathbf{r}_{\mathrm{true},y} - \mathbf{r}_{\mathrm{GNSS},y}|$', ...
        '$|\mathbf{r}_{\mathrm{true},z} - \mathbf{r}_{\mathrm{GNSS},z}|$'};
c_noise = [0.85, 0.33, 0.10];
c_sigma = [0, 0.45, 0.74];
fs = 13;
t_delta = tspan - t0;
for i = 1:3
    subplot(3,1,i); hold on; grid on;
    plot(t_delta, abs(error_GNSS(i,:)), '.', 'Color', c_noise, 'DisplayName', labs{i});
    yline(3*std_pos_m, 'Color', c_sigma, 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', '$+3\sigma$');
    ylabel([labs{i},' [m]'], 'FontSize', fs);
    set(gca, 'FontSize', fs);
    if i==3, xlabel('$t-t_0$ [s]', 'FontSize', fs); end
    legend('Location','best', 'FontSize', fs);
end
sgtitle('Component-wise GNSS position error and $3\sigma$ noise threshold', 'FontSize', fs);
fprintf('--- Point 2 Complete ---\n');
%% STEP #3: SEQUENTIAL FILTERING - GNSS UKF (POINT 3)
fprintf('\n--- POINT 3: UKF GNSS ---\n');

% Step #3.1: Initialize filter state, covariance, and UT tuning parameters
x0_filter = x0;
P0        = diag([25,25,25,0.005,0.005,0.005]);
alpha=0.01; beta=2; kappa=0;
attractor='EARTH'; perturbations='J2';

% Step #3.2: Run Unscented Kalman Filter with Cartesian GNSS measurements
[x_est_gnss, P_est_gnss] = UKF_GNSS(x0_filter, P0, tspan, r_GNSS', ...
                                      alpha, beta, kappa, COV_GNSS, ...
                                      attractor, perturbations);
fprintf('--- Point 3 Complete ---\n');

%% STEP #4: SEQUENTIAL FILTERING - AI-IP UKF (POINT 4)
fprintf('\n--- POINT 4: UKF AI-IP ---\n');

% Step #4.1: Run UKF with conditional updates (active only in acquisition windows)
[x_est_ai, P_est_ai] = UKF_AIIP(x0_filter, P0, tspan, y_aip, ai_visible, ...
                                  alpha, beta, kappa, COV_AIP, attractor, perturbations);
fprintf('--- Point 4 Complete ---\n');

%% STEP #5: PERFORMANCE EVALUATION & COMPARATIVE ANALYSIS (POINT 5)
fprintf('\n--- POINT 5: Compare ---\n');

% Step #5.1: Calculate component-wise residuals and 3-sigma diagonal bounds
e_position_AIP  = x_est_ai(1:3,:)   - r_reference;
e_velocity_AIP  = x_est_ai(4:6,:)   - v_reference;
e_Position_GNSS = x_est_gnss(1:3,:) - r_reference;
e_Velocity_GNSS = x_est_gnss(4:6,:) - v_reference;

ThreeSigma_position_AIP  = nan(3,N); 
ThreeSigma_Velocity_AIP  = nan(3,N);
ThreeSigma_Position_GNSS = nan(3,N); 
ThreeSigma_Velocity_GNSS = nan(3,N);
for j = 1:N
    ThreeSigma_position_AIP(:,j)  = 3*sqrt(diag(P_est_ai(1:3,1:3,j)));
    ThreeSigma_Velocity_AIP(:,j)  = 3*sqrt(diag(P_est_ai(4:6,4:6,j)));
    ThreeSigma_Position_GNSS(:,j) = 3*sqrt(diag(P_est_gnss(1:3,1:3,j)));
    ThreeSigma_Velocity_GNSS(:,j) = 3*sqrt(diag(P_est_gnss(4:6,4:6,j)));
end

% Step #5.2: Compute aggregate norm errors for filter benchmarking
error_position_GNSS = vecnorm(e_Position_GNSS,2,1);
error_velocity_GNSS = vecnorm(e_Velocity_GNSS,2,1);
error_position_AIP = vecnorm(e_position_AIP, 2,1);
error_velocity_AIP = vecnorm(e_velocity_AIP, 2,1);

% --- Result Printing: Points 3, 4, & 5 ---
fprintf('Filter Performance Comparison \n ');
rms_r_g = sqrt(mean(error_position_GNSS.^2));
rms_v_g = sqrt(mean(error_velocity_GNSS.^2));
rms_r_a = sqrt(mean(error_position_AIP.^2));
rms_v_a = sqrt(mean(error_velocity_AIP.^2));

T_perf = table([error_position_GNSS(end); error_position_AIP(end)], [rms_r_g; rms_r_a], ...
               [error_velocity_GNSS(end); error_velocity_AIP(end)], [rms_v_g; rms_v_a], ...
    'VariableNames', {'Final_Pos_Err_km', 'RMS_Pos_Err_km', ...
                      'Final_Vel_Err_kms', 'RMS_Vel_Err_kms'}, ...
    'RowNames', {'GNSS_Point_3', 'AIIP_Point_4'});
disp(T_perf);

fprintf('AI-IP Acquisition Window Analysis \n');
in  = ai_visible;  out = ~ai_visible;
rms_r_in  = sqrt(mean(error_position_AIP(in).^2));
rms_r_out = sqrt(mean(error_position_AIP(out).^2));
rms_v_in  = sqrt(mean(error_velocity_AIP(in).^2));
rms_v_out = sqrt(mean(error_velocity_AIP(out).^2));

T_window_analysis = table([rms_r_in; rms_r_out], [rms_v_in; rms_v_out], ...
    'VariableNames', {'Pos_RMS_km', 'Vel_RMS_kms'}, ...
    'RowNames', {'in_Windows', 'out_Windows'});
disp(T_window_analysis);


fprintf('GNSS Acquisition Window Analysis \n ');
rms_r_in_g  = sqrt(mean(error_position_GNSS(in).^2));
rms_r_out_g = sqrt(mean(error_position_GNSS(out).^2));
rms_v_in_g  = sqrt(mean(error_velocity_GNSS(in).^2));
rms_v_out_g = sqrt(mean(error_velocity_GNSS(out).^2));

T_window_gnss = table([rms_r_in_g; rms_r_out_g], [rms_v_in_g; rms_v_out_g], ...
    'VariableNames', {'Pos_RMS_km', 'Vel_RMS_kms'}, ...
    'RowNames', {'in_Windows', 'out_Windows'});
disp(T_window_gnss);


% Step #5.3: Visual Comparison - Log-scale Norms and Component Residuals
et_ai_srt = sort(ai_epochs_et(:));
d_ai      = diff(et_ai_srt);
brk_ai    = find(d_ai > time_step + 1e-6);
w_bounds  = [0; brk_ai(:); numel(et_ai_srt)];
t_mjd     = tspan / cspice_spd();

% Formatting
fl_r = 1e-7;  fl_v = 1e-9; 
c_gnss_err = [0.8 0 0];    c_gnss_sig = [0 0.5 0];    
c_ai_err   = [0.5 0 0.5];  c_ai_sig   = [0 0.4 0.6];  
c_mkr      = [0 0.9 0.9];  win_col    = [0.9 0.9 0.9]; 
lw_e = 1.2; lw_s = 1.5;
% 

ThreeSigma_Position_GNSS_mean = max(vecnorm(ThreeSigma_Position_GNSS,2,1), fl_r);
ThreeSigma_Velocity_GNSS_mean = max(vecnorm(ThreeSigma_Velocity_GNSS,2,1), fl_v);
ThreeSigma_position_AIP_mean  = max(vecnorm(ThreeSigma_position_AIP, 2,1), fl_r);
ThreeSigma_Velocity_AIP_mean  = max(vecnorm(ThreeSigma_Velocity_AIP, 2,1), fl_v);
% ---  Norm errors (Log Scale) ---
figure('Name',' Norm errors (log)');

t_hrs = (tspan - tspan(1)) / 3600;
win_hrs = (et_ai_srt - tspan(1)) / 3600;

win_col = [0.95, 0.85, 0.45]; alpha_val = 0.5;

for i = 1:2
    subplot(2,1,i);
    hold on; grid on; grid minor; set(gca,'YScale','log');
    
    if i == 1
        cur_err_g = error_position_GNSS; cur_err_a = error_position_AIP; 
        cur_sig_g = ThreeSigma_Position_GNSS_mean; cur_sig_a = ThreeSigma_position_AIP_mean; 
        cur_fl = fl_r; y_lims = [1e-6 1e2]; 
        lbl = '$\|e_r\|$ [km]'; % Updated Label
    else
        cur_err_g = error_velocity_GNSS; cur_err_a = error_velocity_AIP; 
        cur_sig_g = ThreeSigma_Velocity_GNSS_mean; cur_sig_a = ThreeSigma_Velocity_AIP_mean; 
        cur_fl = fl_v; y_lims = [1e-8 1e-1]; 
        lbl = '$\|e_v\|$ [km/s]'; % Updated Label
    end

    % Plot Acquisition Windows Background
    for ii = 2:numel(w_bounds)
        tw1 = win_hrs(w_bounds(ii-1)+1);
        tw2 = win_hrs(w_bounds(ii));
        fill([tw1 tw2 tw2 tw1], [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], ...
             win_col, 'EdgeColor','none', 'FaceAlpha',alpha_val, 'HandleVisibility','off');
    end

    uistack(findobj(gca,'Type','patch'),'bottom');
    
    % Plot Errors and 3-sigma Bounds
    plot(t_hrs, max(cur_err_g,cur_fl), 'Color',c_gnss_err, 'LineWidth',lw_e, 'DisplayName','GNSS');
    plot(t_hrs, max(cur_err_a,cur_fl), 'Color',c_ai_err,   'LineWidth',lw_e, 'DisplayName','AI IP');
    plot(t_hrs, cur_sig_g, '--', 'Color',c_gnss_sig, 'LineWidth',lw_s, 'DisplayName','$3\sigma$ GNSS');
    plot(t_hrs, cur_sig_a, '--', 'Color',c_ai_sig,   'LineWidth',lw_s, 'DisplayName','$3\sigma$ AI IP');

    ylim(y_lims); xlim([t_hrs(1) t_hrs(end)]); ylabel(lbl, 'Interpreter', 'latex'); 
    legend('Location','best', 'Interpreter', 'latex');
end
xlabel('$t - t_0$ [h]', 'Interpreter', 'latex');    
    fprintf('--- Point 5 Complete ---\n');
   %cspice_kclear();



%% =====================================================================
%  HELPER FUNCTIONS
%% =====================================================================
function [xf, tt, xx] = keplerian_propagator(x0, tspan, attractor, effect)
% PURPOSE: 
%   High-level interface for spacecraft state propagation. Selects between 
%   pure Keplerian (Two-Body) or J2-perturbed dynamics based on user input.
%
% INPUTS:
%   x0        : [6x1] Initial state [km; km/s]
%   tspan     : [Nx1] Ephemeris time grid [s]
%   attractor : [string] Central body name for SPICE (e.g., 'EARTH')
%   effect    : [string] Dynamics model ('J2' or 'Keplerian')
%
% OUTPUTS:
%   xf        : [6x1] Final state at tspan(end)
%   tt        : [Nx1] Output time vector
%   xx        : [Nx6] State history [km, km/s]

    % Step #1.1: Retrieve gravitational and geometric constants from SPICE
    mu_prop = cspice_bodvrd(attractor,'GM',1);  mu_prop = mu_prop(1);
    Re_prop = cspice_bodvrd(attractor,'RADII',3); Re_prop = Re_prop(1);

    % Step #1.2: Configure solver tolerances and integrate selected EOMs
    opts = odeset('RelTol',1e-12,'AbsTol',1e-13);
    if strcmp(effect,'J2')
        [tt, xx] = ode113(@(t,x) keplerianJ2_RHS(t,x,mu_prop,Re_prop), ...
                          tspan, x0(:), opts);
    else
        [tt, xx] = ode113(@(t,x) keplerian_RHS(t,x,mu_prop), ...
                          tspan, x0(:), opts);
    end
    xf = xx(end,:)';
end

% -----------------------------------------------------------------------
function dxdt = keplerian_RHS(~, x, mu)
% PURPOSE: 
%   Defines the Right-Hand Side (RHS) for standard Two-Body dynamics.
%
% Step #2.1: Extract position and compute Keplerian acceleration
    r = x(1:3);  rn = norm(r);
    dxdt = [x(4:6); -mu*r/rn^3];
end

% -----------------------------------------------------------------------
function dxdt = keplerianJ2_RHS(t, x, mu, Re)
% PURPOSE: 
%   Defines the RHS for J2-perturbed dynamics. Oblateness is calculated 
%   in the Earth-Fixed frame and rotated back to Inertial (J2000).
%
% Step #3.1: Calculate Central Body acceleration
    j2 = 0.0010826269;
    r = x(1:3);  rn = norm(r);
    dxdt = zeros(6,1);
    dxdt(1:3) = x(4:6);
    dxdt(4:6) = -mu/rn^3 * r;

% Step #3.2: Transform to ECEF and compute J2 perturbation
    R_eci2ecef = cspice_pxform('J2000','IAU_EARTH',t);
    r_ECEF = R_eci2ecef * r;
    a_j2_ecef = 1.5*mu*j2/rn^3*(Re/rn)^2 * ...
                r_ECEF .* (5*(r_ECEF(3)/rn)^2 - [1;1;3]);

% Step #3.3: Map J2 acceleration back to ECI and sum
    dxdt(4:6) = dxdt(4:6) + R_eci2ecef' * a_j2_ecef;
end

% -----------------------------------------------------------------------
function [lon_deg, lat_deg] = compute_ground_track(tspan, trajectory_eci)
% PURPOSE: 
%   Converts an Inertial ECI trajectory into geographic Lon/Lat coordinates.
%
% INPUTS:
%   tspan          : [Nx1] Time vector [s]
%   trajectory_eci : [Nx6] ECI state history [km, km/s]
%
% OUTPUTS:
%   lon_deg, lat_deg : Geographic coordinates [deg]
%
    % Step #5.1: Retrieve Earth shape parameters (Ellipsoidal f)
    radii = cspice_bodvrd('EARTH','RADII',3);
    Re = radii(1);  f = (radii(1)-radii(3))/radii(1);
    n = length(tspan);
    lon_deg = zeros(n,1);  lat_deg = zeros(n,1);

    % Step #5.2: Convert ECI to ECEF and extract geodetic coordinates per epoch
    for i = 1:n
        R = cspice_pxform('J2000','ITRF93', tspan(i));
        p_ecef = R * trajectory_eci(i,1:3)';
        [lon_r, lat_r, ~] = cspice_recgeo(p_ecef, Re, f);
        lon_deg(i) = rad2deg(wrapToPi(lon_r));
        lat_deg(i) = rad2deg(lat_r);
    end
end

% -----------------------------------------------------------------------
function [x_est_history, P_est_history] = UKF_GNSS(x0, P0, tspan, measurements, ...
                                               alpha, beta, kappa, R, ~, ~)
% UKF_GNSS: Implements an Unscented Kalman Filter for satellite navigation.
%
% PURPOSE:
%   Performs state estimation for a spacecraft using GNSS position updates.
%   It employs J2-perturbed dynamics for the prediction step and a standard
%   UKF update for Cartesian ECI measurements.
%
% INPUTS:
%   x0           - [6x1] Initial state vector [position; velocity]
%   P0           - [6x6] Initial state covariance matrix
%   tspan        - [Nx1] Time grid for the estimation process
%   measurements - [Nx3] GNSS Cartesian position measurements at each epoch
%   alpha,beta,k - [1x1] UKF tuning parameters (scaling, distribution, secondary)
%   R            - [3x3] Measurement noise covariance matrix
%
% OUTPUTS:
%   x_est_history   - [6xN] Filtered state history
%   P_est_history   - [6x6xN] Covariance history
%
% Step #7.1: Initialization of UKF tuning weights and time grid
    n = numel(x0);
    N = length(tspan);
    lambda = alpha^2*(n+kappa) - n;
    mu_prop = cspice_bodvrd('EARTH','GM',1); mu_prop = mu_prop(1);

    Wm = [lambda/(n+lambda); repmat(1/(2*(n+lambda)), 2*n, 1)];
    Wc = Wm;
    Wc(1) = Wc(1) + (1 - alpha^2 + beta);

    x_est_history = zeros(n, N);
    P_est_history = zeros(n, n, N);
    xk_plus = x0(:);
    Pk_plus = P0;
    x_est_history(:,1)   = xk_plus;
    P_est_history(:,:,1) = Pk_plus;

    optsODE = odeset('RelTol',1e-13,'AbsTol',[1e-11*ones(3,1);1e-12*ones(3,1)],'MaxStep',10);

    % Step #7.2: Begin recursive UKF loop (Prediction + Update)
    for j = 2:N
        % -- Sigma points generation --
        Pk_plus = 0.5*(Pk_plus + Pk_plus');
        try
            Sk = chol((n+lambda)*Pk_plus,'lower');
        catch
            [U,S] = eig(0.5*((n+lambda)*Pk_plus+((n+lambda)*Pk_plus)'));
            S = max(S,0);  Sk = U*sqrt(S);
        end
        Xi = [xk_plus, xk_plus+Sk, xk_plus-Sk];

        % -- Sigma points propagation (Non-linear J2 dynamics) --
        Re = cspice_bodvrd('EARTH','RADII',3); Re = Re(1);
        Xk = zeros(n, 2*n+1);
        for i = 1:2*n+1
            [~,xx] = ode113(@(t,x) keplerianJ2_RHS(t,x,mu_prop,Re), ...
                            [tspan(j-1), tspan(j)], Xi(:,i), optsODE);
            Xk(:,i) = xx(end,:)';
        end

        % -- Reconstruct a-priori mean and covariance --
        xk_minus = Xk * Wm;
        Pk_minus = zeros(n);
        for i = 1:2*n+1
            d = Xk(:,i) - xk_minus;
            Pk_minus = Pk_minus + Wc(i)*(d*d');
        end
        Pk_minus = 0.5*(Pk_minus + Pk_minus');

        % -- Measurement Update: Map sigma points to Cartesian ECI --
        yk = measurements(j,:)';   
        Y_sigma_points = Xk(1:3,:);
        y_minus = Y_sigma_points * Wm;
        Pyy = zeros(3);  Pxy = zeros(n,3);
        for i = 1:2*n+1
            dx = Xk(:,i)   - xk_minus;
            dy = Y_sigma_points(:,i) - y_minus;
            Pyy = Pyy + Wc(i)*(dy*dy');
            Pxy = Pxy + Wc(i)*(dx*dy');
        end

        % -- Compute Kalman Gain and perform linear state correction --
        Pyy = Pyy + R;
        K = Pxy / Pyy;
        xk_plus = xk_minus + K*(yk - y_minus);
        Pk_plus = Pk_minus - K*Pyy*K';
        Pk_plus = 0.5*(Pk_plus + Pk_plus');

        x_est_history(:,j)   = xk_plus;
        P_est_history(:,:,j) = Pk_plus;
    end
end

% -----------------------------------------------------------------------
function [x_est_history, P_est_history] = UKF_AIIP(x0, P0, tspan, y_aip, ai_visible, ...
                                               alpha, beta, kappa, R, ~, ~)
% UKF_AIIP: Unscented Kalman Filter specialized for Geodetic/AIP sensors.
%
% PURPOSE:
%   Estimates spacecraft state using geodetic measurements (Range, Longitude, 
%   Latitude). Handles the non-linear transformation to the rotating Earth 
%   frame (IAU_EARTH) and applies circular statistics to mitigate 
%   discontinuities in longitude (wrap-around at the date line).
%
% INPUTS:
%   x0           - [6x1] Initial state vector [position; velocity]
%   P0           - [6x6] Initial state covariance matrix
%   tspan        - [Nx1] Time grid for the estimation process
%   y_aip   - [3xN] Measured Geodetic observables: [Range; Lon; Lat]
%   ai_visible   - [1xN] Boolean flag indicating sensor line-of-sight availability
%   alpha,beta,k - [1x1] UKF tuning parameters
%   R            - [3x3] Measurement noise covariance matrix
%
% OUTPUTS:
%   x_est_history   - [6xN] Filtered state history
%   P_est_history   - [6x6xN] Covariance history

% Step #8.1: Initialization of UKF and dynamical parameters
    n = numel(x0);
    N = length(tspan);
    lambda = alpha^2*(n+kappa) - n;
    mu_prop = cspice_bodvrd('EARTH','GM',1); mu_prop = mu_prop(1);

    Wm = [lambda/(n+lambda); repmat(1/(2*(n+lambda)), 2*n, 1)];
    Wc = Wm;
    Wc(1) = Wc(1) + (1 - alpha^2 + beta);

    x_est_history = zeros(n, N);
    P_est_history = zeros(n, n, N);
    xk_plus = x0(:);
    Pk_plus = P0;
    x_est_history(:,1)   = xk_plus;
    P_est_history(:,:,1) = Pk_plus;

    optsODE = odeset('RelTol',1e-13,'AbsTol',[1e-11*ones(3,1);1e-12*ones(3,1)],'MaxStep',10);
    Re = cspice_bodvrd('EARTH','RADII',3); Re = Re(1);

    % Step #8.2: Start UKF recursive loop
    for j = 2:N
        % -- Sigma Point generation and Propagation --
        Pk_plus = 0.5*(Pk_plus + Pk_plus');
        try
            Sk = chol((n+lambda)*Pk_plus,'lower');
        catch
            [U,S] = eig(0.5*((n+lambda)*Pk_plus+((n+lambda)*Pk_plus)'));
            S = max(S,0);  Sk = U*sqrt(S);
        end
        Xi = [xk_plus, xk_plus+Sk, xk_plus-Sk];

        Xk = zeros(n, 2*n+1);
        for i = 1:2*n+1
            [~,xx] = ode113(@(t,x) keplerianJ2_RHS(t,x,mu_prop,Re), ...
                            [tspan(j-1), tspan(j)], Xi(:,i), optsODE);
            Xk(:,i) = xx(end,:)';
        end

        % -- A-priori Mean and Covariance Reconstruction --
        xk_minus = Xk * Wm;
        Pk_minus = zeros(n);
        for i = 1:2*n+1
            d = Xk(:,i) - xk_minus;
            Pk_minus = Pk_minus + Wc(i)*(d*d');
        end
        Pk_minus = 0.5*(Pk_minus + Pk_minus');

        % -- Step #8.3: Verification of measurement availability --
        if ~ai_visible(j) || any(~isfinite(y_aip(:,j)))
            xk_plus = xk_minus;
            Pk_plus = Pk_minus;
            x_est_history(:,j)   = xk_plus;
            P_est_history(:,:,j) = Pk_plus;
            continue;
        end

        % -- Step #8.4: Map Sigma Points to Geodetic Space (Range/Lon/Lat) --
        Y_sigma_points = zeros(3, 2*n+1);
        R_rot = cspice_pxform('J2000','IAU_EARTH', tspan(j));
        for i = 1:2*n+1
            r_ECEF = R_rot * Xk(1:3,i);
            [rho_i, lon_i, lat_i] = cspice_reclat(r_ECEF);
            Y_sigma_points(:,i) = [rho_i; wrapToPi(lon_i); lat_i];
        end

        % -- Step #8.5: Compute circular mean for Longitude to avoid averaging errors --
        y_minus = zeros(3,1);
        y_minus(1) = Y_sigma_points(1,:) * Wm;       
        y_minus(3) = Y_sigma_points(3,:) * Wm;       
        y_minus(2) = angle(sum(Wm(:).' .* exp(1i * Y_sigma_points(2,:))));

        % -- Step #8.6: Compute Covariances using angdiff for longitude residuals --
        Pyy = zeros(3);  Pxy = zeros(n,3);
        for i = 1:2*n+1
            dy = Y_sigma_points(:,i) - y_minus;
            dy(2) = angdiff(y_minus(2), Y_sigma_points(2,i));  
            dx = Xk(:,i) - xk_minus;
            Pyy = Pyy + Wc(i)*(dy*dy');
            Pxy = Pxy + Wc(i)*(dx*dy');
        end
        Pyy = Pyy + R;

        % -- Step #8.7: Compute innovationation and apply final correction --
        y = y_aip(:,j);
        y(2) = wrapToPi(y(2));
        innovation = y - y_minus;
        innovation(2) = angdiff(y_minus(2), y(2));

        K = Pxy / Pyy;
        xk_plus = xk_minus + K * innovation;
        Pk_plus = Pk_minus - K*Pyy*K';
        Pk_plus = 0.5*(Pk_plus + Pk_plus');

        x_est_history(:,j)   = xk_plus;
        P_est_history(:,:,j) = Pk_plus;
    end
end