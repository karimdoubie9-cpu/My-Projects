%% Spacecraft Guidance and Navigation (2025/2026)
% Assignment 2 – Exercise 1: Uncertainty Propagation for INTEGRAL
% Author: Karim Ahmed Mohamed Doubie
%
% Objective:
%   Point 1 – Apsis detection: locate pericenter/apocenter epochs over 5 revolutions.
%   Point 2 – Uncertainty propagation: compare LinCov (STM), Unscented Transform,
%              and Monte Carlo (2500 samples) for mean and covariance evolution.
%   Point 3 – MC convergence: identify Nmin for 1e-3 stability of the running mean.
%   Point 4 – Visualization: perifocal 1-sigma ellipses and worst-case 3-sigma
%              growth (principal eigenvalues) for position and velocity.


%% STEP #0: INITIAL SETTING AND INITIAL DEFINITION
clc; clearvars; close all;
addpath('kernels/');
% Step #0.1: Configure global LaTeX interpreter settings for professional plotting
set(groot, 'DefaultAxesFontSize',                  16);
set(groot, 'DefaultAxesTitleFontSizeMultiplier',    1);
set(groot, 'DefaultTextFontSize',                  12);
set(0,     'DefaultTextInterpreter',         'latex');
set(0,     'DefaultLegendInterpreter',       'latex');
set(0,     'DefaultAxesTickLabelInterpreter','latex');

% Step #0.2: Load SPICE kernels and define initial temporal/orbital parameters
cspice_furnsh('assignment02.tm');

t0_utc = '2025-10-26T01:10:57.769';
et0    = cspice_str2et(t0_utc);

% Initial state in ECI [km, km/s]
x0 = [7896.37847595; -3152.23575783; 149961.25963557;
      -0.6035374932;  0.0159618300;   0.0321154658];

% Step #0.3: Define initial position and velocity covariances
Pr0 = [ 6.7e-3, -2.1e-3,  7.6e-4;
       -2.1e-3,  9.7e-3,  5.3e-4;
        7.6e-4,  5.3e-4,  2.1e-3];
Pv0 = diag([4.7e-7, 3.1e-7, 1.6e-7]);
P0  = blkdiag(Pr0, Pv0);

mu = cspice_bodvrd('EARTH','GM',1);  mu = mu(1);

%% STEP #1: APSIS PASSAGE DETECTION AND NOMINAL TRAJECTORY
% Step #1.1: Compute orbital elements and period at t0
r0_n  = norm(x0(1:3));
v0_n  = norm(x0(4:6));
a_sma = 1 / (2/r0_n - v0_n^2/mu);          % semi-major axis [km]
T_orb = 2*pi * sqrt(a_sma^3 / mu);          % period [s]


% Create the table data
Variable = {'Orbital Period [h]'; 'Semi-major axis [km]'};
Value = [T_orb/3600; a_sma];

% Create the table
T = table(Value, 'RowNames', Variable);

% Display the table
fprintf('\n--- Orbital Parameters ---\n');
disp(T);


% Step #1.2: Propagate for 5.1 periods to detect apocentre/pericentre events
tf = et0 + 5.1*T_orb;
[~, ~, ~, xe, ~, te] = keplerian_STM_propagator(x0, et0, tf, mu);

% Step #1.3: Classify events based on radial distance relative to semi-major axis
r_ev   = vecnorm(xe(1:3,:));                 % 1 x Nevents
flag_apo = r_ev > a_sma;

te_peri = te(~flag_apo);   xe_peri = xe(:, ~flag_apo);
te_apo  = te( flag_apo);   xe_apo  = xe(:,  flag_apo);

% Step #1.4: Selection of the first 5 passage events for each apsis
np = min(5, numel(te_peri));  na = min(5, numel(te_apo));
te_peri = te_peri(1:np);   xe_peri = xe_peri(:, 1:np);
te_apo  = te_apo(1:na);    xe_apo  = xe_apo(:,  1:na);


% --- : Propagated mean state at pericentre passages ---
    fprintf('\n Propagated mean state at pericentre passages\n');   
    % 1. Print the Date Row
    fprintf('%-12s', 'Epoch [UTC]');
    for k = 1:np
        utc_str = cspice_et2utc(te_peri(k), 'ISOC', 3);
        fprintf('  %12s', utc_str(1:10)); % Date
    end       
    % 2. Print the Time Row
    fprintf('\n%-12s', '[UTC]');
    for k = 1:np
        utc_str = cspice_et2utc(te_peri(k), 'ISOC', 3);
        fprintf('  %12s', utc_str(12:19)); % Time
    end     
    % 3. Print the dt [h] Row
    fprintf('\n%-12s', 'dt [h]');
    for k = 1:np
        % Calculate delta time in hours from the start epoch et0
        delta_t_h = (te_peri(k) - et0) / 3600;
        fprintf('  %12.4f', delta_t_h);
    end
    fprintf('\n--------------------------------------------------------------------------\n');  
    % 4. Print the State components
    rowNames = {'rx [km]', 'ry [km]', 'rz [km]', 'vx [km/s]', 'vy [km/s]', 'vz [km/s]'};
    for i = 1:6
        fprintf('%-12s', rowNames{i});
        for k = 1:np
            fprintf('  %12.4f', xe_peri(i,k));
        end
        fprintf('\n');
    end
    fprintf('\n');

% --- Propagated mean state at apocentre passages ---
    fprintf('\n Propagated mean state at apocentre passages\n');
    % 1. Print the Date Row
    fprintf('%-12s', 'Epoch [UTC]');
    for k = 1:na
        utc_str = cspice_et2utc(te_apo(k), 'ISOC', 3);
        fprintf('  %12s', utc_str(1:10));
    end
    % 2. Print the Time Row
    fprintf('\n%-12s', '[UTC]');
    for k = 1:na
        utc_str = cspice_et2utc(te_apo(k), 'ISOC', 3);
        fprintf('  %12s', utc_str(12:19));
    end
    % 3. Print the dt [h] Row
    fprintf('\n%-12s', 'dt [h]');
    for k = 1:na
        % Calculate delta time in hours from the start epoch et0 for apocentres
        delta_t_h_apo = (te_apo(k) - et0) / 3600;
        fprintf('  %12.4f', delta_t_h_apo);
    end
    
    fprintf('\n--------------------------------------------------------------------------\n');
    % 4. Print the State components
    for i = 1:6
        fprintf('%-12s', rowNames{i});
        for k = 1:na
            fprintf('  %12.4f', xe_apo(i,k));
        end
        fprintf('\n');
    end
    fprintf('\n');

% Step #1.5: Build a chronological time grid encompassing all detected events
te_all      = sort([te_peri; te_apo]);
t_span_prop = [et0; te_all];          % (1 + np + na) x 1

% Recover integer indices of peri / apo inside t_span_prop
idx_peri = zeros(1, np);
idx_apo  = zeros(1, na);
for k = 1:np
    [~, idx_peri(k)] = min(abs(t_span_prop - te_peri(k)));
end
for k = 1:na
    [~, idx_apo(k)]  = min(abs(t_span_prop - te_apo(k)));
end

% --- Figure 1: 3-D nominal orbit with apsis markers ---
    [~, ~, xx_nom] = keplerian_propagator(x0, linspace(et0, tf, 3000), mu);
    figure('Name',' Keplerian Orbit');
    % --- Add Earth Sphere ---
    % Earth radius in km
    R_earth = 6378.137; 
    
    % Create sphere mesh
    [sx, sy, sz] = sphere(50);
    
    % Check if the image exists
    if exist('Earth_IMAGE.jpg','file')
        % Image found: map the texture
        earth_img = imread('Earth_IMAGE.jpg');
        surf(sx*R_earth, sy*R_earth, sz*R_earth, ...
            'FaceColor', 'texturemap', 'CData', earth_img, ...
            'EdgeColor', 'none', 'FaceAlpha', 1.0);
    else
        % Image not found: plot a plain sphere of the same dimension
        surf(sx*R_earth, sy*R_earth, sz*R_earth, ...
            'FaceColor', [0, 0.4470, 0.7410], ... % Default MATLAB blue
            'EdgeColor', 'none', 'FaceAlpha', 1.0);
    end
    hold on;
    
    % --- Plot Orbit and Apsis markers ---
    plot3(xx_nom(:,1), xx_nom(:,2), xx_nom(:,3), 'k--', 'LineWidth', 1); 
    plot3(xe_apo(1,:),  xe_apo(2,:),  xe_apo(3,:),  'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', 'Apocenters');
    plot3(xe_peri(1,:), xe_peri(2,:), xe_peri(3,:), 'co', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'c', 'DisplayName', 'Pericenters');
    
    % --- Formatting ---
    axis equal; % Essential to keep Earth spherical
    grid on;
    xlabel('$x$ [km]', 'Interpreter', 'latex');  
    ylabel('$y$ [km]', 'Interpreter', 'latex');  
    zlabel('$z$ [km]', 'Interpreter', 'latex');
    title('Spacecraft trajectory in ECI J2000  Keplerian motion');
    legend({'Earth', 'Orbit','Apocenters','Pericenters'},'Location','best');
    view(3); % Ensure 3D view
%% STEP #2: LINEARIZED COVARIANCE AND UNSCENTED TRANSFORM PROPAGATION
fprintf('\n--- Point 2: Propagating mean and covariance ---\n');

% Step #2.1: Propagation using Linearized Covariance (LinCov)
[x_mean_LC, P_all_LC] = LinCov(x0, P0, t_span_prop, mu);
fprintf('  LinCov complete.\n');

% Step #2.2: Propagation using Unscented Transform (UT) with specific tuning
alpha = 1e-2;  beta = 2;  kappa = 0;  N_st = 6;
[x_mean_UT, P_all_UT] = Unscented_transformation(x0, P0, t_span_prop, mu, ...
                                                   alpha, beta, kappa, N_st);
fprintf('  UT complete.\n');

% Step #2.3: LinCov vs UT mean differences at each apsis  
fprintf('\n--- LinCov vs UT mean differences ---\n');
for k = 1:np
    diff_pos = norm(x_mean_LC(1:3, idx_peri(k)) - x_mean_UT(1:3, idx_peri(k)));
    diff_vel = norm(x_mean_LC(4:6, idx_peri(k)) - x_mean_UT(4:6, idx_peri(k)));
    fprintf('Peri %d: ||LC-UT||_r = %.4f km   ||LC-UT||_v = %.2e km/s\n', ...
            k, diff_pos, diff_vel);
end
for k = 1:na
    diff_pos = norm(x_mean_LC(1:3, idx_apo(k)) - x_mean_UT(1:3, idx_apo(k)));
    diff_vel = norm(x_mean_LC(4:6, idx_apo(k)) - x_mean_UT(4:6, idx_apo(k)));
    fprintf('Apo  %d: ||LC-UT||_r = %.4f km   ||LC-UT||_v = %.2e km/s\n', ...
            k, diff_pos, diff_vel);
end

%% STEP #3: MONTE CARLO SIMULATION AND STATISTICAL VALIDATION
fprintf('\n--- Point 3: Monte Carlo ---\n');

% Step #3.1: Generate and propagate samples
rng(33);          % fixed seed for reproducibility
N_mc = 2500;     % number of samples

fprintf('  Propagating %d samples…\n', N_mc);
[x_mean_MC, P_all_MC, M_samples] = MC(x0, P0, t_span_prop, mu, N_mc);
fprintf('  MC complete.\n');

% Step #3.2: N_min Convergence Analysis at the final detected passage events
tol_rel = 1e-3;
N_vec   = 1:N_mc;

% Sample clouds at last apsis
samp_apo  = squeeze(M_samples(:, :, idx_apo(end)));    % 6 x N_mc
samp_peri = squeeze(M_samples(:, :, idx_peri(end)));   % 6 x N_mc

% Step #3.3: Compute running means and relative error against N_mc reference
run_mean_apo  = cumsum(samp_apo,  2) ./ N_vec;   % 6 x N_mc
run_mean_peri = cumsum(samp_peri, 2) ./ N_vec;

% Reference = full-N_mc sample mean
ref_apo  = run_mean_apo(:,  end);
ref_peri = run_mean_peri(:, end);

% Component-wise relative error
err_apo  = abs(run_mean_apo  - ref_apo)  ./ (abs(ref_apo)  + eps);
err_peri = abs(run_mean_peri - ref_peri) ./ (abs(ref_peri) + eps);

% Determine smallest N where all components satisfy tolerance
Nmin_Apogee = find(all(err_apo  < tol_rel, 1), 1, 'first');
Nmin_Perigee = find(all(err_peri < tol_rel, 1), 1, 'first');
if isempty(Nmin_Apogee), Nmin_Apogee = N_mc; end
if isempty(Nmin_Perigee), Nmin_Perigee = N_mc; end
N_min = max(Nmin_Apogee, Nmin_Perigee);

% Define the data
Variable = {'tol_rel'; 'N_min,Apogee'; 'N_min,Perigee'; 'N_min,MC'};
Value = [tol_rel; Nmin_Apogee; Nmin_Perigee; N_min];
% Create the table
T = table(Variable, Value);
% Display the table
fprintf('\n--- Parameters Table ---\n');
disp(T);

% --- Figure – Running MC mean of STATE components ---
figure('Name','MC running mean – state');
% Combine labels for all 6 state components
lbl_state = {'$r_x$ [km]', '$r_y$ [km]', '$r_z$ [km]', ...
             '$v_x$ [km/s]', '$v_y$ [km/s]', '$v_z$ [km/s]'};
% Define colors: Green and Purple
col_apo =[0 0.5 0]; 
col_peri = [0.5 0 0.5];
% Loop through all 6 components (1-3 for position, 4-6 for velocity)
for i = 1:6
    subplot(3,2,i);  
    hold on;  grid on;
    % Plot Apo and Peri (indices naturally match 1:6 for run_mean matrices)
    plot(N_vec, run_mean_apo(i,:),  '-',  'Color', col_apo,  'LineWidth',1.2, 'DisplayName','Apo');
    plot(N_vec, run_mean_peri(i,:), '--', 'Color', col_peri, 'LineWidth',1.2, 'DisplayName','Peri');
    xline(N_min, 'k:', 'LineWidth',1.2, 'HandleVisibility','off');
    % Set labels and limits
    ylabel(lbl_state{i}, 'Interpreter', 'latex');  
    xlabel('$N_{\mathrm{MC}}$', 'Interpreter', 'latex'); % Moved inside loop so it appears on all subplots 
    xlim([1 N_mc]);    
    % Add legend only to the first subplot
    if i == 1
        legend('Interpreter', 'latex', 'Location','best');
    end
end
% Set a single overarching title
sgtitle('Estimated mean state at the last apocentre and pericentre (MC)');

% Step #3.4: Construct Perifocal frame for visualization of in-plane uncertainty
ecc_v = ((v0_n^2 - mu/r0_n)*x0(1:3) - dot(x0(1:3),x0(4:6))*x0(4:6)) / mu;
P_hat = ecc_v / norm(ecc_v);                        % pericenter direction
h_v   = cross(x0(1:3), x0(4:6));
W_hat = h_v / norm(h_v);                            % orbital normal
Q_hat = cross(W_hat, P_hat);                        % in-plane, perp to P
R_pf  = [P_hat, Q_hat]';   % 2 x 3  –  projects ECI pos → [P; Q]



% --- Figure 4 – Orbital plane at ALL APOCENTRES ---
c_MC = [0 0.45 0.74];   c_LC = [0.85 0.33 0.10];   c_UT = [0.47 0.67 0.19];
% % --- Figure 6 – Time evolution of 3σ_r and 3σ_v ---
[s3r_LC, s3v_LC] = sigma3_Kepler(P_all_LC, t_span_prop, 'max');
[s3r_UT, s3v_UT] = sigma3_Kepler(P_all_UT, t_span_prop, 'max');
[s3r_MC, s3v_MC] = sigma3_Kepler(P_all_MC, t_span_prop, 'max');

dt_h_all  = (t_span_prop - et0) / 3600;    % Δt from t0 [h]
dt_apo    = dt_h_all(idx_apo);
dt_peri   = dt_h_all(idx_peri);
LW = 1.5;

figure('Name',' 3sigma evolution')

subplot(2,2,1);  hold on;  grid on;
plot(dt_apo, s3r_LC(idx_apo), '-o','Color',c_LC,'LineWidth',LW,'MarkerSize',5,'DisplayName','LinCov');
plot(dt_apo, s3r_UT(idx_apo), '-o','Color',c_UT,'LineWidth',LW,'MarkerSize',5,'DisplayName','UT');
plot(dt_apo, s3r_MC(idx_apo), '-*','Color',c_MC,'LineWidth',LW,'MarkerSize',7,'DisplayName','MC');
ylabel('$3\sigma_r$ [km]');  title('Apocenters');
xlabel('$t - t_0$ [h]');
legend('Location','northwest');

subplot(2,2,2);  hold on;  grid on;
plot(dt_peri, s3r_LC(idx_peri), '-o','Color',c_LC,'LineWidth',LW,'MarkerSize',5);
plot(dt_peri, s3r_UT(idx_peri), '-o','Color',c_UT,'LineWidth',LW,'MarkerSize',5);
plot(dt_peri, s3r_MC(idx_peri), '-*','Color',c_MC,'LineWidth',LW,'MarkerSize',7);
ylabel('$3\sigma_r$ [km]');  title('Pericenters');
xlabel('$t - t_0$ [h]');

subplot(2,2,3);  hold on;  grid on;
plot(dt_apo, s3v_LC(idx_apo), '-o','Color',c_LC,'LineWidth',LW,'MarkerSize',5);
plot(dt_apo, s3v_UT(idx_apo), '-o','Color',c_UT,'LineWidth',LW,'MarkerSize',5);
plot(dt_apo, s3v_MC(idx_apo), '-*','Color',c_MC,'LineWidth',LW,'MarkerSize',7);
ylabel('$3\sigma_v$ [km/s]');  xlabel('$t - t_0$ [h]');

subplot(2,2,4);  hold on;  grid on;
plot(dt_peri, s3v_LC(idx_peri), '-o','Color',c_LC,'LineWidth',LW,'MarkerSize',5);
plot(dt_peri, s3v_UT(idx_peri), '-o','Color',c_UT,'LineWidth',LW,'MarkerSize',5);
plot(dt_peri, s3v_MC(idx_peri), '-*','Color',c_MC,'LineWidth',LW,'MarkerSize',7);
ylabel('$3\sigma_v$ [km/s]');  xlabel('$t - t_0$ [h]');

sgtitle(['$3\sqrt{\max(\lambda_i(P_{r,i}))}$ and ' ...
         '$3\sqrt{\max(\lambda_i(P_{v,i}))}$ for MC, LinCov, UT']);

% --- Figure  – Orbital plane at ALL APOCENTRES ---
c_MC = [0 0.45 0.74];
c_LC = [0.85 0.33 0.10];
c_UT = [0.47 0.67 0.19];

% Distinct colours for MC sample clouds – one per apsis event
cmap_apo = lines(na);   

figure('Name',' Orbital plane – apocenters');
hold on;  grid on;  box on;

% Scatter MC samples – each apocentre in its own colour
h_smp_apo = gobjects(na,1);
for k = 1:na
    smp = R_pf * squeeze(M_samples(1:3, :, idx_apo(k)));   % 2 x N_mc
    h_smp_apo(k) = scatter(smp(1,:), smp(2,:), 1.5, ...
                            cmap_apo(k,:), '.', ...
                            'DisplayName', sprintf('MC samples – apo %d', k));
end

% 1σ covariance ellipses + mean markers for each apocentre
h_ell_apo = [];
h_mLC_apo = [];  h_mUT_apo = [];  h_mMC_apo = [];
for k = 1:na
    ki  = idx_apo(k);
    mMC = R_pf * x_mean_MC(1:3, ki);
    mLC = R_pf * x_mean_LC(1:3, ki);
    mUT = R_pf * x_mean_UT(1:3, ki);
    PMC = R_pf * P_all_MC(1:3,1:3,ki) * R_pf';
    PLC = R_pf * P_all_LC(1:3,1:3,ki) * R_pf';
    PUT = R_pf * P_all_UT(1:3,1:3,ki) * R_pf';

    hh = plotcov1sig(PMC, mMC', 'Color', c_MC, 'LineWidth', 1.2);
    if k == 1, h_ell_apo = hh; end
    plotcov1sig(PLC, mLC', 'Color', c_LC, 'LineWidth', 1.2);
    plotcov1sig(PUT, mUT', 'Color', c_UT, 'LineWidth', 1.2);

    hh = plot(mLC(1), mLC(2), 'o', 'Color', c_LC, ...
              'MarkerSize', 5, 'MarkerFaceColor', 'w');
    if k == 1, h_mLC_apo = hh; end

    hh = plot(mUT(1), mUT(2), 'o', 'Color', c_UT, ...
              'MarkerSize', 5, 'MarkerFaceColor', 'w');
    if k == 1, h_mUT_apo = hh; end

    hh = plot(mMC(1), mMC(2), 'o', 'Color', c_MC, ...
              'MarkerSize', 5, 'MarkerFaceColor', 'w');
    if k == 1, h_mMC_apo = hh; end
end


xlabel('$x_p$ [km]', 'Interpreter', 'latex');
ylabel('$y_p$ [km]', 'Interpreter', 'latex');
title('Uncertainty distribution projected onto the orbital plane at apocentre passages');

% Legend order: LinCov mean | UT mean | MC mean | MC samples | 1sigma ellipses
legend([h_mLC_apo, h_mUT_apo, h_mMC_apo, h_smp_apo(1), h_ell_apo], ...
       {'LinCov mean', 'UT mean', 'MC mean', ...
        'MC samples', '$1\sigma$ ellipses (MC/LinCov/UT)'}, ...
       'Location', 'best');

% --- Figure 5 – Orbital plane at ALL PERICENTRES ---
cmap_peri = lines(np);  

figure('Name',' Orbital plane – pericenters');
hold on;  grid on;  box on;

h_smp_peri = gobjects(np,1);
for k = 1:np
    smp = R_pf * squeeze(M_samples(1:3, :, idx_peri(k)));
    h_smp_peri(k) = scatter(smp(1,:), smp(2,:), 1.5, ...
                             cmap_peri(k,:), '.', ...
                             'DisplayName', sprintf('MC samples – peri %d', k));
end

h_ell_peri = [];
h_mLC_peri = [];  h_mUT_peri = [];  h_mMC_peri = [];
for k = 1:np
    ki  = idx_peri(k);
    mMC = R_pf * x_mean_MC(1:3, ki);
    mLC = R_pf * x_mean_LC(1:3, ki);
    mUT = R_pf * x_mean_UT(1:3, ki);
    PMC = R_pf * P_all_MC(1:3,1:3,ki) * R_pf';
    PLC = R_pf * P_all_LC(1:3,1:3,ki) * R_pf';
    PUT = R_pf * P_all_UT(1:3,1:3,ki) * R_pf';

    hh = plotcov1sig(PMC, mMC', 'Color', c_MC, 'LineWidth', 1.2);
    if k == 1, h_ell_peri = hh; end
    plotcov1sig(PLC, mLC', 'Color', c_LC, 'LineWidth', 1.2);
    plotcov1sig(PUT, mUT', 'Color', c_UT, 'LineWidth', 1.2);

    hh = plot(mLC(1), mLC(2), 'o', 'Color', c_LC, ...
              'MarkerSize', 5, 'MarkerFaceColor', 'w');
    if k == 1, h_mLC_peri = hh; end

    hh = plot(mUT(1), mUT(2), 'o', 'Color', c_UT, ...
              'MarkerSize', 5, 'MarkerFaceColor', 'w');
    if k == 1, h_mUT_peri = hh; end

    hh = plot(mMC(1), mMC(2), 'o', 'Color', c_MC, ...
              'MarkerSize', 5, 'MarkerFaceColor', 'w');
    if k == 1, h_mMC_peri = hh; end
end

xlabel('$x_p$ [km]', 'Interpreter', 'latex');
ylabel('$y_p$ [km]', 'Interpreter', 'latex');
title('Uncertainty distribution projected onto the orbital plane at pericentre passages');

legend([h_mLC_peri, h_mUT_peri, h_mMC_peri, h_smp_peri(1), h_ell_peri], ...
       {'LinCov mean', 'UT mean', 'MC mean', ...
        'MC samples', '$1\sigma$ ellipses (MC/LinCov/UT)'}, ...
       'Location', 'best');




%% ================================================================
%  HELPER FUNCTIONS
%% ================================================================

function h = plotcov1sig(P, mu_pt, varargin)
% PURPOSE: 
%   Generates a 2D geometric 1-sigma confidence ellipse based on a 
%   provided covariance matrix and mean state. Used for visualizing 
%   in-plane uncertainty distribution.
%
% INPUTS:
%   P       : [2x2] Covariance sub-matrix (e.g., position P-Q plane)
%   mu_pt   : [2x1] Mean state vector for centering
%   varargin: Plotting options (Color, LineWidth, etc.)
%
% OUTPUTS:
%   h       : Handle to the plotted line object

% Step #1.1: Generate points for a unit circle
angles      = linspace(0, 2*pi, 200);
unit_circle = [cos(angles); sin(angles)];

% Step #1.2: Perform Eigen-decomposition to find principal axes
[V, D]      = eig(P);

% Step #1.3: Scale and rotate the circle to match the 1-sigma distribution
ellipse     = V * sqrt(abs(D)) * unit_circle;

% Step #1.4: Translate to mean and plot
h = plot(ellipse(1,:) + mu_pt(1), ellipse(2,:) + mu_pt(2), varargin{:});
end

% ----------------------------------------------------------------

function [x_mean, P] = LinCov(x0, P0, t_span, mu)
% PURPOSE: 
%   Propagates state uncertainty using the Linearized Covariance (LinCov) 
%   method via the State Transition Matrix (STM). This assumes the 
%   system remains close to the nominal trajectory (linear regime).
%
% INPUTS:
%   x0      : [6x1] Initial state [km, km/s]
%   P0      : [6x6] Initial covariance [km^2, km^2/s^2]
%   t_span  : [Mx1] Epoch grid for evaluation [s]
%   mu      : [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUTS:
%   x_mean  : [6xM] Propagated nominal state
%   P       : [6x6xM] Propagated covariance tensors

% Step #2.1: Integrate nominal state and augmented STM components
sol    = keplerian_STM_propagator(x0, t_span(1), t_span(end), mu);
xx_aug = deval(sol, t_span);
x_mean = xx_aug(1:6, :);

% Step #2.2: Linearly map the covariance across the time grid
M      = length(t_span);
P      = zeros(6,6,M);
for i = 1:M
    % Extract PHI(t, t0) from augmented state
    PHI         = reshape(xx_aug(7:end, i), 6, 6);
    % Apply P(t) = PHI * P0 * PHI'
    P(:,:,i)    = PHI * P0 * PHI';
end
end

% ----------------------------------------------------------------

function [x_mean, P] = Unscented_transformation(x0, P0, t_span, mu, ...
                                                  alpha, beta, kappa, N)
% PURPOSE: 
%   Implements the Scaled Unscented Transform for uncertainty 
%   propagation. It captures non-linear effects by propagating a 
%   discrete set of sigma points through the full dynamics.
%
% INPUTS:
%   x0, P0  : Initial mean and covariance
%   t_span  : Evaluation epochs [s]
%   alpha   : Spread parameter (usually small, e.g., 1e-3)
%   beta    : Prior knowledge (2 is optimal for Gaussian)
%   kappa   : Secondary scaling parameter (usually 0)
%   N       : Dimension of the state (6)
%
% OUTPUTS:
%   x_mean  : [6xM] Propagated mean
%   P       : [6x6xM] Propagated covariance

% Step #3.1: Calculate scaling factor and generate weights
lambda = alpha^2 * (N + kappa) - N;
W_m        = ones(2*N+1,1) / (2*(N+lambda));
W_c        = W_m;
W_m(1)     = lambda / (N + lambda);
W_c(1)     = W_m(1) + (1 - alpha^2 + beta);

% Step #3.2: Generate 2N+1 sigma points via Matrix Square Root
try
    S = chol((N + lambda)*P0, 'lower');
catch
    S = sqrtm((N + lambda)*P0);
end
sigma_pts = [x0, x0 + S, x0 - S]; 

% Step #3.3: Propagate every sigma point individually
num_ep = length(t_span);
prop   = zeros(N, 2*N+1, num_ep);
for i = 1:(2*N+1)
    [~, ~, xx] = keplerian_propagator(sigma_pts(:,i), t_span, mu);
    prop(:,i,:) = xx';
end

% Step #3.4: Reconstruct mean and covariance from weighted samples
x_mean = zeros(N, num_ep);
P      = zeros(N, N, num_ep);
for t = 1:num_ep
    Y           = squeeze(prop(:,:,t));          % [N x 2N+1]
    x_mean(:,t) = Y * W_m;
    for j = 1:(2*N+1)
        d          = Y(:,j) - x_mean(:,t);
        P(:,:,t)   = P(:,:,t) + W_c(j)*(d*d');
    end
end
end

% ----------------------------------------------------------------

function [x_mean, P, propagated_samples] = MC(x0, P0, t_span, mu, n_mc)
% PURPOSE: 
%   Performs Monte Carlo uncertainty propagation by sampling the 
%   initial distribution and propagating each particle through the 
%   non-linear EOMs to provide a baseline for accuracy comparison.
%
% INPUTS:
%   x0, P0 : Initial state and covariance
%   t_span : Evaluation epochs [s]
%   n_mc   : Number of samples
%
% OUTPUTS:
%   x_mean             : [6xM] Sample mean
%   P                  : [6x6xM] Sample covariance
%   propagated_samples : [6xN_mc x M] Raw sample history

% Step #4.1: Generate initial random sample cloud
initial_samples = mvnrnd(x0', P0, n_mc)';    

% Step #4.2: Parallel/Serial propagation of the entire cloud
num_ep             = length(t_span);
propagated_samples = zeros(6, n_mc, num_ep);
for i = 1:n_mc
    [~, ~, xx] = keplerian_propagator(initial_samples(:,i), t_span, mu);
    propagated_samples(:,i,:) = xx';
end

% Step #4.3: Calculate sample statistics for every epoch
x_mean = zeros(6, num_ep);
P      = zeros(6, 6, num_ep);
for t = 1:num_ep
    cloud       = squeeze(propagated_samples(:,:,t));
    x_mean(:,t) = mean(cloud, 2);
    P(:,:,t)    = cov(cloud');
end
end

% ----------------------------------------------------------------

function [xf, tt, xx] = keplerian_propagator(x0, tspan, mu)
% PURPOSE: High-fidelity integration of the two-body problem equations 
% of motion.
%
% INPUTS:
%   x0    : [6x1] Initial state
%   tspan : [Mx1] Integration time span [s]
%   mu    : [1x1] Gravity [km^3/s^2]
%
% OUTPUTS:
%   xf    : [6x1] Final state
%   tt    : [Mx1] Output time grid
%   xx    : [Mx6] Output state history

% Step #5.1: Configure ODE solver and integrate dynamics
opts    = odeset('RelTol',1e-12,'AbsTol',1e-12);
[tt,xx] = ode113(@(t,x) keplerian_rhs(t,x,mu), tspan, x0, opts);
xf      = xx(end,:)';
end

% ----------------------------------------------------------------

function [sol, xf, PHIf, xe, PHIe, te] = keplerian_STM_propagator(x0, t0, tf, mu)
% PURPOSE: 
%   Integrates the spacecraft state along with the State Transition 
%   Matrix (STM). Includes an event function to capture apsis 
%   (pericenter/apocenter) crossings.
%
% INPUTS:
%   x0 : Initial state
%   t0 : Start time [s]
%   tf : End time [s]
%   mu : Gravity
%
% OUTPUTS:
%   sol  : ODE solution structure
%   PHIf : Final STM [6x6]
%   xe   : States at detected events
%   te   : Times at detected events

% Step #6.1: Augment state with Identity matrix for STM initialization
Phi0         = eye(6);
x0_aug       = [x0; Phi0(:)];

% Step #6.2: Integrate with Apsis event detection
opts         = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@detect_apo_peri);
sol          = ode113(@(t,x) keplerian_STM_rhs(t,x,mu), [t0,tf], x0_aug, opts);

% Step #6.3: Extract final state and STM
xx_f = deval(sol, tf);
xf   = xx_f(1:6);
PHIf = reshape(xx_f(7:end), 6, 6);

% Step #6.4: Extract event data
te  = sol.xe';        
ye  = sol.ye;         
xe  = ye(1:6,:);
num_ev = length(te);
PHIe   = zeros(6,6,num_ev);
for i = 1:num_ev
    PHIe(:,:,i) = reshape(ye(7:end,i), 6, 6);
end
end

% ----------------------------------------------------------------
% 
function dxdt = keplerian_rhs(~, x, mu)
% AIM: Computes acceleration for unperturbed Keplerian two-body dynamics.
% INPUTS: 
%   x  : [6x1] State vector [pos; vel] [km, km/s].
%   mu : [1x1] Gravitational parameter [km^3/s^2].
% OUTPUT: 
%   dxdt : [6x1] State derivative [vel; acc] [km/s, km/s^2].

    r=x(1:3);
    dxdt = [x(4:6); -mu*r/(sqrt(dot(r,r))*dot(r,r))];
end


% ----------------------------------------------------------------
% 
function dx = keplerian_STM_rhs(~, x, mu)
% PURPOSE: 
%   Calculates the derivative of the augmented state, including the 
%   Jacobian of the gravitational field for STM propagation.
%
% Step #8.1: Extract state and reshape current STM
r    = x(1:3);
rn   = norm(r);
Phi  = reshape(x(7:end), 6, 6);

% Step #8.2: Compute Gravity Gradient Matrix (G) and Jacobian (A)
A    = [zeros(3),                                                   eye(3);
         (3*mu/rn^5)*(r*r') - (mu/rn^3)*eye(3),                  zeros(3)];

% Step #8.3: Solve Variational Equation: dPhi/dt = A * Phi
dPhi = A * Phi;
dx   = zeros(42,1);
dx(1:6)   = [x(4:6); -mu*r/rn^3];
dx(7:end) = dPhi(:);
end



% ----------------------------------------------------------------


function [val, stop, dir] = detect_apo_peri(~, x)
% PURPOSE: Triggers at apsis passages where radial velocity vanishes.
% Condition: d|r|/dt = (r·v)/|r| = 0  ↔  r·v = 0
    val  = dot(x(1:3), x(4:6));   % dot via transpose-multiply
    stop = 0;                    % continue integration
    dir  = 0;                    % both minima and maxima
end
% ----------------------------------------------------------------

function [s3_pos, s3_vel] = sigma3_Kepler(P_all, T_grid, mode)
% PURPOSE: 
%   Extracts a scalar 3-sigma uncertainty metric from the 
%   covariance tensors by finding the principal eigenvalue (worst-case).
%
% INPUTS:
%   P_all  : [6x6xM] Covariance tensors
%   T_grid : Time epochs
%   mode   : 'max' for worst-case, 'min' for best-case
%
% OUTPUTS:
%   s3_pos : [1xM] 3-sigma position envelope [km]
%   s3_vel : [1xM] 3-sigma velocity envelope [km/s]

% Step #10.1: Calculate eigenvalues of position/velocity sub-matrices per epoch
M      = length(T_grid);
s3_pos = zeros(1,M);
s3_vel = zeros(1,M);
for i = 1:M
    er = eig(P_all(1:3,1:3,i));
    ev = eig(P_all(4:6,4:6,i));
    
    % Step #10.2: Extract square root of variance along principal axis
    if strcmp(mode,'max')
        s3_pos(i) = 3*sqrt(max(er));
        s3_vel(i) = 3*sqrt(max(ev));
    else
        s3_pos(i) = 3*sqrt(min(er));
        s3_vel(i) = 3*sqrt(min(ev));
    end
end
end