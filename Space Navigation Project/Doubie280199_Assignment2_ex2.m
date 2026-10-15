%% Spacecraft Guidance and Navigation (2025/2026)
% Assignment 2 – Exercise 2: Orbit Determination for INTEGRAL
% Author: Karim Ahmed Mohamed Doubie
%
% Objective:
%   Point 1 – Compute visibility windows for TROLL and KOUROU ground stations.
%   Point 2 – Generate SGP4 truth trajectory and synthesize noisy measurements
%              (Az, El, Range).
%   Point 3 – Batch least-squares OD using TROLL only (angles-only and
%              angles+range, Keplerian and J2): cases (i) to (iv).
%   Point 4 – Extended visibility: repeat OD including KOUROU measurements:
%              cases (v) and (vi).
%   Point 5 – Propagate all six estimated states over 6 orbital periods and
%              assess performance against SGP4 ground truth.


%% Step #0: Initialization and environment setup
clc; clearvars; close all;

addpath('kernels/');
addpath('sgp4/');


% Set the default font size for all axes (includes tick labels and titles)
set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1);
% Optional: Set the font size for all other text objects (legends, text annotations)
set(groot, 'DefaultTextFontSize', 12);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');

cspice_furnsh('assignment02.tm');

rng(55);
format long
 %  Given Data
    noise_az_deg  = 125e-3;%deg
    noise_el_deg  = 125e-3; %deg
    noise_rho_km  = 0.01;%km
    noise_az_rad  = noise_az_deg * pi/180;%rad
    noise_el_rad  = noise_el_deg * pi/180; %rad

%% Step#1: VISIBILITY WINDOWS
fprintf('--- Point 1: Visibility Windows ---\n\n');

% Step #1.1: Define time epochs and initial spacecraft state
utc_start = '2025-10-26T01:10:57.769';   epoch_start = cspice_str2et(utc_start);
utc_end   = '2025-10-28T10:00:00.000';   epoch_end   = cspice_str2et(utc_end);
utc_ref   = '2025-10-27T07:10:57.769';   epoch_ref   = cspice_str2et(utc_ref);

pos_init = [7896.37847595; -3152.23575783; 149961.25963557];
vel_init = [-0.6035374932;  0.0159618300;   0.0321154658];
state_init = [pos_init; vel_init];

% Step #1.2: Retrieve constants and propagate initial trajectory guesses
grav_param = cspice_bodvrd('EARTH','GM',1);   grav_param = grav_param(1);
earth_radius = cspice_bodvrd('EARTH','RADII',3); earth_radius = earth_radius(1);

J2_coeff = 0.0010826269;
tgrid_pre = epoch_start:60:epoch_ref;
[~, ~, traj_kep] = keplerian_propagator(state_init, tgrid_pre, 'EARTH', 'NONE');
[~, ~, traj_j2]  = keplerian_propagator(state_init, tgrid_pre, 'EARTH', 'J2');
guess_kep = traj_kep(end,:)';
guess_j2  = traj_j2(end,:)';

% Step #1.3: Propagate observation trajectory and perform visibility check
tgrid_obs = epoch_ref:60:epoch_end;
[~, epoch_col, traj_obs] = keplerian_propagator(guess_kep, tgrid_obs, 'EARTH', 'NONE');
epoch_row    = epoch_col.';
state_sat_eci = traj_obs.';

ground_stations = {'TROLL','KOUROU'};
size_GS= numel(ground_stations);
elev_mask_deg = 0;
elev_mask_rad = elev_mask_deg * pi/180;

vis_data = struct();
for station = 1:size_GS
    stn_name = ground_stations{station};
    % Compute antenna pointing vectors
    [azimuth, elevation, slant_range, ~] = antenna_pointing_rad(stn_name, epoch_row, state_sat_eci);
    above_mask = elevation > elev_mask_rad;
    idx_above = find(above_mask);

    % Step #1.4: Identify visibility windows
    if isempty(idx_above)
        win_start_idx = []; win_end_idx = [];
        fprintf('\nVisibility windows @ %s: none\n', stn_name);
    else
        d             = diff([0 above_mask 0]);
        win_start_idx = find(d ==  1);
        win_end_idx   = find(d == -1) - 1;
        num_windows = numel(win_start_idx);
        UTC_in  = arrayfun(@(i) cspice_et2utc(epoch_row(win_start_idx(i)),'ISOC',3), (1:num_windows)', 'UniformOutput', false);
        UTC_out = arrayfun(@(i) cspice_et2utc(epoch_row(win_end_idx(i)),'ISOC',3),   (1:num_windows)', 'UniformOutput', false);

        Tbl_vis = table((1:num_windows)', UTC_in, UTC_out, ...
            'VariableNames', {'Window', 'Start_UTC', 'End_UTC'});

        fprintf('\nVisibility windows at %s (elev > %.1f deg):\n', stn_name, elev_mask_deg);
        disp(Tbl_vis);
    end
    
    % Step #1.5: Store data and generate visualization plots
    vis_data.(stn_name).azimuth    = azimuth;
        az_in_deg    = azimuth   * 180/pi;
    vis_data.(stn_name).elevation  = elevation;
        el_in_deg    = elevation * 180/pi;
    vis_data.(stn_name).slant_rng  = slant_range;
        rho_km       = slant_range;
    vis_data.(stn_name).win_start  = win_start_idx;
    vis_data.(stn_name).win_end    = win_end_idx;
    vis_data.(stn_name).epochs     = epoch_row;
        t_from_ref   = (epoch_row - epoch_ref)/3600;

    clr_data = [0.20, 0.60, 0.20];
    clr_mask = [0.40, 0.00, 0.60];

    figure('Name',sprintf('Visibility_%s',stn_name));
    sgtitle(['Visibility Profile: ',stn_name]);
    
    subplot(3,1,1)
    plot(t_from_ref(above_mask), az_in_deg(above_mask), 'Color', clr_data, 'LineWidth', 2.2);
    grid on;
    xlabel('$t-t_i$ [h]'); ylabel('Azimuth [$^\circ$]');
    title('Azimuth Profile');
    
    subplot(3,1,2)
    plot(t_from_ref(above_mask), el_in_deg(above_mask), 'Color', clr_data, 'LineWidth', 2.2); hold on
    yline(elev_mask_deg,'--', 'Color', clr_mask, 'LineWidth', 2.0, 'DisplayName','Elevation Mask');
    grid on;
    xlabel('$t-t_i$ [h]'); ylabel('Elevation [$^\circ$]');
    title('Elevation Profile');
    
    subplot(3,1,3)
    plot(t_from_ref(above_mask), rho_km(above_mask), 'Color', clr_data, 'LineWidth', 2.2);
    grid on;
    xlabel('$t-t_i$ [h]'); ylabel('Range [km]');
    title('Range Profile');
    
    % --- Figure: Ground Track / Sky Map (Azimuth vs. Elevation) ---
    figure('Name', sprintf('Azimuth_vs_Elevation_%s', stn_name));
    plot(az_in_deg(above_mask), el_in_deg(above_mask), 'Color', clr_data, 'LineWidth', 2.2); hold on
    yline(elev_mask_deg, '--', 'Color', clr_mask, 'LineWidth', 2.0, 'DisplayName', 'Elevation Mask');
    grid on;
    xlabel('Azimuth [$^\circ$]'); ylabel('Elevation [$^\circ$]');
    title(['Az-El profile over TROLL station: ', stn_name]);
    axis tight;
end

%%  Step# 2: SGP4 TRUTH + NOISY MEASUREMENTS
fprintf('\n===================== POINT 2: SGP4 TRUTH + NOISY MEASUREMENTS =====================\n');

% Step #2.1: Initialize SGP4 propagator and reference epoch state
run_type='u'; ops_mode='a'; wgs_const=72;
dt_betweenum_obs = 60.0;
tle_line1 = '1 27540U 02048A   25298.18996762  .00000023  00000-0  00000-0 0  9991';
tle_line2 = '2 27540  91.0900 358.3898 8628526 267.2169  62.6507  0.37900391 25750';
sat_rec = twoline2rv(tle_line1, tle_line2, run_type, 'e', ops_mode, wgs_const);

[yr,mo,dy,hr,mn,sc] = invjday(sat_rec.jdsatepoch, sat_rec.jdsatepochf);
tle_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f',[yr,mo,dy,hr,mn,sc]);
tle_epoch_et  = cspice_str2et(tle_epoch_str);


Tbl_noise = table(noise_az_deg, noise_el_deg, noise_rho_km, ...
    'VariableNames', {'σ_Az [deg]', 'σ_El [deg]', 'σ_ρ [km]'});
disp(Tbl_noise);

Tbl_sat = table(sat_rec.satnum, {tle_epoch_str}, {utc_ref}, ...
    'VariableNames', {'Sat_ID', 'TLE_Reference_UTC', 'ti_Reference_UTC'});
disp(Tbl_sat);

% Step #2.2: Define frame conversion constants and reference truth state
nutation_dpsi = -0.118679 * pi/(180*3600);
nutation_deps = -0.009495 * pi/(180*3600);
tsince_ref = (epoch_ref - tle_epoch_et)/60;
[~, r_ref_teme, v_ref_teme] = sgp4(sat_rec, tsince_ref);
years_ref = cspice_unitim(epoch_ref, 'ET', 'TDT') / cspice_jyear();
time_ref   = years_ref / 100;
accel_teme = [0;0;0];
[r_ref_eci, v_ref_eci, ~] = teme2eci(r_ref_teme(:), v_ref_teme(:), accel_teme, time_ref, nutation_dpsi, nutation_deps);
x_truth_ref = [r_ref_eci; v_ref_eci];

% Step #2.3: Generate noisy measurements for each ground station
meas_store = struct();
for station = 1:numel(ground_stations)
   stn_name   = ground_stations{station};
    w_start    = vis_data.(stn_name).win_start;
    w_end      = vis_data.(stn_name).win_end;
    epoch_grid = vis_data.(stn_name).epochs;
    
    et_obs = [];
    jWin = 1;
    while jWin <= numel(w_start)
        t_in_win  = epoch_grid(w_start(jWin));
        t_out_win = epoch_grid(w_end(jWin));
        et_obs = [et_obs, t_in_win:dt_betweenum_obs:t_out_win];
        jWin = jWin + 1;
    end
    
    num_obs     = size(et_obs,2);
    x_truth_stn = zeros(3, num_obs);
    v_truth_stn = zeros(3, num_obs);
    
    j1 = 1;
    while j1 <= num_obs
        dt_min = (et_obs(j1) - tle_epoch_et)/60;
        [~, r_teme_k, v_teme_k] = sgp4(sat_rec, dt_min);
        years_k = cspice_unitim(et_obs(j1), 'ET', 'TDT') / cspice_jyear();
        time_k   = years_k / 100;
        [r_eci_k, v_eci_k, ~] = teme2eci(r_teme_k(:), v_teme_k(:), accel_teme, time_k, nutation_dpsi, nutation_deps);
        x_truth_stn(:,j1) = r_eci_k;
        v_truth_stn(:,j1) = v_eci_k;
        j1 = j1 + 1;
    end
    
    state_truth_stn = [x_truth_stn; v_truth_stn];
    [az_ref_stn, el_ref_stn, rho_ref_stn, ~] = antenna_pointing_rad(stn_name, et_obs, state_truth_stn);
    
    % Add noise
    el_noisy  = NaN(1, num_obs);
    az_noisy  = NaN(1, num_obs);
    rho_noisy = NaN(1, num_obs);
    
    j = 1;
    while j <= num_obs
        az_k  = az_ref_stn(j)  + noise_az_rad  * randn();
        el_k  = el_ref_stn(j)  + noise_el_rad  * randn();
        rho_k = rho_ref_stn(j) + noise_rho_km  * randn();
        if el_k > elev_mask_rad
            az_noisy(j)  = mod(az_k, 2*pi);
            el_noisy(j)  = el_k;
            rho_noisy(j) = rho_k;
        end
        j = j + 1;
    end
    
    % Step #2.4: Calculate residuals and analyze measurement performance
    passed_idx   = ~isnan(az_noisy); count_samples = sum(passed_idx);
    dRho = NaN(1, num_obs);
    dAz  = NaN(1, num_obs);
    dEl  = NaN(1, num_obs);
    
    j = 1;
    while j <= num_obs
        if passed_idx(j)
            wrap_az = az_noisy(j) - az_ref_stn(j);
            dAz(j)  = angle(exp(1i * wrap_az));
            dEl(j)  = el_noisy(j)  - el_ref_stn(j);
            dRho(j) = rho_noisy(j) - rho_ref_stn(j);
        end
        j = j + 1;
    end

    dAz_deg   = dAz  * 180/pi;dEl_deg   = dEl  * 180/pi;
    bound_az  = 3*noise_az_deg;bound_el  = 3*noise_el_deg;
    bound_rho = 3*noise_rho_km; t_obs_hr  = (et_obs - epoch_ref)/3600;

    fprintf('====================================================================================\n');
    Tbl_cnt = table({stn_name}, num_obs, count_samples, ...
        'VariableNames', {'Station', 'Candidate_Epochs', 'Valid_Measurements'});
    disp(Tbl_cnt);

       if count_samples > 0
            sigma_dRho = std(dRho(passed_idx));
            sigma_dEl  = std(dEl_deg(passed_idx));
            sigma_dAz  = std(dAz_deg(passed_idx));
            mean_dRho = mean(dRho(passed_idx));         
            mean_dEl  = mean(dEl_deg(passed_idx)); 
            mean_dAz  = mean(dAz_deg(passed_idx));
    
            norm_dRho = dRho(passed_idx)/noise_rho_km;
            norm_dEl  = dEl_deg(passed_idx)/noise_el_deg;
            norm_dAz  = dAz_deg(passed_idx)/noise_az_deg;
    
            Means_v   = [mean_dRho; mean_dEl; mean_dAz];
            Stds_v    = [sigma_dRho; sigma_dEl; sigma_dAz];
            MaxAbs_v  = [max(abs(dRho(passed_idx))); max(abs(dEl_deg(passed_idx))); max(abs(dAz_deg(passed_idx)))];
            Thresh_v  = [bound_rho; bound_el; bound_az];
            Outlier_v = [100*mean(abs(norm_dRho)>3); 100*mean(abs(norm_dEl)>3); 100*mean(abs(norm_dAz)>3)];
    
            Tbl_res = table(Means_v, Stds_v, MaxAbs_v, Thresh_v, Outlier_v, ...
                'VariableNames', {'μ', 'σ', 'max_|Res|', '3σ_Bound', 'Outliers [%]'}, ...
                'RowNames', {'Range [km]', 'Elevation [deg]', 'Azimuth [deg]'});
            disp(Tbl_res);
        end
    
    % Step #2.5: Plot residuals and normalized residual histograms
    clr_res   = [0.4940 0.1840 0.5560];
    clr_hist_rho = [0.3010 0.7450 0.9330];
    clr_hist_az  = [0.6350 0.0780 0.1840];
    clr_hist_el  = [0.4660 0.6740 0.1880];

    figure('Name',sprintf('Residuals_%s',stn_name));
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    sgtitle(['Measurement Residuals Analysis: ',stn_name]);

    subplot(3,1,1)
    plot(t_obs_hr, dRho, '.', 'Color', clr_res, 'DisplayName','$\Delta Range$'); hold on
    yline(bound_rho,'k--','DisplayName','$\pm 3\sigma$'); 
    yline(-bound_rho,'k--','HandleVisibility','off');
    yline(0,'r-','HandleVisibility','off'); 
    grid on;
    ylabel('$\Delta Range$ [km]'); 
    title('Range Residual Evolution'); 
    legend('Location','best');
    
    subplot(3,1,2)
    plot(t_obs_hr, dAz_deg, '.', 'Color', clr_res, 'DisplayName','$\Delta Az$'); hold on
    yline(bound_az,'k--','DisplayName','$\pm 3\sigma$'); 
    yline(-bound_az,'k--','HandleVisibility','off');
    yline(0,'r-','HandleVisibility','off'); 
    grid on;
    ylabel('$\Delta Az~[^\circ]$'); 
    title('Azimuth Residual Evolution'); 
    legend('Location','best');
    
    subplot(3,1,3)
    plot(t_obs_hr, dEl_deg, '.', 'Color', clr_res, 'DisplayName','$\Delta El$'); hold on
    yline(bound_el,'k--','DisplayName','$\pm 3\sigma$'); 
    yline(-bound_el,'k--','HandleVisibility','off');
    yline(0,'r-','HandleVisibility','off'); 
    grid on;
    xlabel('$t-t_i$ [h]'); 
    ylabel('$\Delta El~[^\circ]$'); 
    title('Elevation Residual Evolution'); 
    legend('Location','best');

   if count_samples > 0
           hist_edges = -6:0.29:6;
    figure('Name',sprintf('HistNorm_%s',stn_name));
    sgtitle(['Normalized Residual Distribution at ', stn_name], 'Interpreter', 'latex');
    
    subplot(3,1,1)
    histogram(norm_dRho,hist_edges,'FaceColor',clr_hist_rho,'FaceAlpha',0.7,'EdgeAlpha',0.2,'DisplayName','$(y_{\rho,\mathrm{meas}} - y_{\rho,\mathrm{true}})/\sigma_\rho$'); hold on
    xline(3,'k--','DisplayName','$\pm 3\sigma$'); xline(-3,'k--','HandleVisibility','off'); grid on;
    ylabel('Count'); title('Range Normalized Residuals'); legend('Location','best');
    
    subplot(3,1,2)
    histogram(norm_dEl,hist_edges,'FaceColor',clr_hist_el,'FaceAlpha',0.7,'EdgeAlpha',0.2,'DisplayName','$(y_{El,\mathrm{meas}} - y_{El,\mathrm{true}})/\sigma_{El}$'); hold on
    xline(3,'k--','DisplayName','$\pm 3\sigma$'); xline(-3,'k--','HandleVisibility','off'); grid on;
    ylabel('Count'); title('Elevation Normalized Residuals'); legend('Location','best');
    
    subplot(3,1,3)
    histogram(norm_dAz,hist_edges,'FaceColor',clr_hist_az,'FaceAlpha',0.7,'EdgeAlpha',0.2,'DisplayName','$(y_{Az,\mathrm{meas}} - y_{Az,\mathrm{true}})/\sigma_{Az}$'); hold on
    xline(3,'k--','DisplayName','$\pm 3\sigma$'); xline(-3,'k--','HandleVisibility','off'); grid on;
    ylabel('Count');
    title('Azimuth Normalized Residuals'); legend('Location','best');
    end

    % Step #2.6: Final temporal verification and data storage
    t_hr_obs_valid = et_obs(passed_idx)/cspice_spd();
    figure('Name',['Temporal Verification: ',stn_name]);
    subplot(3,1,1); hold on; grid on;
    plot(t_obs_hr, az_ref_stn*180/pi,'k','LineWidth',1.5,'DisplayName','True');
    plot(t_obs_hr(passed_idx), az_noisy(passed_idx)*180/pi,'b.','MarkerSize',8,'DisplayName','Noisy');
    plot(t_obs_hr, (az_ref_stn + noise_az_rad)*180/pi,'r--');
    plot(t_obs_hr, (az_ref_stn - noise_az_rad)*180/pi,'r--','DisplayName','$\pm\sigma$');
    title(['Azimuth @ ',stn_name]); ylabel('[deg]'); legend('Location','best');
    subplot(3,1,2); hold on; grid on;
    plot(t_obs_hr, el_ref_stn*180/pi,'k','LineWidth',1.5,'DisplayName','True');
    plot(t_obs_hr(passed_idx), el_noisy(passed_idx)*180/pi,'b.','MarkerSize',8,'DisplayName','Noisy');
    plot(t_obs_hr, (el_ref_stn + noise_el_rad)*180/pi,'r--');
    plot(t_obs_hr, (el_ref_stn - noise_el_rad)*180/pi,'r--','DisplayName','$\pm\sigma$');
    title(['Elevation @ ',stn_name]); ylabel('[deg]'); legend('Location','best');
    subplot(3,1,3); hold on; grid on;
    plot(t_obs_hr, rho_ref_stn,'k','LineWidth',1.5,'DisplayName','True');
    plot(t_obs_hr(passed_idx), rho_noisy(passed_idx),'b.','MarkerSize',8,'DisplayName','Noisy');
    plot(t_obs_hr, rho_ref_stn + noise_rho_km,'r--');
    plot(t_obs_hr, rho_ref_stn - noise_rho_km,'r--','DisplayName','$\pm\sigma$');
    title(['Range @ ',stn_name]); ylabel('[km]'); xlabel('$t-t_i$ [h]'); legend('Location','best');

    meas_store.(stn_name).et_obs  = et_obs;
    meas_store.(stn_name).az_noisy   = az_noisy;
    meas_store.(stn_name).el_noisy   = el_noisy;
    meas_store.(stn_name).rho_noisy  = rho_noisy;
    meas_store.(stn_name).passed_idx   = passed_idx;
end
fprintf('--- Point 2 Complete ---\n');
%% Step#3+4: ORBIT DETERMINATION (6 cases)
fprintf('\n--- POINT 3+4: Orbit Determination ---\n');
format long

% Step #3.1: Define navigation epoch and optimization settings
nav_epoch_et = epoch_ref;
lsq_opts     = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter');

% Step #3.2: Compute weighting matrices for residuals (Angles-only and Full)
R_ang_only = diag([noise_az_rad^2, noise_el_rad^2]);
W_ang_only = inv(sqrtm(R_ang_only));
R_full     = diag([noise_az_rad^2, noise_el_rad^2, noise_rho_km^2]);
W_full     = inv(sqrtm(R_full));

% Step #3.3: Extract and format valid measurements for TROLL and KOUROU
% TROLL valid measurement arrays
valid_T       = meas_store.TROLL.passed_idx;
epochs_T      = meas_store.TROLL.et_obs(valid_T);
az_T          = meas_store.TROLL.az_noisy(valid_T);
el_T          = meas_store.TROLL.el_noisy(valid_T);
rho_T         = meas_store.TROLL.rho_noisy(valid_T);
obs_ang_T     = [az_T(:), el_T(:)];
obs_full_T    = [az_T(:), el_T(:), rho_T(:)];

% KOUROU valid measurement arrays
valid_K       = meas_store.KOUROU.passed_idx;
epochs_K      = meas_store.KOUROU.et_obs(valid_K);
az_K          = meas_store.KOUROU.az_noisy(valid_K);
el_K          = meas_store.KOUROU.el_noisy(valid_K);
rho_K         = meas_store.KOUROU.rho_noisy(valid_K);
obs_full_K    = [az_K(:), el_K(:), rho_K(:)];

% Step #3.4: Perform Orbit Determination for all 6 cases (i through vi)
    % Case (i): TROLL angles-only, Kepler
    cost_i   = @(x) get_station_residuals(x, nav_epoch_et, epochs_T, W_ang_only, obs_ang_T, 'TROLL', 'Kepler', grav_param);
    res0_i   = cost_i(guess_kep);
    [x_od_i, rn_i, r_i, ef_i, ~,~, J_i] = lsqnonlin(cost_i, guess_kep, [],[],lsq_opts);
    J_i = full(J_i);
    n_st = numel(guess_kep); m_i = numel(r_i); 
    s2_i = rn_i/(m_i-n_st);
    P_od_i   = s2_i * inv(J_i.'*J_i);
    sig_r_i  = sqrt(trace(P_od_i(1:3,1:3)));
    sig_v_i  = sqrt(trace(P_od_i(4:6,4:6)));
    sig_a_i  = a_covariance_to_std(x_od_i, P_od_i, grav_param);
    
    lbl_i = 'Case (i): TROLL, angles-only, Kepler';
    met_i = ["Position Error (||dr||) [km]"; "Velocity Error (||dv||) [km/s]"; ...
             "Uncertainty (sigma_r) [km]"; "Uncertainty (sigma_v) [km/s]"; "Uncertainty (sigma_a) [km]"];
    val_i = [norm(x_od_i(1:3) - x_truth_ref(1:3)); norm(x_od_i(4:6) - x_truth_ref(4:6)); ...
             sig_r_i; sig_v_i; sig_a_i];
    fprintf('\n=== %s ===\n', lbl_i);
    disp(table(met_i, val_i, 'VariableNames', {'Metric', 'Value'}));
    plot_od_residuals('Case (i): TROLL ang-only (Kepler)', res0_i, r_i, W_ang_only, R_ang_only, ...
         epochs_T, [], noise_az_rad, noise_el_rad, noise_rho_km);
    % Case (ii): TROLL angles-only, J2
    cost_ii  = @(x) get_station_residuals(x, nav_epoch_et, epochs_T, W_ang_only, obs_ang_T, 'TROLL', 'J2', grav_param);
    res0_ii  = cost_ii(guess_j2);
    [x_od_ii, rn_ii, r_ii, ef_ii, ~,~, J_ii] = lsqnonlin(cost_ii, guess_j2, [],[],lsq_opts);
    J_ii = full(J_ii);
    m_ii = numel(r_ii);
    s2_ii = rn_ii/(m_ii-n_st);
    P_od_ii  = s2_ii * inv(J_ii.'*J_ii);
    sig_r_ii = sqrt(trace(P_od_ii(1:3,1:3)));
    sig_v_ii = sqrt(trace(P_od_ii(4:6,4:6)));
    sig_a_ii = a_covariance_to_std(x_od_ii, P_od_ii, grav_param);
    
    lbl_ii = 'Case (ii): TROLL, angles-only, J2';
    val_ii = [norm(x_od_ii(1:3) - x_truth_ref(1:3)); norm(x_od_ii(4:6) - x_truth_ref(4:6)); ...
              sig_r_ii; sig_v_ii; sig_a_ii];
    fprintf('\n=== %s ===\n', lbl_ii);
    disp(table(met_i, round(val_ii,4), 'VariableNames', {'Metric', 'Value'}));
    plot_od_residuals('Case (ii): TROLL ang-only (J2)', res0_ii, r_ii, W_ang_only, R_ang_only, ...
        epochs_T, [], noise_az_rad, noise_el_rad, noise_rho_km);


    % Case (iii): TROLL angles+range, Kepler
    cost_iii = @(x) get_station_residuals(x, nav_epoch_et, epochs_T, W_full, obs_full_T, 'TROLL', 'Kepler', grav_param);
    res0_iii = cost_iii(guess_kep);
    [x_od_iii, rn_iii, r_iii, ef_iii, ~,~, J_iii] = lsqnonlin(cost_iii, guess_kep, [],[],lsq_opts);
    J_iii = full(J_iii);
    m_iii = numel(r_iii);
    s2_iii = rn_iii/(m_iii-n_st);
    P_od_iii  = s2_iii * inv(J_iii.'*J_iii);
    sig_r_iii = sqrt(trace(P_od_iii(1:3,1:3)));
    sig_v_iii = sqrt(trace(P_od_iii(4:6,4:6)));
    sig_a_iii = a_covariance_to_std(x_od_iii, P_od_iii, grav_param);
    
    lbl_iii = 'Case (iii): TROLL, angles + range, Kepler';
    val_iii = [norm(x_od_iii(1:3) - x_truth_ref(1:3)); norm(x_od_iii(4:6) - x_truth_ref(4:6)); ...
               sig_r_iii; sig_v_iii; sig_a_iii];
    fprintf('\n=== %s ===\n', lbl_iii);
    disp(table(met_i, round(val_iii,4), 'VariableNames', {'Metric', 'Value'}));
    plot_od_residuals('Case (iii): TROLL ang+range (Kepler)', res0_iii, r_iii, W_full, R_full, ...
        epochs_T, [], noise_az_rad, noise_el_rad, noise_rho_km);

    % Case (iv): TROLL angles+range, J2
    cost_iv  = @(x) get_station_residuals(x, nav_epoch_et, epochs_T, W_full, obs_full_T, 'TROLL', 'J2', grav_param);
    res0_iv  = cost_iv(guess_j2);
    [x_od_iv, rn_iv, r_iv, ef_iv, ~,~, J_iv] = lsqnonlin(cost_iv, guess_j2, [],[],lsq_opts);
    J_iv = full(J_iv);
    m_iv = numel(r_iv);
    s2_iv = rn_iv/(m_iv-n_st);
    P_od_iv  = s2_iv * inv(J_iv.'*J_iv);
    sig_r_iv = sqrt(trace(P_od_iv(1:3,1:3)));
    sig_v_iv = sqrt(trace(P_od_iv(4:6,4:6)));
    sig_a_iv = a_covariance_to_std(x_od_iv, P_od_iv, grav_param);
    
    lbl_iv = 'Case (iv): TROLL, angles + range, J2';
    val_iv = [norm(x_od_iv(1:3) - x_truth_ref(1:3)); norm(x_od_iv(4:6) - x_truth_ref(4:6)); ...
              sig_r_iv; sig_v_iv; sig_a_iv];
    fprintf('\n=== %s ===\n', lbl_iv);
    disp(table(met_i, round(val_iv,4), 'VariableNames', {'Metric', 'Value'}));
    plot_od_residuals('Case (iv): TROLL ang+range (J2)', res0_iv, r_iv, W_full, R_full, ...
        epochs_T, [], noise_az_rad, noise_el_rad, noise_rho_km);

    cost_v   = @(x) get_multi_station_residuals(x, nav_epoch_et, epochs_T, obs_full_T, epochs_K, obs_full_K, W_full, 'Kepler', grav_param);
    res0_v   = cost_v(guess_kep);
    [x_od_v, rn_v, r_v, ef_v, ~,~, J_v] = lsqnonlin(cost_v, guess_kep, [],[],lsq_opts);
    J_v = full(J_v);
    m_v = numel(r_v);
    s2_v = rn_v/(m_v-n_st);
    P_od_v   = s2_v * inv(J_v.'*J_v);
    sig_r_v  = sqrt(trace(P_od_v(1:3,1:3)));
    sig_v_v  = sqrt(trace(P_od_v(4:6,4:6)));
    sig_a_v  = a_covariance_to_std(x_od_v, P_od_v, grav_param);
    
    lbl_v = 'Case (v): TROLL + KOUROU, angles + range, Kepler';
    val_v = [norm(x_od_v(1:3) - x_truth_ref(1:3)); norm(x_od_v(4:6) - x_truth_ref(4:6)); ...
             sig_r_v; sig_v_v; sig_a_v];
    fprintf('\n=== %s ===\n', lbl_v);
    disp(table(met_i, round(val_v,4), 'VariableNames', {'Metric', 'Value'}));
    plot_od_residuals('Case (v): TR+KO ang+range (Kepler)', res0_v, r_v, W_full, R_full, ...
        epochs_T, epochs_K, noise_az_rad, noise_el_rad, noise_rho_km);


    % Case (vi): TROLL+KOUROU angles+range, J2
    cost_vi  = @(x) get_multi_station_residuals(x, nav_epoch_et, epochs_T, obs_full_T, epochs_K, obs_full_K, W_full, 'J2', grav_param);
    res0_vi  = cost_vi(guess_j2);
    [x_od_vi, rn_vi, r_vi, ef_vi, ~,~, J_vi] = lsqnonlin(cost_vi, guess_j2, [],[],lsq_opts);
    J_vi = full(J_vi);
    m_vi = numel(r_vi);
    s2_vi = rn_vi/(m_vi-n_st);
    P_od_vi  = s2_vi * inv(J_vi.'*J_vi);
    sig_r_vi = sqrt(trace(P_od_vi(1:3,1:3)));
    sig_v_vi = sqrt(trace(P_od_vi(4:6,4:6)));
    sig_a_vi = a_covariance_to_std(x_od_vi, P_od_vi, grav_param);
    
    lbl_vi = 'Case (vi): TROLL + KOUROU, angles + range, J2';
    val_vi = [norm(x_od_vi(1:3) - x_truth_ref(1:3)); norm(x_od_vi(4:6) - x_truth_ref(4:6)); ...
              sig_r_vi; sig_v_vi; sig_a_vi];
    fprintf('\n=== %s ===\n', lbl_vi);
    disp(table(met_i, round(val_vi,4), 'VariableNames', {'Metric', 'Value'}));
    plot_od_residuals('Case (vi): TR+KO ang+range (J2)', res0_vi, r_vi, W_full, R_full, ...
        epochs_T, epochs_K, noise_az_rad, noise_el_rad, noise_rho_km);

% Step #3.5: Summarize navigation solutions
X_nav = [x_od_i, x_od_ii, x_od_iii, x_od_iv, x_od_v, x_od_vi];
nav_case_names = ["Case (i)"; "Case (ii)"; "Case (iii)"; "Case (iv)"; "Case (v)"; "Case (vi)"];

Tbl_nav = table(nav_case_names, X_nav(1,:)', X_nav(2,:)', X_nav(3,:)', ...
              X_nav(4,:)', X_nav(5,:)', X_nav(6,:)', ...
              'VariableNames', {'Case', 'rx[km]', 'ry[km]', 'rz[km]', 'vx[km/s]', 'vy[km/s]', 'vz[km/s]'});
fprintf('\n=== NAVIGATION SOLUTIONS AT t_i ===\n');
disp(Tbl_nav);

%%  Step#5: PERFORMANCE ASSESSMENT
fprintf('\n===================== POINT 5: PERFORMANCE vs SGP4 =====================\n');

% Step #5.1: Define propagation time span and calculate orbital period
dt_step   = 60;
N_periods = 6;

r_ref = norm(x_truth_ref(1:3));
v_ref = norm(x_truth_ref(4:6));
a_ref = 1/(2/r_ref - v_ref^2/grav_param);     % vis-viva: a = 1/(2/r - v²/μ)
T_orbit = 2*pi*sqrt(a_ref^3/grav_param);        % Kepler's third law

t_prop = epoch_ref : dt_step : (epoch_ref + N_periods*T_orbit);
n_steps = numel(t_prop);
t_hr_prop  = (t_prop - epoch_ref)/3600;

% Step #5.2: Propagate SGP4 truth trajectory over the defined time span
pos_sgp4 = zeros(3, n_steps); vel_sgp4 = zeros(3, n_steps);
j = 1;
while j <= n_steps
    dt_min = (t_prop(j) - tle_epoch_et)/60;
    [~, r_teme_p, v_teme_p] = sgp4(sat_rec, dt_min);
    tdt_years = cspice_unitim(t_prop(j), 'ET', 'TDT') / cspice_jyear();
    time_p     = tdt_years / 100;
    [r_eci_p, v_eci_p, ~] = teme2eci(r_teme_p(:), v_teme_p(:), accel_teme, time_p, nutation_dpsi, nutation_deps);
    pos_sgp4(:,j) = r_eci_p;
    vel_sgp4(:,j) = v_eci_p;
    j = j + 1;
end
state_sgp4_ref = [pos_sgp4; vel_sgp4];

% Step #5.3: Propagate determined states using J2 dynamics for comparison
dynamics_j2 = @(t,x) keplerianJ2_rhs(t, x, grav_param, earth_radius, J2_coeff);
prop_opts  = odeset('RelTol',1e-9,'AbsTol',[1e-7*ones(3,1);1e-10*ones(3,1)]);

[~, prop_i]   = ode113(dynamics_j2, t_prop(:), x_od_i,   prop_opts);
[~, prop_ii]  = ode113(dynamics_j2, t_prop(:), x_od_ii,  prop_opts);
[~, prop_iii] = ode113(dynamics_j2, t_prop(:), x_od_iii, prop_opts);
[~, prop_iv]  = ode113(dynamics_j2, t_prop(:), x_od_iv,  prop_opts);
[~, prop_v]   = ode113(dynamics_j2, t_prop(:), x_od_v,   prop_opts);
[~, prop_vi]  = ode113(dynamics_j2, t_prop(:), x_od_vi,  prop_opts);

prop_i   = prop_i.';   prop_ii  = prop_ii.';  prop_iii = prop_iii.';
prop_iv  = prop_iv.';  prop_v   = prop_v.';   prop_vi  = prop_vi.';

% Step #5.4: Calculate position/velocity error norms relative to SGP4
err_r_i   = vecnorm(prop_i(1:3,:)   - state_sgp4_ref(1:3,:));
err_v_i   = vecnorm(prop_i(4:6,:)   - state_sgp4_ref(4:6,:));
err_r_ii  = vecnorm(prop_ii(1:3,:)  - state_sgp4_ref(1:3,:));
err_v_ii  = vecnorm(prop_ii(4:6,:)  - state_sgp4_ref(4:6,:));
err_r_iii = vecnorm(prop_iii(1:3,:) - state_sgp4_ref(1:3,:));
err_v_iii = vecnorm(prop_iii(4:6,:) - state_sgp4_ref(4:6,:));
err_r_iv  = vecnorm(prop_iv(1:3,:)  - state_sgp4_ref(1:3,:));
err_v_iv  = vecnorm(prop_iv(4:6,:)  - state_sgp4_ref(4:6,:));
err_r_v   = vecnorm(prop_v(1:3,:)   - state_sgp4_ref(1:3,:));
err_v_v   = vecnorm(prop_v(4:6,:)   - state_sgp4_ref(4:6,:));
err_r_vi  = vecnorm(prop_vi(1:3,:)  - state_sgp4_ref(1:3,:));
err_v_vi  = vecnorm(prop_vi(4:6,:)  - state_sgp4_ref(4:6,:));

% Step #5.5: Generate performance metrics table and comparative plots

perf_labels = ["(i) TROLL Ang (Kep)"; "(ii) TROLL Ang (J2)"; "(iii) TROLL All (Kep)"; ...
               "(iv) TROLL All (J2)"; "(v) TR+KO All (Kep)"; "(vi) TR+KO All (J2)"];
cell_r = {err_r_i, err_r_ii, err_r_iii, err_r_iv, err_r_v, err_r_vi};
cell_v = {err_v_i, err_v_ii, err_v_iii, err_v_iv, err_v_v, err_v_vi};

max_dr  = cellfun(@max, cell_r)';
rms_dr  = cellfun(@(x) sqrt(mean(x.^2)), cell_r)';
fin_dr  = cellfun(@(x) x(end), cell_r)';
max_dv  = cellfun(@max, cell_v)';
rms_dv  = cellfun(@(x) sqrt(mean(x.^2)), cell_v)';
fin_dv  = cellfun(@(x) x(end), cell_v)';

Tbl_perf = table(perf_labels, max_dr, rms_dr, fin_dr, max_dv, rms_dv, fin_dv, ...
    'VariableNames', {'Case', 'Max_dr[km]', 'RMS_dr[km]', 'Final_dr[km]', ...
                      'Max_dv[km/s]', 'RMS_dv[km/s]', 'Final_dv[km/s]'});
fprintf('\n=== PERFORMANCE SUMMARY vs SGP4 TRUTH ===\n');
disp(Tbl_perf);

plot_palette = [0.93, 0.11, 0.14; 0.12, 0.47, 0.71; 0.44, 0.62, 0.29; ...
                0.26, 0.26, 0.26; 0.94, 0.60, 0.07; 0.58, 0.44, 0.86];
t_mjd2k = (t_prop - cspice_str2et('2000-01-01T12:00:00'))/86400;

figure('Name','OD Log Performance');
case_legend_str = {'i  TROLL ang-only Kepler', 'ii TROLL ang-only J2', ...
'iii TROLL ang+range Kepler', 'iv  TROLL ang+range J2', ...
'v  TROLL+KOUROU ang+range Kepler', 'vi TROLL+KOUROU ang+range J2'};
subplot(2,1,1); hold on; grid on;
for ii = 1:6
    plot(t_hr_prop, cell_r{ii}, 'Color', plot_palette(ii,:), 'LineWidth', 1.5);
end
set(gca, 'YScale', 'log');
ylabel('$\log_{10}(\|\Delta r\|)$ [km]');
%title('Error Norms: $\varepsilon_{pos}$');
title('Position Error Norm: $\|\Delta\mathbf{r}\|$');
legend(case_legend_str, 'Location', 'best');
subplot(2,1,2); hold on; grid on;
for ii = 1:6
    plot(t_hr_prop, cell_v{ii}, 'Color', plot_palette(ii,:), 'LineWidth', 1.5);
end
set(gca, 'YScale', 'log');
xlabel('$t - t_i$ [h]');
ylabel('$\log_{10}(\|\Delta v\|)$ [km/s]');
%title('Error Norms: $\varepsilon_{vel}$');
title('Velocity Error Norm: $\|\Delta\mathbf{v}\|$');

% Step #5.6: 3D Visualization of orbit comparison
traj3d_i   = prop_i(1:3,:).';
traj3d_ii  = prop_ii(1:3,:).';
traj3d_iii = prop_iii(1:3,:).';
traj3d_iv  = prop_iv(1:3,:).';
traj3d_v   = prop_v(1:3,:).';
traj3d_vi  = prop_vi(1:3,:).';
pos_sgp4_3d = state_sgp4_ref(1:3,:);

figure('Name','3D Orbit Comparison');
hold on; grid on; axis equal;
R_E_plot = cspice_bodvrd('EARTH','RADII',3); R_E_plot = R_E_plot(1);
[Xs, Ys, Zs] = sphere(60);
if exist('Earth_IMAGE.jpg','file')
    earth_img = imread('Earth_IMAGE.jpg');
    surf(Xs*R_E_plot, Ys*R_E_plot, Zs*R_E_plot, ...
        'FaceColor','texturemap','CData',earth_img,'EdgeColor','none');
    lighting gouraud; camlight headlight; material dull;
else
    surf(Xs*R_E_plot, Ys*R_E_plot, Zs*R_E_plot, ...
        'FaceColor',[0.5 0.7 0.9],'FaceAlpha',0.5,'EdgeColor','none');
end
h_sgp4 = plot3(pos_sgp4_3d(1,:), pos_sgp4_3d(2,:), pos_sgp4_3d(3,:), ...
    'k-',  'LineWidth', 2.5, 'DisplayName','SGP4 Truth');
h_ci   = plot3(traj3d_i(:,1),  traj3d_i(:,2),  traj3d_i(:,3),  '--', ...
    'Color',[1.0 0.6 0.7], 'LineWidth',1.2, 'DisplayName','(i) Angles (Kep)');
h_cii  = plot3(traj3d_ii(:,1), traj3d_ii(:,2), traj3d_ii(:,3), '-',  ...
    'Color',[0.8 0.0 0.0], 'LineWidth',1.0, 'DisplayName','(ii) Angles (J2)');
h_civ  = plot3(traj3d_iv(:,1), traj3d_iv(:,2), traj3d_iv(:,3), '-',  ...
    'Color',[0.0 0.0 0.8], 'LineWidth',1.2, 'DisplayName','(iv) TROLL All (J2)');
h_cvi  = plot3(traj3d_vi(:,1), traj3d_vi(:,2), traj3d_vi(:,3), '-',  ...
    'Color',[0.0 0.6 0.0], 'LineWidth',2.0, 'DisplayName','(vi) TR+KO (J2) [Best]');

legend([h_sgp4, h_ci, h_cii, h_civ, h_cvi], 'Location','bestoutside');
title('Orbit Propagation Comparison @EARTH J2000');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
view(3);
fprintf('--- Point 5 Complete ---\n');
 %cspice_kclear();

%% =====================================================================
%  Helper FUNCTIONS:
%% =====================================================================
function dxdt = keplerianJ2_rhs(t, x, mu, Re, J2)
% AIM: Computes the 6-state time derivative including J2 Earth perturbation.
% INPUTS: 
%   t  : [1x1] Current time (ET) [s].
%   x  : [6x1] State vector [pos_x, pos_y, pos_z, vel_x, vel_y, vel_z] [km, km/s].
%   mu : [1x1] Earth gravitational parameter [km^3/s^2].
%   Re : [1x1] Earth equatorial radius [km].
%   J2 : [1x1] J2 perturbation coefficient [dimensionless].
% OUTPUT: 
%   dxdt : [6x1] State derivative [vel_x, ..., acc_z] [km/s, km/s^2].

    r_ECI=x(1:3);
    a_kep = -mu*r_ECI/(norm(r_ECI))^3;
    R_eci2ecef   = cspice_pxform('J2000','IAU_EARTH',t);
    r_ECEF =R_eci2ecef*r_ECI;
    a_ECEF = 1.5*J2*mu*Re^2/((sqrt(dot(r_ECEF,r_ECEF)))^5)*[r_ECEF(1)*(5*((r_ECEF(3))^2)/dot(r_ECEF,r_ECEF)-1);
                                   r_ECEF(2)*(5*((r_ECEF(3))^2)/dot(r_ECEF,r_ECEF)-1);
                                   r_ECEF(3)*(5*((r_ECEF(3))^2)/dot(r_ECEF,r_ECEF)-3)];
    dxdt = zeros(6,1);
    dxdt(1:3) = x(4:6);
    dxdt(4:6) = a_kep +R_eci2ecef'*a_ECEF;
end

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



function [xf, tt, xx] = keplerian_propagator(x0, tspan, attractor, effect, opts)
% AIM: Integrates orbit from initial condition x0 over tspan.
% INPUTS: 
%   x0        : [6x1] Initial state [km, km/s].
%   tspan     : [1xN] Vector of time points for integration [s].
%   attractor : [char] Name of central body (e.g., 'EARTH').
%   effect    : [char] Dynamics model: 'NONE' (Kep) or 'J2'.
% OUTPUTS: 
%   xf        : [6x1] Final propagated state [km, km/s].
%   tt        : [Nx1] Time vector [s].
%   xx        : [Nx6] State history [km, km/s].
    mu_p = cspice_bodvrd(attractor,'GM',1); mu_p=mu_p(1);
    Re_p = cspice_bodvrd(attractor,'RADII',3); Re_p=Re_p(1);

    % FIX: match the tighter tolerances used in the OD cost function
    if nargin < 5 || isempty(opts)
        opts = odeset('RelTol',1e-10,'AbsTol',[1e-6*ones(3,1);1e-9*ones(3,1)]);
    end
    if strcmp(effect,'J2')
        J2l=0.0010826269;
        [tt,xx] = ode113(@(t,x) keplerianJ2_rhs(t,x,mu_p,Re_p,J2l), tspan, x0(:), opts);
    else
        [tt,xx] = ode113(@(t,x) keplerian_rhs(t,x,mu_p), tspan, x0(:), opts);
    end
    xf = xx(end,:)';
end

function [azimuth, elevation, range, range_rate] = antenna_pointing_rad(station_name, et_vec, x_satellite_eci)
% AIM: Computes topocentric observables for a spacecraft relative to a ground station.
% INPUTS:
%   station_name : [char]   Name of the ground station (defined in kernels).
%   et_vec       : [1xN]    Vector of ephemeris times [s].
%   rv_sat_eci   : [6xN]    Spacecraft ECI state (position [km], velocity [km/s]).
% OUTPUTS:
%   azimuth      : [1xN]    Azimuth angle [rad].
%   elevation    : [1xN]    Elevation angle [rad].
%   range        : [1xN]    Slant range [km].
%   range_rate   : [1xN]    Range rate [km/s].

    et_vec      = et_vec(:).';
    n_points    = numel(et_vec);

    range       = zeros(1, n_points);
    range_rate  = zeros(1, n_points);
    azimuth     = zeros(1, n_points);
    elevation   = zeros(1, n_points);


    topo_frame  = [station_name, '_TOPO'];

    for st_idx = 1:n_points
        et = et_vec(st_idx);

        [x_station_eci, ~] = cspice_spkezr(station_name, et, 'J2000', 'NONE', 'EARTH');
        x_rel_eci = x_satellite_eci(:, st_idx) - x_station_eci;

        R_eci2topo  = cspice_sxform('J2000', topo_frame, et);
        rv_rel_topo = R_eci2topo * x_rel_eci;

        required_state   = cspice_xfmsta(rv_rel_topo, 'RECTANGULAR', 'LATITUDINAL', 'EARTH');

        range(st_idx)      = required_state(1);
        azimuth(st_idx)    = mod(required_state(2), 2*pi);
        elevation(st_idx)  = required_state(3);
        range_rate(st_idx) = required_state(4);
    end
end



function res_vec = get_station_residuals(x_estimate, et_start, et_meas_vec, ...
        W_matrix, meas_data, station_name, dynamic_model, mu)
% AIM: Builds the weighted residual vector for a single ground station OD problem.
% INPUTS:
%   x_estimate    : [6x1]   Current state estimate [km; km/s].
%   et_start      : [1x1]   Navigation epoch (integration origin) [s ET].
%   et_meas_vec   : [Nx1]   Observation epochs [s ET].
%   W_matrix      : [PxP]   Measurement weight matrix (inv of noise sqrt).
%   meas_data     : [NxP]   Observed Az [rad], El [rad], (Range [km]).
%   station_name  : [char]  Ground station identifier string.
%   dynamic_model : [char]  Force model selector: 'Kepler' or 'J2'.
%   mu            : [1x1]   Earth gravitational parameter [km^3/s^2].
% OUTPUT:
%   res_vec       : [N*Px1] Stacked weighted residuals.

    n_epochs   = numel(et_meas_vec);
    n_meas_typ = size(meas_data, 2);

    earth_radii = cspice_bodvrd('EARTH', 'RADII', 3);
    Re_eq       = earth_radii(1);
    J2_val      = 0.0010826269;

    ode_settings = odeset('RelTol', 1e-12, 'AbsTol', [1e-6*ones(3,1); 1e-9*ones(3,1)]);

    if strcmp(dynamic_model, 'J2')
        eom_handle = @(t, state) keplerianJ2_rhs(t, state, mu, Re_eq, J2_val);
    else
        eom_handle = @(t, state) keplerian_rhs(t, state, mu);
    end

    epoch_span    = [et_start; et_meas_vec(:)];
    [~, traj_all] = ode113(eom_handle, epoch_span, x_estimate(:), ode_settings);
    traj_eci      = traj_all(2:end, :).';

    [az_computed, el_computed, rho_computed, ~] = ...
        antenna_pointing_rad(station_name, et_meas_vec, traj_eci);

    obs_predicted = zeros(n_epochs, n_meas_typ);
    obs_predicted(:, 1) = az_computed(:);
    obs_predicted(:, 2) = el_computed(:);
    if n_meas_typ == 3
        obs_predicted(:, 3) = rho_computed(:);
    end

    weighted_res = zeros(n_epochs, n_meas_typ);
    for ep = 1:n_epochs
        innov    = obs_predicted(ep, :) - meas_data(ep, :);
        innov(1) = angle(exp(1i * innov(1)));
        weighted_res(ep, :) = (W_matrix * innov.').';
    end

    res_vec = weighted_res(:);
end


function res_combined = get_multi_station_residuals(x_estimate, et_start, ...
        et_meas_troll, meas_data_troll, et_meas_kourou, meas_data_kourou, ...
        W_matrix, dynamic_model, mu)
% AIM: Concatenated weighted residuals from TROLL and KOUROU for dual-station OD.
% INPUTS:
%   x_estimate       : [6x1]     Current state estimate [km; km/s].
%   et_start         : [1x1]     Navigation epoch (integration origin) [s ET].
%   et_meas_troll    : [1xN_T]   TROLL observation epochs [s ET].
%   meas_data_troll  : [N_T x P] TROLL observations: Az [rad], El [rad], (Range [km]).
%   et_meas_kourou   : [1xN_K]   KOUROU observation epochs [s ET].
%   meas_data_kourou : [N_K x P] KOUROU observations: Az [rad], El [rad], (Range [km]).
%   W_matrix         : [PxP]     Measurement weight matrix (inv of noise sqrt).
%   dynamic_model    : [char]    Force model selector: 'Kepler' or 'J2'.
%   mu               : [1x1]     Earth gravitational parameter [km^3/s^2].
% OUTPUT:
%   res_combined     : [(N_T+N_K)*P x 1] Concatenated weighted residuals.

    earth_radii = cspice_bodvrd('EARTH', 'RADII', 3);
    Re_eq       = earth_radii(1);
    J2_val      = 0.0010826269;
    n_obs_types = size(meas_data_troll, 2);

    ode_settings = odeset('RelTol', 1e-12, 'AbsTol', [1e-6*ones(3,1); 1e-9*ones(3,1)]);

    if strcmp(dynamic_model, 'J2')
        eom_handle = @(t, state) keplerianJ2_rhs(t, state, mu, Re_eq, J2_val);
    else
        eom_handle = @(t, state) keplerian_rhs(t, state, mu);
    end

    % --- TROLL block ---
    n_troll   = numel(et_meas_troll);
    block_T   = zeros(n_troll, n_obs_types);
    if n_troll > 0
        tspan_troll   = [et_start; et_meas_troll(:)];
        [~, traj_T]   = ode113(eom_handle, tspan_troll, x_estimate(:), ode_settings);
        prop_troll     = traj_T(2:end, :).';

        [az_T, el_T, rho_T, ~] = antenna_pointing_rad('TROLL', et_meas_troll, prop_troll);
        pred_T = [az_T(:), el_T(:)];
        if n_obs_types == 3
            pred_T(:, 3) = rho_T(:); 
        end

        for ep = 1:n_troll
            innov    = pred_T(ep, :) - meas_data_troll(ep, :);
            innov(1) = angle(exp(1i * innov(1)));
            block_T(ep, :) = (W_matrix * innov.').';
        end
    end

    % --- KOUROU block ---
    n_kourou  = numel(et_meas_kourou);
    block_K   = zeros(n_kourou, n_obs_types);
    if n_kourou > 0
        tspan_kourou  = [et_start; et_meas_kourou(:)];
        [~, traj_K]   = ode113(eom_handle, tspan_kourou, x_estimate(:), ode_settings);
        prop_kourou    = traj_K(2:end, :).';

        [az_K, el_K, rho_K, ~] = antenna_pointing_rad('KOUROU', et_meas_kourou, prop_kourou);
        pred_K = [az_K(:), el_K(:)];
        if n_obs_types == 3
            pred_K(:, 3) = rho_K(:); 
        end

        for ep = 1:n_kourou
            innov    = pred_K(ep, :) - meas_data_kourou(ep, :);
            innov(1) = angle(exp(1i * innov(1)));
            block_K(ep, :) = (W_matrix * innov.').';
        end
    end

    res_combined = [block_T(:); block_K(:)];
end


function plot_od_residuals(case_name, res_init, res_final, W, R, ...
                            t_T, t_K, sigma_az, sigma_el, sigma_range)
% AIM: Visualizes the evolution of weighted residuals for OD diagnostics.
% INPUTS:
%   case_name   : [char] Label for the plot.
%   res_init    : [M x 1] Initial weighted residual vector.
%   res_final   : [M x 1] Final weighted residual vector.
%   W           : [PxP] Weighting matrix used in LSQ.
%   R           : [PxP] Covariance matrix (R = inv(W'*W)).
%   t_T / t_K   : [Nx1] Epochs for TROLL and KOUROU [s].
%   sigma_az/el : [1x1] Standard deviation for angles [rad].
%   sigma_range : [1x1] Standard deviation for range [km].

    m       = size(W, 1);
    sqrtR   = sqrtm(R);
    N_T = numel(t_T); N_K = numel(t_K);

    function uw = unweight(r_vec, N, m_sz, sqR)
        mat = reshape(r_vec(1:N*m_sz), m_sz, N);
        uw  = (sqR * mat).';
    end

    uw_T_i = unweight(res_init(1:N_T*m),  N_T, m, sqrtR);
    uw_T_f = unweight(res_final(1:N_T*m), N_T, m, sqrtR);

    if N_K > 0
        uw_K_i = unweight(res_init( N_T*m+1:end), N_K, m, sqrtR);
        uw_K_f = unweight(res_final(N_T*m+1:end), N_K, m, sqrtR);
    end

    uw_T_i(:,1:2) = uw_T_i(:,1:2) * 180/pi;
    uw_T_f(:,1:2) = uw_T_f(:,1:2) * 180/pi;
    if N_K > 0
        uw_K_i(:,1:2) = uw_K_i(:,1:2) * 180/pi;
        uw_K_f(:,1:2) = uw_K_f(:,1:2) * 180/pi;
    end

    sig_vals = [sigma_az*180/pi,  sigma_el*180/pi];
    if m == 3, sig_vals(3) = sigma_range; end

    c_T_init  = [1.00 0.55 0.20]; c_T_final = [0.80 0.15 0.05];
    c_K_init  = [0.95 0.75 0.10]; c_K_final = [0.00 0.45 0.74];

    titles = {'Azimuth Residuals','Elevation Residuals','Range Residuals'};
    ylbls  = {'[deg]','[deg]','[km]'};

    figure('Name',['Residuals: ',case_name]);
    for k = 1:m
        subplot(m, 1, k); hold on; grid on; box on;
        plot(t_T, uw_T_i(:,k), 'o', 'MarkerSize',6, 'MarkerEdgeColor',c_T_init,  'MarkerFaceColor','none', 'DisplayName','TROLL (init)');
        plot(t_T, uw_T_f(:,k), 'o', 'MarkerSize',6, 'MarkerEdgeColor',c_T_final, 'MarkerFaceColor',c_T_final, 'DisplayName','TROLL (final)');
        if N_K > 0
            plot(t_K, uw_K_i(:,k), 'o', 'MarkerSize',5, 'MarkerEdgeColor',c_K_init,  'MarkerFaceColor','none', 'DisplayName','KOUROU (init)');
            plot(t_K, uw_K_f(:,k), 'o', 'MarkerSize',5, 'MarkerEdgeColor',c_K_final, 'MarkerFaceColor',c_K_final, 'DisplayName','KOUROU (final)');
        end
        yline( 3*sig_vals(k), 'k--', 'LineWidth',1.2, 'DisplayName','$\pm 3\sigma$');
        yline(-3*sig_vals(k), 'k--', 'LineWidth',1.2, 'HandleVisibility','off');
        yline(0,            'k-',  'LineWidth',0.6, 'HandleVisibility','off');
        ylabel(ylbls{k}); title(titles{k});
        if k == 1, legend('Location','best'); end
        if k == m, xlabel('Ephemeris Time [s]'); end
    end
    % sgtitle(['Unweighted Residuals – ', case_name]);
end



function sigma_a = a_covariance_to_std(x_estimate, cov_matrix, mu)
% AIM: Maps state covariance to the standard deviation of the semi-major axis
%      via linear covariance propagation using a central finite-difference gradient.
% INPUTS:
%   x_estimate : [6x1] State estimate [km, km/s].
%   cov_matrix : [6x6] State covariance matrix [km^2, (km/s)^2].
%   mu         : [1x1] Gravitational parameter [km^3/s^2].
% OUTPUT:
%   sigma_a    : [1x1] Standard deviation of semi-major axis [km].

    n_state = numel(x_estimate);
    grad_a  = zeros(1, n_state);

    % Nominal semi-major axis via vis-viva
    r_nom = norm(x_estimate(1:3));
    v_nom = norm(x_estimate(4:6));
    a_nom = 1/(2/r_nom - v_nom^2/mu);

    % Central finite-difference gradient (Remark #1 step size)
    for k = 1:n_state
        h = sqrt(eps) * max(1, abs(x_estimate(k)));

        x_fwd = x_estimate; x_fwd(k) = x_fwd(k) + h;
        x_bwd = x_estimate; x_bwd(k) = x_bwd(k) - h;

        a_fwd = 1/(2/norm(x_fwd(1:3)) - norm(x_fwd(4:6))^2/mu);
        a_bwd = 1/(2/norm(x_bwd(1:3)) - norm(x_bwd(4:6))^2/mu);

        grad_a(k) = (a_fwd - a_bwd) / (2*h);
    end

    sigma_a = sqrt(abs(grad_a * cov_matrix * grad_a.'));
end