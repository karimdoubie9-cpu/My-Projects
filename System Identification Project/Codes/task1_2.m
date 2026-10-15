clear; clc; close all;

% --- COLOR CONFIGURATION (RGB) ---
palette = orderedcolors("gem");
c_raw   = [0.75, 0.75, 0.75]; % Gray (Raw FFT)
c_full  = palette(1,:); % Blue (Full Model)
c_welch = palette(2,:); % Orange (Welch Reference)
c_OE   = palette(3,:); % yellow (Manual Output Error Model)

fprintf('\n============================================================\n');
fprintf(' 1. TRAINING DATA LOADING\n');
fprintf('============================================================\n');

% Load training dataset with longitudinal dynamics signals

% DATA LOADED                       % Case for sensitivity analysis
load('data_sim.mat');           % Nominal case, equilibrated sensors
% load('data_sim_better_q.mat');  % Better pitch rate sensor
% load('data_sim_better_ax.mat'); % Better longitudinal acceleration sensor
if exist('ans', 'var'), sim_results = ans; end

Mtot = sim_results.Mtot; % Total pitching moment input
q    = sim_results.q;    % Pitch rate output
ax   = sim_results.ax;   % Longitudinal acceleration output
sample_time = 0.004;     % Sampling period [s]
t = 0:sample_time:(length(Mtot)-1)*sample_time;
N = length(t);

fprintf('Data loaded: %d samples | Ts: %.3f s\n', N, sample_time);

%% TIME-DOMAIN VISUALIZATION: SYSTEM SIGNALS
% Setup parameters for time-history plots
c_exc   = [0.85, 0.33, 0.10]; % Excitation color
c_input = [0.20, 0.20, 0.20]; % Total input color
c_out   = [0.00, 0.45, 0.74]; % Output signals color

figure('Name', 'System Signals (4-Plot)', 'Color', 'w');

% --- Subplot 1: Pure Excitation ---
subplot(2,2,1);
plot(ExcitationM(:,1),ExcitationM(:,2), 'Color', c_exc, 'LineWidth', 1.0);
ylabel(' [Nm]', 'Interpreter', 'latex');
xlabel('Time [s]', 'Interpreter', 'latex');
title('\textbf{1. Pure Excitation Signal}', 'Interpreter', 'latex');
xlim([t(1) t(end)]);
grid on; grid minor; set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');

% --- Subplot 2: Total Input (Control + Excitation) ---
subplot(2,2,2);
plot(t, Mtot, 'Color', c_input, 'LineWidth', 1.0);
grid on; grid minor; set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
ylabel(' [Nm]', 'Interpreter', 'latex');
xlabel('Time [s]', 'Interpreter', 'latex');
title('\textbf{2. Total Input}', 'Interpreter', 'latex');
xlim([t(1) t(end)]);

% --- Subplot 3: Output Pitch Rate ---
subplot(2,2,3);
plot(t, q, 'Color', c_out, 'LineWidth', 1.0);
grid on; grid minor; set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
ylabel(' [rad/s]', 'Interpreter', 'latex');
xlabel('Time [s]', 'Interpreter', 'latex');
title('\textbf{3. Output: Pitch Rate}', 'Interpreter', 'latex');
xlim([t(1) t(end)]);

% --- Subplot 4: Output Longitudinal Acceleration ---
subplot(2,2,4);
plot(t, ax, 'Color', c_out, 'LineWidth', 1.0);
grid on; grid minor; set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
ylabel(' [m/s$^2$]', 'Interpreter', 'latex');
xlabel('Time [s]', 'Interpreter', 'latex');
title('\textbf{4. Output: Longitudinal Acceleration}', 'Interpreter', 'latex');
xlim([t(1) t(end)]);

%% ========================================================================
%% 2. SPECTRAL ANALYSIS: COHERENCE EVALUATION
%% ========================================================================
fprintf('\n============================================================\n');
fprintf(' 2. SPECTRAL ANALYSIS (COHERENCE)\n');
fprintf('============================================================\n');

% Frequency resolution and DFT point selection
Nf = 2^ceil(log2(N));
[~, f_full] = dft(Mtot,t,Nf);

% Welch Method configuration: 16 segments for variance reduction
Nwin = 16;
lwindow = round(N/Nwin);

% Computation of Power Spectral Densities (PSD) and Cross-Spectral Densities (CSD)
Suu = cpsd(Mtot,Mtot,lwindow,[],f_full,1/sample_time);
Suq = cpsd(q,Mtot,lwindow,[],f_full,1/sample_time);
Sua = cpsd(ax,Mtot,lwindow,[],f_full,1/sample_time);
Sqq = cpsd(q,q,lwindow,[],f_full,1/sample_time);
Saa = cpsd(ax,ax,lwindow,[],f_full,1/sample_time);

idx_pos = f_full>=0;
f = f_full(idx_pos);

% Squared Coherence Function estimation to assess data reliability
gamma_q  = abs(Suq(idx_pos)).^2 ./ (Suu(idx_pos).*Sqq(idx_pos));
gamma_ax = abs(Sua(idx_pos)).^2 ./ (Suu(idx_pos).*Saa(idx_pos));

% Plot Coherence results
figure('Name','Coherence', 'Color', 'w', 'Position', [100 100 800 400]);
subplot(1,2,1);
plot(f, gamma_q, 'Color', c_full, 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'FontSize', 12); grid on; grid minor; ylim([0 1.1]);
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('Coherence  [-]', 'Interpreter', 'latex');
title('\textbf{Coherence: Pitch Rate ($q$)}', 'Interpreter', 'latex');

subplot(1,2,2);
plot(f, gamma_ax, 'Color', c_full, 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'FontSize', 12); grid on; grid minor; ylim([0 1.1]);
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('Coherence  [-]', 'Interpreter', 'latex');
title('\textbf{Coherence: Acceleration ($a_x$)}', 'Interpreter', 'latex');

%% ========================================================================
%% 3. GREY-BOX PARAMETER ESTIMATION
%% ========================================================================
fprintf('\n============================================================\n');
fprintf(' 3. GREY-BOX IDENTIFICATION\n');
fprintf('============================================================\n');

% Selection of frequency range with high coherence for identification
f_range_idx = f > 0.76 & f < 1.85;
f_ident = f(f_range_idx);

% Frequency-domain data preparation (DFT)
[Mtot_f,~] = dft(Mtot,t,Nf);
[q_f,~]    = dft(q,t,Nf);
[ax_f,~]   = dft(ax,t,Nf);

Mtot_f = Mtot_f(idx_pos); Mtot_f = Mtot_f(f_range_idx);
q_f    = q_f(idx_pos);    q_f    = q_f(f_range_idx);
ax_f   = ax_f(idx_pos);   ax_f   = ax_f(f_range_idx);

% Create identification data object (Frequency domain)
data_id = iddata([q_f ax_f], Mtot_f, sample_time, ...
'Frequency',f_ident,'FrequencyUnit','Hz');
data_id.InputName = 'Mtot';
data_id.OutputName = {'q','ax'};

% Initialize Grey-Box model with physical aerodynamic derivatives
ssmodel = @UAVmodel;
params = {'Xu',0;'Xq',0;'Mu',0;'Mq',0;'Xd',0;'Md',0};
sys_ini = idgrey(ssmodel,params,'c');

% Use greyest
opt = greyestOptions('Display','off','SearchMethod','lm');
sys_grey = greyest(data_id,sys_ini,opt);

% --- Statistical Uncertainty Report ---
param_vals = sys_grey.Report.Parameters.ParVector;
cov_mat = getcov(sys_grey);
if isempty(cov_mat), cov_mat = zeros(length(param_vals)); end
stddev = sqrt(diag(cov_mat));
unc_rel = (stddev ./ abs(param_vals)) * 100;
param_names = {'Xu', 'Xq', 'Mu', 'Mq', 'Xd', 'Md'};

fprintf('%-10s | %-12s | %-15s\n', 'PARAMETER', 'VALUE', 'UNCERTAINTY (%)');
fprintf('--------------------------------------------\n');
for i = 1:length(param_vals)
fprintf('%-10s | %12.4f | %14.2f %%\n', param_names{i}, param_vals(i), unc_rel(i));
end
fprintf('--------------------------------------------\n');

%% ========================================================================
%% 4. MANUAL OUTPUT ERROR (OE) IDENTIFICATION (FULL FREQUENCY DOMAIN)
%% ========================================================================
fprintf('\n============================================================\n');
fprintf(' 4. MANUAL OUTPUT ERROR (OE) IDENTIFICATION (FULL FREQUENCY DOMAIN)\n');
fprintf('============================================================\n');


% OE method

% sampling frequency [Hz]
fs = 1/sample_time;

% Subtracting means from original time history data
%filtering
fc = 10;                % Cutoff frequency [Hz]
fn = fs / 2;           % Nyquist frequency
order = 1;              % 1st Order

% Design the filter coefficients (Butterworth)
[b_lpf, a_lpf] = butter(order, fc/fn);
ax_pre_filter=ax;

% Apply Zero-Phase filtering
% We overwrite 'ax' so the rest of the script uses the filtered version
ax = filtfilt(b_lpf, a_lpf, ax);

x_withoutmean = Mtot-mean(Mtot);
y1_withoutmean = q-mean(q);
y2_withoutmean = ax-mean(ax);

% grid on;
nx = length(x_withoutmean);
na = 12; % Averaging factor
window = hanning(floor(nx/na));
noverlap = floor(length(window)/2); % 50% overlap
nfft = 2048;

% 2.3 Compute Spectra
[Suu, f_hz] = pwelch(x_withoutmean, window, noverlap, nfft, fs);
[Suq, ~] = cpsd(y1_withoutmean, x_withoutmean, window, noverlap, nfft, fs);
[Sqq, ~] = pwelch(y1_withoutmean, window, noverlap, nfft, fs);
[Suax, ~] = cpsd(y2_withoutmean, x_withoutmean, window, noverlap, nfft, fs);
[Saxax, ~] = pwelch(y2_withoutmean, window, noverlap, nfft, fs);

% Convert frequency to rad/s for plotting consistency with reference
f_rads = f_hz * 2 * pi;

% 2.4 Compute FRF
Gq_np = Suq ./ Suu;   % Non-parametric FRF for pitch rate
Gax_np = Suax ./ Suu; % Non-parametric FRF for acceleration

% 2.5 Compute Coherence for Quality Assessment
gamma_sq_q = (abs(Suq).^2) ./ (Suu .* Sqq);
gamma_sq_ax = (abs(Suax).^2) ./ (Suu .* Saxax);


theta_true = [ -0.1068; 0.1192; -5.9755; -2.6478; -10.1647; 450.71 ];
sys_true = idgrey('uav_longitudinal_grey', theta_true, 'c');




% 4. Data Preparation for Identification
omega = f_rads;

% --- MODIFIED: Do not limit frequency domain, use the full spectrum > 0.1 rad/s ---
valid_idx = find(omega > 2 & omega < 12);

% Extract data for the valid range
omega_id = omega(valid_idx);
H_q_id   = Gq_np(valid_idx);   % Measured FRF for Pitch Rate
H_ax_id  = Gax_np(valid_idx);  % Measured FRF for Acceleration (optional validation)

% 6. Run Output Error Identification
disp('Starting Identification...');


guess = 2*[Xu Xq Mu Mq Xd Md] + 0.1 * (-0.5 + rand(1,6)).*[Xu Xq Mu Mq Xd Md];
G_th = @(theta) get_drone_model(theta); % SIMO TF


%
% Ensure Inputs are Column Vectors first [K x 1]
H_q_id = H_q_id(:);
H_ax_id = H_ax_id(:);

% Create [2 x K] Matrix using DOT TRANSPOSE (.')
Y_m = [H_q_id.'; H_ax_id.'];

R = diag([(noise.ang_rate_stand_dev)^2 (noise.acc_stand_dev)^2]);
[theta_est, Cov_theta, J_opt] = outputError_FD2(Y_m, omega_id, G_th, guess',R);

sys_OE = idgrey('uav_longitudinal_grey', theta_est, 'c');

unc_rel = 100 * sqrt(diag(Cov_theta)) ./ abs(theta_est);
disp('Identification complete.')
disp('Identified Parameters:');
fprintf('%-10s | %-12s | %-15s\n', 'PARAMETER', 'VALUE', 'UNCERTAINTY (%)');
fprintf('--------------------------------------------\n');
for i = 1:length(param_vals)
fprintf('%-10s | %12.4f | %14.2f %%\n', param_names{i}, theta_est(i), unc_rel(i));
end
fprintf('--------------------------------------------\n');



% Comparative FIT analysis within the identification range
[~, fit_full_vec] = compare(data_id, sys_grey);
[~, fit_OE_vec]  = compare(data_id, sys_OE);
fit_q_full  = fit_full_vec(1);
fit_ax_full = fit_full_vec(2);
fit_q_OE   = fit_OE_vec(1);
fit_ax_OE  = fit_OE_vec(2);


%% ========================================================================
%% 6. CROSS-VALIDATION WITH INDEPENDENT DATASET
%% ========================================================================
fprintf('\n============================================================\n');
fprintf(' 6. VALIDATION (Using iddata/compare)\n');
fprintf('============================================================\n');

% Load 3-2-1-1 validation dataset

% DATA LOADED                            % Case for sensitivity analysis
load('data_sim_3211.mat');           % Nominal case, equilibrated sensors
% load('data_sim_3211_better_q.mat');  % Better pitch rate sensor
% load('data_sim_3211_better_ax.mat'); % Better longitudinal acceleration sensor
if exist('ans', 'var'), simv = ans; end

Mtot_v = simv.Mtot;
q_v    = simv.q;
ax_v   = simv.ax;
t_v = 0:sample_time:(length(Mtot_v)-1)*sample_time;
Nv = length(t_v);
Nf_v = 2^ceil(log2(Nv));

% Preparation of Validation Data in Frequency Domain
[~,f_full_v] = dft(Mtot_v,t_v,Nf_v);
idx_pos_v = f_full_v>=0;
f_v = f_full_v(idx_pos_v);

% Define validation frequency range
f_range_idx_v = f_v>0.3 & f_v<1.7;
f_ident_v = f_v(f_range_idx_v);

[Mtot_fv,~] = dft(Mtot_v,t_v,Nf_v);
[q_fv,~]    = dft(q_v,t_v,Nf_v);
[ax_fv,~]   = dft(ax_v,t_v,Nf_v);

Mtot_fv = Mtot_fv(idx_pos_v); Mtot_fv = Mtot_fv(f_range_idx_v);
q_fv    = q_fv(idx_pos_v);    q_fv    = q_fv(f_range_idx_v);
ax_fv   = ax_fv(idx_pos_v);   ax_fv   = ax_fv(f_range_idx_v);

data_val = iddata([q_fv ax_fv], Mtot_fv, sample_time, ...
'Frequency',f_ident_v,'FrequencyUnit','Hz');

% Calculate goodness-of-fit (FIT) for validation
[~, fit_full_vec] = compare(data_val, sys_grey);
[~, fit_OE_vec]  = compare(data_val, sys_OE);
fit_q_full  = fit_full_vec(1);
fit_ax_full = fit_full_vec(2);
fit_q_OE   = fit_OE_vec(1);
fit_ax_OE  = fit_OE_vec(2);


fprintf('%-15s | %-12s | %-12s\n', 'OUTPUT', 'GREYEST', 'MANUAL OE');
fprintf('----------------|--------------|------------\n');
fprintf('%-15s | %10.2f %% | %8.2f %%\n', 'Pitch Rate (q)', fit_q_full, fit_q_OE);
fprintf('%-15s | %10.2f %% | %8.2f %%\n', 'Accel (ax)', fit_ax_full, fit_ax_OE);
fprintf('----------------|--------------|------------\n');


%% ========================================================================
%% 7. FREQUENCY RESPONSE FUNCTION (FRF) DATA PREPARATION
%% ========================================================================

% Empirical Transfer Function Estimation (ETFE) via Welch for comparison
window = round(Nv/8);
noverlap = round(window/2);
[Hq_w,f_w]  = tfestimate(Mtot_v,q_v,window,noverlap,Nf_v,1/sample_time);
[Hax_w,~]   = tfestimate(Mtot_v,ax_v,window,noverlap,Nf_v,1/sample_time);

idx_fit = f_w>0.3 & f_w<1.7;
f_plot = f_w(idx_fit);

% Evaluate model frequency responses
resp_full = freqresp(sys_grey, f_plot, 'Hz');
resp_OE  = freqresp(sys_OE, f_plot, 'Hz');
Gq_full  = squeeze(resp_full(1,1,:));
Gax_full = squeeze(resp_full(2,1,:));
Gq_OE   = squeeze(resp_OE(1,1,:));
Gax_OE  = squeeze(resp_OE(2,1,:));

% Raw FFT processing for full-spectrum visualization
F_raw = (0:Nf_v-1)'*(1/sample_time)/Nf_v;
U_fft = fft(Mtot_v,Nf_v);
Q_fft = fft(q_v,Nf_v);
A_fft = fft(ax_v,Nf_v);
Hq_raw  = Q_fft ./ U_fft;
Hax_raw = A_fft ./ U_fft;

% Indices for extended plotting range
idx_plot_raw   = F_raw > 0.1 & F_raw < 100;
idx_plot_welch = f_w   > 0.1 & f_w   < 100;
f_raw_p = F_raw(idx_plot_raw);
f_w_p   = f_w(idx_plot_welch);

% Extended frequency response evaluation
resp_full_ext = freqresp(sys_grey, f_w_p, 'Hz');
resp_OE_ext  = freqresp(sys_OE, f_w_p, 'Hz');
Gq_full_ext  = squeeze(resp_full_ext(1,1,:));
Gax_full_ext = squeeze(resp_full_ext(2,1,:));
Gq_OE_ext   = squeeze(resp_OE_ext(1,1,:));
Gax_OE_ext  = squeeze(resp_OE_ext(2,1,:));

% FRF of the true system
resp_true_ext  = freqresp(sys_true, f_w_p, 'Hz');
Gq_true_ext  = squeeze(resp_full_ext(1,1,:));
Gax_true_ext = squeeze(resp_full_ext(2,1,:));

%% ========================================================================
%% 8. VALIDATION PLOTS (BODE MAGNITUDE & PHASE)
%% ========================================================================

% --- VALIDATION PLOT: PITCH RATE () ---
figure('Name', 'Validation - q', 'Color', 'w', 'Position', [100 100 800 600])
subplot(2,1,1)
plot(f_raw_p, 20*log10(abs(Hq_raw(idx_plot_raw))), 'Color', c_raw, 'LineWidth', 0.5); hold on
plot(f_w_p,   20*log10(abs(Hq_w(idx_plot_welch))), 'Color', c_welch, 'LineWidth', 2.0)
plot(f_w_p,   20*log10(abs(Gq_true_ext)), '-',    'Color', 'red',  'LineWidth', 2)
plot(f_w_p,   20*log10(abs(Gq_full_ext)), '--',    'Color', c_full,  'LineWidth', 2)
plot(f_w_p,   20*log10(abs(Gq_OE_ext)),  '-.',    'Color', c_OE,   'LineWidth', 2)
grid on; grid minor; set(gca, 'XScale', 'log', 'FontSize', 12); xlim([0.1 100]);
ylabel('Magnitude [dB]', 'Interpreter', 'latex');
legend('Raw FFT', 'Welch (Ref)', 'True system', 'Greyest', 'Manual Output Error', 'Location','SouthWest', 'Interpreter','latex');
title(['\textbf{Pitch Rate ($q$) - best FIT.: ', num2str(fit_q_full, '%.2f'), ' \%}'], 'Interpreter', 'latex');

subplot(2,1,2)
plot(f_raw_p, angle(Hq_raw(idx_plot_raw))*180/pi, 'Color', c_raw, 'LineWidth', 0.5); hold on
plot(f_w_p,   angle(Hq_w(idx_plot_welch))*180/pi, 'Color', c_welch, 'LineWidth', 2.0)
plot(f_w_p,   angle(Gq_true_ext)*180/pi, '-',    'Color', 'red',  'LineWidth', 2)
plot(f_w_p,   angle(Gq_full_ext)*180/pi, '--',    'Color', c_full,  'LineWidth', 2)
plot(f_w_p,   angle(Gq_OE_ext)*180/pi,  '-.',    'Color', c_OE,   'LineWidth', 2)
grid on; grid minor; set(gca, 'XScale', 'log', 'FontSize', 12); xlim([0.1 100]);
ylabel('Phase [deg]', 'Interpreter', 'latex'); xlabel('Frequency [Hz]', 'Interpreter', 'latex');
yticks(-180:90:180); ylim([-200 200]);

% --- VALIDATION PLOT: LONGITUDINAL ACCELERATION () ---
figure('Name', 'Validation - ax', 'Color', 'w', 'Position', [150 150 800 600])
subplot(2,1,1)
plot(f_raw_p, 20*log10(abs(Hax_raw(idx_plot_raw))), 'Color', c_raw, 'LineWidth', 0.5); hold on
plot(f_w_p,   20*log10(abs(Hax_w(idx_plot_welch))), 'Color', c_welch, 'LineWidth', 2.0)
plot(f_w_p,   20*log10(abs(Gax_true_ext)), '-',    'Color', 'red',  'LineWidth', 2)
plot(f_w_p,   20*log10(abs(Gax_full_ext)), '--',    'Color', c_full,  'LineWidth', 2)
plot(f_w_p,   20*log10(abs(Gax_OE_ext)),  '-.',    'Color', c_OE,   'LineWidth', 2)
grid on; grid minor; set(gca, 'XScale', 'log', 'FontSize', 12); xlim([0.1 100]);
ylabel('Magnitude [dB]', 'Interpreter', 'latex');
legend('Raw FFT',  'Welch (Ref)', 'True system', 'Greyest', 'Manual Output Error', 'Location','SouthWest', 'Interpreter','latex');
title(['\textbf{Acceleration ($a_x$) - best FIT.: ', num2str(fit_ax_full, '%.2f'), ' \%}'], 'Interpreter', 'latex');

subplot(2,1,2)
plot(f_raw_p, angle(Hax_raw(idx_plot_raw))*180/pi, 'Color', c_raw, 'LineWidth', 0.5); hold on
plot(f_w_p,   angle(Gax_true_ext)*180/pi, '-',    'Color', 'red',  'LineWidth', 2)
plot(f_w_p,   angle(Hax_w(idx_plot_welch))*180/pi, 'Color', c_welch, 'LineWidth', 2.0)
plot(f_w_p,   angle(Gax_full_ext)*180/pi, '--',    'Color', c_full,  'LineWidth', 2)
plot(f_w_p,   angle(Gax_OE_ext)*180/pi,  '-.',    'Color', c_OE,   'LineWidth', 2)
grid on; grid minor; set(gca, 'XScale', 'log', 'FontSize', 12); xlim([0.1 100]);
ylabel('Phase [deg]', 'Interpreter', 'latex'); xlabel('Frequency [Hz]', 'Interpreter', 'latex');
yticks(-180:90:180); ylim([-200 200]);

fprintf('\n============================================================\n');
fprintf(' ANALYSIS COMPLETED\n');
fprintf('============================================================\n');


%% AUXILIARY FUNCTIONS

function [G,f] = dft(g,t,Nf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCREET FOURIER TRANSFORM
% g = input signal
% t = input time axis
% Nf = number of frequency samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(exist('Nf'))
    Nf = length(t);
end

if size(g,1)==1
    g = g.';
end

G = fft(g,Nf,1); % transform along columns by default
G = fftshift(G,1);

% frequency axis
dt = t(2)-t(1);
if mod(Nf,2) % odd number of samples
    f = (-(Nf-1)/2:(Nf-1)/2)/Nf/dt;
else   % even number of samples
    f = (-Nf/2:Nf/2-1)/Nf/dt;
end

% fft delay compensation
phase = -2*pi*f(:)*t(1);
w = exp(1i*phase(:))*ones(1,size(g,2));
G = G.*w;
end