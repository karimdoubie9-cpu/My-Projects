%% ========================================================================
%  STAGE 2: STATISTICAL ANALYSIS AND PLOTTING (Task 3.2)
%  ========================================================================
close all
clear 
clc
% Load Monte Carlo results
load('MC_results_Perturbed1.mat')

% Extract variables from saved structure
theta_hat_mc   = MC_results.theta_hat_mc;
sigma_theta=MC_results.sigma_theta_mc;
results        = MC_results.results_struct;
numSimulations = MC_results.numSimulations;
% Empirical standard deviation from MC runs
std_emp = std(theta_hat_mc, 0, 2);          % 6x1, sample std across runs
mean_theta = mean(theta_hat_mc, 2);          % 6x1, MC mean

% Relative std as %, normalized by MC mean
U_rel_emp = 100 * std_emp ./ abs(mean_theta)
%%
% Define parameter names for plotting with LaTeX formatting
param_names_plain = {'Xu', 'Xq', 'Mu', 'Mq', 'X_delta', 'M_delta'};
param_names_latex = {'$X_u$ $[1/s]$', '$X_q$ $[m/s]$', '$M_u$ $[rads/m]$', ...
                     '$M_q$ $[1/s]$', '$X_{\delta}$ $[m/s^2]$', '$M_{\delta}$ $[rad/s^2]$'};
param_indices = [1, 2, 3, 4, 5, 6]; % Indices for all 6 parameters

% Get the true values for plotting the reference line
theta_true=[ -0.1068; 0.1192; -5.9755; -2.6478; -10.164; 450.0 ];
% -> Step #2.1: Histograms of Estimated Parameters
    figure('Name', 'Monte Carlo Parameter Distributions');
    sgtitle('Monte Carlo Parameter Distributions', 'FontWeight', 'bold');
    
    for i = 1:6
        j = param_indices(i); 
        subplot(3, 2, i);
        hold on;
        
        histogram(theta_hat_mc(j, :), 'Normalization', 'pdf', ...
                  'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
        
        mean_val = mean(theta_hat_mc(j, :));
        xline(mean_val, '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
        
        title(param_names_latex{i}, 'Interpreter', 'latex');
        xlabel(param_names_latex{i}, 'Interpreter', 'latex');
        ylabel('Probability Density');
        grid on;
        hold off;
      end


% -> Step #2.3: Histograms of Pole Distribution
    all_poles_mc = zeros(numSimulations, 3);
    
    for i = 1:numSimulations
        sys_i = idgrey('uav_longitudinal_grey', theta_hat_mc(:,i), 'c');
        poles = sort(pole(sys_i), 'descend', 'ComparisonMethod', 'abs');
        all_poles_mc(i, :) = poles.';
    end
    
    figure('Name', 'Monte Carlo Pole Distributions');
    sgtitle('Monte Carlo Pole Distributions', 'FontWeight', 'bold');
    
    for i = 1:3
        subplot(3, 2, 2*i - 1);
        hold on;
        real_parts = real(all_poles_mc(:, i));
        histogram(real_parts, 'Normalization', 'pdf', ...
                  'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
        xline(mean(real_parts), '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
        grid on;
        title(['Re($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        xlabel(['Re($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        ylabel('Probability Density');
        hold off;
    
        subplot(3, 2, 2*i);
        hold on;
        imag_parts = imag(all_poles_mc(:, i));
        histogram(imag_parts, 'Normalization', 'pdf', ...
                  'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
        xline(mean(imag_parts), '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
        grid on;
        title(['Im($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        xlabel(['Im($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        ylabel('Probability Density');
        hold off;
    end

% -> Step #2.4: Histograms of Zero Distribution
    all_zeros_q_mc = zeros(numSimulations, 2);
    all_zeros_ax_mc = zeros(numSimulations, 3);
    
    for i = 1:numSimulations
        sys_i = idgrey('uav_longitudinal_grey', theta_hat_mc(:,i), 'c');
        zeros_q = sort(zpk(sys_i).Z{1}, 'descend', 'ComparisonMethod', 'abs');
        zeros_ax = sort(zpk(sys_i).Z{2}, 'descend', 'ComparisonMethod', 'abs');
        all_zeros_q_mc(i, :) = zeros_q.';
        all_zeros_ax_mc(i, :) = zeros_ax.';
    end
    
    figure('Name', 'Monte Carlo Zero Distributions (q)');
    sgtitle('Monte Carlo Zero Distributions (q)', 'FontWeight', 'bold');
    
    for i = 1:2
        subplot(2, 2, 2*i - 1);
        hold on;
        real_parts = real(all_zeros_q_mc(:, i));
        histogram(real_parts, 'Normalization', 'pdf', ...
                  'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
        xline(mean(real_parts), '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
        grid on;
        title(['Re($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        xlabel(['Re($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        ylabel('Probability Density');
        hold off;
    
        subplot(2, 2, 2*i);
        hold on;
        imag_parts = imag(all_zeros_q_mc(:, i));
        histogram(imag_parts, 'Normalization', 'pdf', ...
                  'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
        xline(mean(imag_parts), '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
        grid on;
        title(['Im($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        xlabel(['Im($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        ylabel('Probability Density');
        hold off;
    end

    figure('Name', 'Monte Carlo Zero Distributions (ax)');
    sgtitle('Monte Carlo Zero Distributions (ax)', 'FontWeight', 'bold');
    
    for i = 1:3
        subplot(3, 2, 2*i - 1);
        hold on;
        real_parts = real(all_zeros_ax_mc(:, i));
        histogram(real_parts, 'Normalization', 'pdf', ...
                  'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
        xline(mean(real_parts), '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
        grid on;
        title(['Re($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        xlabel(['Re($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        ylabel('Probability Density');
        hold off;
    
        subplot(3, 2, 2*i);
        hold on;
        imag_parts = imag(all_zeros_ax_mc(:, i));
        histogram(imag_parts, 'Normalization', 'pdf', ...
                  'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
        xline(mean(imag_parts), '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
        grid on;
        title(['Im($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        xlabel(['Im($p_', num2str(i), '$)'], 'Interpreter', 'latex');
        ylabel('Probability Density');
        hold off;
    end

%%
%->Step#2.3: Convergence of the Means (Running average) 
    k = 1:numSimulations;
    cumulative_Mean_theta = cumsum(theta_hat_mc, 2) ./ k;
    
    figure('Name', 'Monte Carlo Convergence of the Mean');
    sgtitle('Convergence of the Mean of Parameter Estimates', 'FontWeight', 'bold');
    
    for j = 1:6
        subplot(3, 2, j);
        plot(k, cumulative_Mean_theta(j,:), 'LineWidth', 1.5);
        title(['Mean of ', param_names_plain{j}]);
        xlabel('Number of MC runs');
        grid on;
        xlim([1, numSimulations]);
    end

%->Step#2.4: Dispersion of poles and zero
    % Store all poles and zeros from all N models
    all_poles = [];
    all_zeros_q = [];
    all_zeros_ax = [];
    
    for i = 1:numSimulations
        % Recreate the idgrey object for the i-th parameter set
        sys_i = idgrey('uav_longitudinal_grey', theta_hat_mc(:,i), 'c');
        % Get poles and zeros
        poles = pole(sys_i);
        zeros_q = zpk(sys_i).Z{1};
        zeros_ax = zpk(sys_i).Z{2};

        all_poles = [all_poles; poles(:)];
        all_zeros_q = [all_zeros_q; zeros_q(:)];
        all_zeros_ax = [all_zeros_ax; zeros_ax(:)];
    end

    % poles plot
    figure('Name', 'Monte Carlo Pole Dispersion');
    hold on; grid on;

    % Plot the poles and zeros from all runs
    plot(real(all_poles), imag(all_poles), 'rx', 'MarkerSize', 8, 'DisplayName', 'Estimated Poles');
    
    % For comparison, plot the true poles and zeros
    sys_true = idgrey('uav_longitudinal_grey', theta_true, 'c');
    true_poles = pole(sys_true);
    plot(real(true_poles), imag(true_poles), 'b*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'True Poles');

    % Plot dashed grey real and imaginary axes lines through origin
    xline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'HandleVisibility','off');   % real axis
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'HandleVisibility','off');   % imaginary axis

    legend('Location', 'best');
    title('Dispersion of Poles');
    xlabel('Real Axis');
    ylabel('Imaginary Axis');
    axis equal;


    % zeros plot
    figure('Name', 'Monte Carlo Zero Dispersion');
    hold on; grid on;

    % Plot the poles and zeros from all runs
    plot(real(all_zeros_q), imag(all_zeros_q), 'rx', 'MarkerSize', 8, 'DisplayName', 'Estimated q Zeros');
    plot(real(all_zeros_ax), imag(all_zeros_ax), 'gx', 'MarkerSize', 8, 'DisplayName', 'Estimated ax Zeros');
    
    % For comparison, plot the true poles and zeros
    true_zeros_q = zpk(sys_true).Z{1};
    true_zeros_ax = zpk(sys_true).Z{2};
    plot(real(true_zeros_q), imag(true_zeros_q), 'b*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'True zeros');
    plot(real(true_zeros_ax), imag(true_zeros_ax), 'b*', 'MarkerSize', 12, 'LineWidth', 2, 'HandleVisibility','off');

    % Plot dashed grey real and imaginary axes lines through origin
    xline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'HandleVisibility','off');   % real axis
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'HandleVisibility','off');   % imaginary axis

    legend('Location', 'best');
    title('Dispersion of Zeros');
    xlabel('Real Axis');
    ylabel('Imaginary Axis');
    axis equal;

%->Step#2.5: Dispersion of FRF 
    figure('Name', 'Monte Carlo FRF Dispersion');
    sgtitle('FRF Dispersion', 'FontWeight', 'bold');
    w_plot = logspace(-1, 2, 500); % Frequency vector for Bode plots [rad/s]
    
    % --- Calculate Bode Responses for a subset of models to keep plot clean ---
    num_plot_frfs = min(50, numSimulations); % Plot up to 50 random FRFs
    rand_indices = randperm(numSimulations, num_plot_frfs);
    bode_mag_q  = zeros(length(w_plot), num_plot_frfs);
    bode_phase_q = zeros(length(w_plot), num_plot_frfs);
    bode_mag_ax  = zeros(length(w_plot), num_plot_frfs);
    bode_phase_ax = zeros(length(w_plot), num_plot_frfs);
    
    for i = 1:num_plot_frfs
        idx = rand_indices(i);
        sys_i = idgrey('uav_longitudinal_grey', theta_hat_mc(:,idx), 'c');
        [mag, ph] = bode(sys_i, w_plot);
        
        bode_mag_q(:,i)  = 20*log10(squeeze(mag(1,1,:)));
        bode_phase_q(:,i)= squeeze(ph(1,1,:));
        bode_mag_ax(:,i) = 20*log10(squeeze(mag(2,1,:)));
        bode_phase_ax(:,i)= squeeze(ph(2,1,:));
    end
    
    % --- Calculate Nominal and True Model Responses ---
    theta_nominal = mean(theta_hat_mc, 2);
    sys_nominal = idgrey('uav_longitudinal_grey', theta_nominal, 'c');
    sys_true = idgrey('uav_longitudinal_grey', theta_true, 'c');
    
    [mag_nom, ph_nom]   = bode(sys_nominal, w_plot);
    [mag_true, ph_true] = bode(sys_true, w_plot);
    
    % --- Plot Dispersion for Pitch Rate (q) ---
    subplot(2,2,1);
    plot(w_plot, bode_mag_q, 'Color', [0.8 0.8 0.8]); hold on;
    h_nom = plot(w_plot, 20*log10(squeeze(mag_nom(1,1,:))), 'b', 'LineWidth', 2);
    h_true = plot(w_plot, 20*log10(squeeze(mag_true(1,1,:))), 'r--', 'LineWidth', 2);
    set(gca, 'XScale', 'log'); grid on;
    title('FRF Dispersion (q)'); ylabel('Magnitude [dB]');
    legend([h_nom, h_true], 'Nominal Model', 'True System');
    
    subplot(2,2,3);
    plot(w_plot, bode_phase_q, 'Color', [0.8 0.8 0.8]); hold on;
    plot(w_plot, squeeze(ph_nom(1,1,:)), 'b', 'LineWidth', 2);
    plot(w_plot, squeeze(ph_true(1,1,:)), 'r--', 'LineWidth', 2);
    set(gca, 'XScale', 'log'); grid on;
    xlabel('Frequency [rad/s]'); ylabel('Phase [deg]');
    
    % --- Plot Dispersion for Acceleration (ax) ---
    subplot(2,2,2);
    plot(w_plot, bode_mag_ax, 'Color', [0.8 0.8 0.8]); hold on;
    plot(w_plot, 20*log10(squeeze(mag_nom(2,1,:))), 'b', 'LineWidth', 2);
    plot(w_plot, 20*log10(squeeze(mag_true(2,1,:))), 'r--', 'LineWidth', 2);
    set(gca, 'XScale', 'log'); grid on;
    title('FRF Dispersion (a_x)');
    
    subplot(2,2,4);
    plot(w_plot, bode_phase_ax, 'Color', [0.8 0.8 0.8]); hold on;
    plot(w_plot, squeeze(ph_nom(2,1,:)), 'b', 'LineWidth', 2);
    plot(w_plot, squeeze(ph_true(2,1,:)), 'r--', 'LineWidth', 2);
    set(gca, 'XScale', 'log'); grid on;
    xlabel('Frequency [rad/s]');
    
    fprintf('\n--- End of Statistical Analysis ---\n');