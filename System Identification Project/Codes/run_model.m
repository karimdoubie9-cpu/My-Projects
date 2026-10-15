function estimated_model = run_model(out)

%Purpose: This function encapsulates a single, complete grey-box identification trial.
        % It takes the raw output from a Simulink simulation, performs all necessary data 
        % processing and estimation steps, and returns the final identified model.
   
%Input:  
    % Type: A Simulink.SimulationOutput object.
    % Structure: This object contains the results logged from the sim command. The function specifically expects it to have the following fields, which are timeseries objects:
                % out.Mtot1: The total pitching moment time-series.
                % out.q1: The pitch rate time-series.
                % out.ax1: The longitudinal acceleration time-series
                % run_model - Performs a single grey-box identification run.
                % This function takes the raw output structure from a Simulink simulation,
                % processes the data, and returns the estimated idgrey model.

%Output: estimated_model
% Type: An idgrey model object from the System Identification Toolbox.
% Content: This object represents the final identified grey-box model of the UAV's longitudinal dynamics. It contains:
        % The estimated values for the six physical parameters (Xu, Xq, Mu, Mq, X_delta, M_delta).
        % The associated parameter covariance matrix, which quantifies the
        % uncertainty of the estimate
     
% --- Step 1: Data Extraction ---
    % t = out.q1.Time;
    % Ts = t(2) - t(1);
    % fs = 1/Ts;
    % u_raw = out.Mtot1.Data;
    % q_raw = out.q1.Data;
    % ax_raw = out.ax1.Data;

    t = out.q.Time;
    Ts = t(2) - t(1);
    fs = 1/Ts;
    u_raw = out.Mtot.Data;
    q_raw = out.q.Data;
    ax_raw = out.ax.Data;



% --- Step 2: Frequency-Domain Data Preparation (Using RAW FFT) ---

% Compute the full FFT of the zero-mean time-series signals
u_FFT = fft(u_raw - mean(u_raw));
q_FFT = fft(q_raw - mean(q_raw));
ax_FFT = fft(ax_raw - mean(ax_raw));

% Create the corresponding one-sided frequency vector
f_hz = linspace(0, fs/2, floor(length(ax_FFT)/2) + 1);
f_rads = 2 * pi * f_hz;

% Select the frequency range of interest for identification
idx_id = (f_rads >= 5) & (f_rads <= 12);
omega_id = f_rads(idx_id);

% Prepare the FFT data for the 'iddata' object
% Note: The first half of the FFT is used for real signals.
N_full = length(u_FFT);
Y_full = [q_FFT(1:floor(N_full/2)+1), ax_FFT(1:floor(N_full/2)+1)];
U_full = u_FFT(1:floor(N_full/2)+1);

Y = [Y_full(idx_id, 1).'; Y_full(idx_id, 2).']; % Correctly form the Y matrix

% Create the iddata object
data = iddata(Y.', U_full(idx_id), Ts, 'Frequency', omega_id);
data.InputName = {'delta'};
data.OutputName = {'Angular rate', 'Acc'};

% --- Step 3: Perform Parametric Grey-Box Model Identification ---

% Define Grey-Box Model and a randomized initial guess
theta0 = [ -0.1068; 0.1192; -5.9755; -2.6478; -10.164; 450.0 ] .* (1 + 0.2 * (rand(6,1)-0.5));
sys_grey = idgrey('uav_longitudinal_grey', theta0, 'c');

% Configure and Run the 'greyest' Estimator
opt = greyestOptions;
opt.EnforceStability = false;
opt.Display = 'off'; % Keep this 'off' to prevent flooding the command window
opt.SearchMethod = 'auto';

% Perform the estimation
estimated_model = greyest(data, sys_grey, opt);

end % End of function