%% SCRIPT: Task 3.1 - Monte Carlo Simulation for a Specific Noise Condition
    % This script runs the Monte Carlo simulation (Task 3.1) specifically for
    % the "best" noise configuration identified in Task 2 ('Perturbed 1').
    % It modifies the simulation parameters ONCE before the loop begins.
    %N.B: in step1.2 the seed of all white noise blocks in the simulink model is changed to sim_seed
            %  whose value at each iteration will be changed automotically  
clear;
clc;
close all;
addpath(genpath(fullfile(pwd, 'simulator_antx_identified-single_axis')));

%% ========================================================================
%  STAGE 1: SETUP AND SIMULATION PARAMETER MODIFICATION
%  ========================================================================

% --- Step 1.1: Define Simulation Parameters ---
    numSimulations = 1000; % Justification: Sufficient for stable statistics
    % Define the desired noise condition from Task 2
    % (Perturbed 1: Low Accel Noise)
    target_acc_std_dev = 0.35;      % [m/s^2]
    target_ang_rate_std_dev = deg2rad(1.3); % [rad/s]

% --- Step 1.2: Modify the 'main_quad_ANTX.m' script ---
    % This is done once before the loop for efficiency.
    fprintf('Modifying simulation script for the target noise condition...\n');
    main_script_content = fileread('main_quad_ANTX.m');
    % Replace the noise parameter lines using regular expressions
    modified_script = regexprep(main_script_content, ...
        'noise.acc_stand_dev\s*=\s*noise.Enabler\s*\*\s*[\d\.]*;\s*%', ...
        sprintf('noise.acc_stand_dev = noise.Enabler * %.4f; %%', target_acc_std_dev));
    modified_script = regexprep(modified_script, ...
        'noise.ang_rate_stand_dev\s*=\s*noise.Enabler\s*\*\s*deg2rad\([\d\.]*\);\s*%', ...
        sprintf('noise.ang_rate_stand_dev = noise.Enabler * %.4f; %%', target_ang_rate_std_dev));
    % Write the temporary script that will be called inside the loop
    fid = fopen('main_quad_ANTX_MC.m', 'w');
    fprintf(fid, '%s', modified_script);
    fclose(fid);
    fprintf('Temporary simulation script ''main_quad_ANTX_MC.m'' created.\n\n');


    % --- Step 1.2: Automically configure simulink seeds and set them 
                    % to sim_seed which will changed in the main for loop ---
        % This section finds all Random Number blocks and sets their Seed parameter 
        % to the variable name 'sim_seed'. Run this ONCE.
        model_name = 'Simulator_Single_Axis';
        load_system(model_name);
        % Find all blocks of type 'RandomNumber' (standard White Noise blocks)
        noise_blocks = find_system(model_name, 'LookUnderMasks', 'on', 'BlockType', 'RandomNumber');
        % Loop through found blocks and set the Seed parameter
        for k = 1:length(noise_blocks)
            % Set the Seed parameter to an expression that adds sim_seed and the block index k
            % Use MATLAB-expression syntax so Simulink evaluates it at runtime
            set_param(noise_blocks{k}, 'Seed', sprintf('sim_seed + %d', k));
            fprintf('Updated Seed for block: %s\n', noise_blocks{k});
        end
        % Save and close the system so changes persist for the simulation run
        save_system(model_name);
        close_system(model_name);
%% ========================================================================
%  STAGE 2: MONTE CARLO SIMULATION LOOP
%  ========================================================================
    % Pre-allocate storage for results
    results = struct('output', cell(numSimulations, 1));
    theta_hat_mc = zeros(6, numSimulations);
    sigma_theta_mc = zeros(6, numSimulations);
    
    fprintf('====== STARTING MONTE CARLO SIMULATION (%d RUNS) ======\n', numSimulations);
    prog_bar = waitbar(0, 'Starting Monte Carlo Simulation...');
    
    for i = 1:numSimulations
        % Set the seed for the Simulink noise blocks.
        sim_seed = i;
        % Run the MODIFIED main_quad_ANTX_MC.m. This ensures each simulation run uses
        % the target noise standard deviations.
        run('main_quad_ANTX_MC.m'); 
        % Call the identification function to process the simulation output
        results(i).output = run_model(out);
        % Store the key results for analysis
        theta_hat_mc(:,i) = results(i).output.ParameterVector;
        P_theta = getcov(results(i).output);
        sigma_theta_mc(:,i) = sqrt(diag(P_theta));
        % Update progress
        waitbar(i/numSimulations, prog_bar, sprintf('Running Identification %d of %d...', i, numSimulations));
    end
    close(prog_bar);
    delete('main_quad_ANTX_MC.m'); % Clean up
    fprintf('\n====== MONTE CARLO SIMULATION COMPLETE ======\n');
    %% Save MC results in MC_results_Perturbed1.mat
    MC_results.numSimulations = numSimulations;
    MC_results.theta_hat_mc   = theta_hat_mc;
    MC_results.sigma_theta_mc = sigma_theta_mc;
    MC_results.results_struct = results;
    MC_results.target_acc_std_dev      = target_acc_std_dev;
    MC_results.target_ang_rate_std_dev = target_ang_rate_std_dev;
    save('MC_results_Perturbed1.mat', 'MC_results');
    fprintf('Monte Carlo results saved to MC_results_Perturbed1.mat\n');

% Remove the folder and all its subfolders from MATLAB path
folderToRemove = fullfile('simulator_antx_identified-single_axis');
rmpath(genpath(folderToRemove));