% MonteCarlo simulation of a SRM to compute burning time uncertainty
pkg load parallel

% Input
ThreadsToUse = 10;
N = 75;                     % Number of samples for each variable
TotalSamples = N^3;      % Total number of Monte Carlo simulations

% Chamber diameter
dc = 0.500e-3;
Inc_dc = (10/3)*1e-6;
% Nozzle throat diameter
dtht = 0.03303e-3;
Inc_dt = (1/3)*1e-6;
% Vieille’s law "a" [mm/s/bar^n]
a = 1.7270;
Inc_a = 0.0184;
% Vieille’s law "n"
n = 0.3821;
Inc_n = 0.0027;
% Cstar
Cstar = 1513.1382;
Inc_Cs = 14.7999;

% Initialize random generator
rng("shuffle");

% Generate random samples for each parameter
a_vet    = a     + Inc_a  * randn(N, 1);
n_vet    = n     + Inc_n  * randn(N, 1);
dtht_vet = dtht  + Inc_dt * randn(N, 1);
dc_vet   = dc    + Inc_dc * randn(N, 1);
Cs_vet   = Cstar + Inc_Cs * randn(N, 1);

% Smart sampling using random index selection (not full cartesian)
idx = randi(N, TotalSamples, 5);  % Each row is a random combo of indices

% Prepare matrix of input combinations
Comb = [a_vet(idx(:,1)), n_vet(idx(:,2)), dtht_vet(idx(:,3)), dc_vet(idx(:,4)), Cs_vet(idx(:,5))];

% Parallel simulation using pararrayfun
output = pararrayfun(ThreadsToUse, ...
    @(i) EndBurnerSIMP(Comb(i,1), Comb(i,2), Comb(i,3), Comb(i,4), Comb(i,5)), ...
    1:TotalSamples, "UniformOutput", false);

% Combine the outputs
result = cell2mat(output);           % result is [TotalSamples * 4 x 1]
result = reshape(result, 4, [])';    % Reshape to [TotalSamples x 4]

% Now compute the statistics properly
means = mean(result, 1);             % Average across rows
stds = std(result, 0, 1);            % Standard deviation across rows

% Print results
fprintf('Burn Time Avg: %.2f Std: %.3f s\n', means(1), stds(1));
fprintf('Thrust Avg: %.2f Std: %.3f mN\n', means(2).*1e3, stds(2).*1e3);
fprintf('Total Impulse Avg: %.2f Std: %.3f mNs\n', means(3).*1e3, stds(3).*1e3);
fprintf('Specific Impulse Avg: %.2f Std: %.3f s\n', means(4), stds(4));