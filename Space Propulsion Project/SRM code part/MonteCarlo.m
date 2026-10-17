% MonteCarlo simulation of a SRM to compute burning time uncertainty
clear ;
close all ;
clc ;

% Input :
% Number of samples per random variable
N = 50 ;  
% Number of random variables 
NRND = 5 ;
% Chamber diameter
dc = 0.500e-3 ;
Inc_dc = (10/3)*1e-6 ;
% Nozzle throat diameter
dtht = 0.03303e-3 ;
Inc_dt = (1/3)*1e-6 ;
% Vieille’s law "a" [mm/s/bar^n]
a = 1.7270 ;
Inc_a = 0.0184 ;
% Vieille’s law "n"
n = 0.3821 ;
Inc_n = 0.0027 ;
% Cstar
Cstar = 1513.1382 ;
Inc_Cs = 14.7999 ;
% plot update rate
update_after = 1000 ;

% Random Vector Creation 
% I set a different seed for each process
rng("shuffle");     
a_vet = a + Inc_a * randn(N, 1);
rng("shuffle");   
n_vet = n + Inc_n * randn(N, 1);
rng("shuffle");   
dtht_vet = dtht + Inc_dt * randn(N, 1);
rng("shuffle");   
dc_vet = dc + Inc_dc * randn(N, 1);
rng("shuffle");   
Cs_vet = Cstar + Inc_Cs * randn(N, 1);

% Generate random index
rndInd = randperm(N^NRND) ;

% Create a matrix of all the possible combinations with random position
Comb = - ones(N^NRND,NRND) ;
ind = 0 ;
for i = 1 : N
    for j = 1 : N
        for k = 1 : N
            for y = 1 : N
                for z = 1 : N
                    ind = ind + 1 ;
                    Comb(rndInd(ind), :) = [a_vet(i), n_vet(j), dtht_vet(k), dc_vet(y), Cs_vet(z)] ;
                end
            end
        end
    end
end

% Initialize figure
figure;
subplot(2,1,1);
h_mean = plot(NaN, NaN, 'b-');
ylabel('Mean');
xlabel('Iteration');

subplot(2,1,2);
h_std = plot(NaN, NaN, 'r-');
ylabel('Std Dev');
xlabel('Iteration');

% Simulate Engine Burn
tburn = - ones(N^3,1) ;
F = tburn ;
I = tburn ;
Isp = tburn ;
means = zeros(floor((N^3)/1000),1) ;
stds = means ;
for i = 1 : N^(NRND-1)
    [tburn(i), F(i), I(i), Isp(i)] = EndBurnerSIM(Comb(i,1), Comb(i,2), Comb(i,3), Comb(i,4), Comb(i,5)) ;

    if mod(i, update_after) == 0
        idx = i / update_after ;
        means(idx) = mean(F(1:i)) ;
        stds(idx) = std(F(1:i)) ;
        
        % Update plot
        set(h_mean, 'XData', (1:idx)*update_after, 'YData', means(1:idx)) ;
        set(h_std, 'XData', (1:idx)*update_after, 'YData', stds(1:idx)) ;
        drawnow ;  % Force plot update
    end
end

%% Print results
fprintf('Burn Time Avg: %.2f Std: %.3f s\n', mean(tburn), std(tburn));
fprintf('Thrust Avg: %.2f Std: %.3f mN\n', mean(F).*1e3, std(F).*1e3);
fprintf('Total Impulse Avg: %.2f Std: %.3f mNs\n', mean(I).*1e3, std(I).*1e3);
fprintf('Specific Impulse Avg: %.2f Std: %.3f s\n', mean(Isp), std(Isp));

%{
N = 50 ;
Burn Time Avg: 0.40 Std: 0.012 s
Thrust Avg: 5.01 Std: 0.205 mN
Total Impulse Avg: 2.02 Std: 0.038 mNs
Specific Impulse Avg: 198.47 Std: 2.143 s
%}
%{
N = 75 ;
Burn Time Avg: 0.40 Std: 0.011 s
Thrust Avg: 4.99 Std: 0.185 mN
Total Impulse Avg: 2.02 Std: 0.033 mNs
Specific Impulse Avg: 198.83 Std: 2.155 s
%}
%{
N = 100 ;
Burn Time Avg: 0.41 Std: 0.012 s
Thrust Avg: 5.01 Std: 0.209 mN
Total Impulse Avg: 2.03 Std: 0.038 mNs
Specific Impulse Avg: 198.91 Std: 2.101 s
%}
%{
N = 100 ;
Burn Time Avg: 0.41 Std: 0.012 s
Thrust Avg: 4.98 Std: 0.187 mN
Total Impulse Avg: 2.02 Std: 0.032 mNs
Specific Impulse Avg: 198.64 Std: 1.867 s
%}
%{
N = 100 ;
Burn Time Avg: 0.40 Std: 0.012 s
Thrust Avg: 5.02 Std: 0.190 mN
Total Impulse Avg: 2.02 Std: 0.032 mNs
Specific Impulse Avg: 198.70 Std: 1.921 s
%}


