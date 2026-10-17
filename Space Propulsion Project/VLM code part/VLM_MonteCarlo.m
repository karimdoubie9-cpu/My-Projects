% Nominal values
clear
clc
close all
param.Mm = 18.01528; % Molecular mass, water
param.R = 8314.46/param.Mm; % Gas constant
param.g0 = 9.81;
param.rho = 1000; % Density of liquid water
param.gamma = 1.32; % Specifc heat ratio, water vapor

% Geometry
param.Dinj = 28e-6; % Injector orifice diameter, NOMINAL
param.Dt = 141e-6 ; % Throat diameter, NOMINAL
param.De = 373e-6; % Throat diameter, NOMINAL
param.Ainj = 0.25*pi*param.Dinj^2;
param.At = 0.25*pi*param.Dt^2;
param.Ae = 0.25*pi*param.De^2;
param.alpha = deg2rad(20); % Divergence half-angle

% Heating
param.wtot = 11.26; % Power [W], NOMINAL

% Injection
param.Ptank = 234.3e3; % Pressure of liquid water in tank
param.Ttank = 273+60; % Temperature of liquid water in tank
param.Cd = 0.7; % Discharge Coefficient of the injector
param.plenumVel = 0.36; % Velocity in feedline, for dynamic pressure calc.

disp("Nominal Results:")
VLM_simulator(param)
disp(VLM_simulator(param).thrust*1000)
%%%%%%%%%%%%%%%%% Stochastic simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input :
N = 30; % Size of sampling vectors. Num of iterations will be N^3
inc_Dt = 0.01*param.Dt;
inc_De = 0.050*param.De;
inc_Dinj = 0.05*param.Dinj;
V_inc=0.001; % V
i_inc=0.001; % A
inc_power = sqrt((param.wtot/5)^2*V_inc^2+5^2*i_inc^2) ;


% Setup plots
figure;
subplot(2,1,1);
h_mean = plot(NaN, NaN, 'b-');
ylabel('Mean');
xlabel('Iteration');

subplot(2,1,2);
h_std = plot(NaN, NaN, 'r-');
ylabel('Std Dev');
xlabel('Iteration');
update_after = 200;

% I set a different seed for each process with rng
rng(5);     
Dt_vec = param.Dt + inc_Dt * randn(N, 1);
rng(3);     
De_vec = param.De + inc_De * randn(N, 1);
rng(2);
Dinj_vec = param.Dinj + inc_Dinj * randn(N, 1);
rng(1);
wtot_vec = param.wtot + inc_power * randn(N,1);
rndInd = randperm(N^4) ;

% Create a matrix of all the possible combinations and randomize position
Comb = - ones(N^4,4) ;
ind = 0 ;
for i = 1 : N
    for j = 1 : N
        for k = 1 : N
            for l = 1 : N
                ind = ind + 1 ;
                Comb(rndInd(ind),:) = [Dt_vec(i) De_vec(j) Dinj_vec(k) wtot_vec(l)] ;
            end
        end
    end
end

% Simulate burn
thrust_result = - ones(N^4,1) ;
means = zeros(floor((N^4)/1000),1) ;
stds = means ;
for i = 1 : N^4
    param.At = 0.25*pi*Comb(i,1)^2;
    param.Ae = 0.25*pi*Comb(i,2)^2;
    param.Ainj = 0.25*pi*Comb(i,3)^2;
    param.wtot = Comb(i,4);
    thrust_result(i) = VLM_simulator(param).Isp ;

    if mod(i, update_after) == 0
        idx = i / update_after ;
        means(idx) = mean(thrust_result(1:i)) ;
        stds(idx) = std(thrust_result(1:i)) ;

        % Update plot
        set(h_mean, 'XData', (1:idx)*update_after, 'YData', means(1:idx)) ;
        set(h_std, 'XData', (1:idx)*update_after, 'YData', stds(1:idx)) ;
        drawnow ;  % Force plot update
    end
end

mean(thrust_result) 
std(thrust_result)