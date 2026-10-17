% MonteCarlo simulation of a SRM to compute burning time uncertainty
pkg load parallel

% Input :
ThreadsToUse = 10;
N = 75 ;
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

% Prepare Inputs for parallel processing
% Extract 3 column vectors from Comb and convert them to cell arrays
in1 = num2cell(Comb(:,1) ) ;
in2 = num2cell(Comb(:,2) ) ;
in3 = num2cell(Comb(:,3) ) ;
in4 = num2cell(Comb(:,4) ) ;
in5 = num2cell(Comb(:,5) ) ;

% Simulate Engine Burn
output = parcellfun(ThreadsToUse, @(in1,in2,in3,in4,in5) EndBurnerSIMP(in1,in2,in3,in4,in5)...
    , in1, in2, in3, in4, in5, "UniformOutput", false);

% Convert the output cell array to a numeric array                                                     
result = cell2mat(output) ;

% Compute final Mean and STD
means = mean(result);
stds = std(result);

% Print results
fprintf('Burn Time Avg: %.2f Std: %.3f s\n', means(1), stds(1));
fprintf('Thrust Avg: %.2f Std: %.3f mN\n', means(2).*1e3, stds(2).*1e3);
fprintf('Total Impulse Avg: %.2f Std: %.3f mNs\n', means(3).*1e3, stds(3).*1e3);
fprintf('Specific Impulse Avg: %.2f Std: %.3f s\n', means(4), stds(4));