function [theta_hat,Cov_theta, J_hat] = outputError_FD2( g_m, omega_vec, G_th, theta_0,R)

%% OUTPUTERROR_FD output error method in frequency domain

% Input
%   g_m         measured FRF data
%   omega_vec   frequency grid [rad/s]
%   G_th        tf model parametrized in theta
%   theta_0     parameter vector initial guess
% 
% Output
%   theta_hat   optimal parameter vector
%   J_hat       optimal cost function value

% Author: 
%   Simone Panza <simone.panza@polimi.it>
% 
% changelog
%   v1
%       first version

% parameters for convergence criteria
TRESHOLD_THETA = 1e-3;   % absolute treshold on the parameter vector
TRESHOLD_J = 1e-3;       % absolute threshold on cost function
Nmax = 50;               % max number of iterations

assert(iscolumn(theta_0), 'Parameter vector must be a column vector.');

K = length(omega_vec);
nth = length(theta_0);

% initialization
FLAG_CONVERGENCE = false;  % do at least one iteration
theta = theta_0;           % init condition
c = 0;                     % iteration counter

% initial cost function
R_new = [R/2 zeros(2,2); zeros(2,2) R/2];
Rinv=inv(R_new);
[Rinv,J, e] = costFunction(g_m, omega_vec, G_th, theta,Rinv);


% covariance matrix --> TODO: estimate it with sample variance estimator
%Rinv = eye(2);

% debug
%fprintf('Initial guess\t cost function value %f ',J);
% for ii = 1:nth
%     fprintf('\ttheta%d %f ', ii, theta(ii));
% end
% fprintf('\n');

% NEWTON RAPSON METHOD
while ~FLAG_CONVERGENCE
    
    % compute gradient and hessian of J wrt theta

    % compute sensitivity of the model FRF with respect to parameters 
    s = computeFRFsensitivity(omega_vec, G_th, theta);
    
    % compute jacobian and hessian
    G = zeros(nth, 1);
    H = zeros(nth, nth);
    
    for kk = 1:K
        dydtheta = zeros(4, nth);

        % compute sensitivity of the model FRF with respect to parameters
        for ii = 1:nth
            dydtheta(:, ii) = s(ii).dy(:, kk);
        end

        t = -e(:,kk)' * Rinv * dydtheta;

        G = G + t';  % jacobian
        
        H = H + dydtheta' * Rinv * dydtheta;  % hessian
    end
    
    % compute the increment (Newton-Raphson)
    delta_theta = -H\G;
    
    % update the estimate
    theta_new = theta + delta_theta;
    
    % compute cost function
    [Rinv,J_new, e_new] = costFunction(g_m, omega_vec, G_th, theta_new,Rinv);
    %[Rinv,J, e] = costFunction(g_m, omega_vec, G_th, theta,Rinv);

    % check convergence
    if (norm(delta_theta) < TRESHOLD_THETA ...
            || norm(J_new - J) < TRESHOLD_J ...
            || c > Nmax)
        FLAG_CONVERGENCE = true;
    end
    
    % update
    theta = theta_new;
    c = c + 1;
    J = J_new;
    e = e_new;
    
    % DEBUG
    %fprintf('Iteration #%d\t cost function value %f ',c,J);
    % for ii=1:nth
    %     fprintf('\ttheta%d %f ', ii, theta(ii));
    % end
    % fprintf('\n');

end

% best estimate
theta_hat = theta;
J_hat = J;

dof = K * 4 - nth; 
sigma_sq = (2 * J) / dof; 
Cov_theta = sigma_sq * inv(H);
end

