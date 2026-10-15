function ad = adaptive(UAV, par_g, opts)
%ADAPTIVE  Predictor-based MRAC augmentation parameters (Task 3.2 / 3.3).
%
%   ad = design.adaptive(UAV, g)
%
%   Uncertain translational dynamics assumed for the augmentation:
%
%       p_dot = v
%       v_dot = -g*e3 + (1/m)*( Lambda_f*u + f_e )
%
%   with Lambda_f a diagonal control-effectiveness matrix and f_e a constant
%   exogenous force.  The predictor implemented in 'adaptive controller/
%   Predictor dynamics' is
%
%       p_hat_dot = v_hat + L_p*(p - p_hat)
%       v_hat_dot = -g*e3 + (1/m)*( Lambda_hat*u + theta_hat ) + L_v*(v - v_hat)
%
%   where u is the ACHIEVABLE force u_h = u*min(1, T_MAX/||u||) (saturation
%   hedging, block 'achievable u'). The update laws ('theta_hat', 'lambda hat
%   dot', 'Lambda_dynamic') are, with s = (v - v_hat)/m,
%
%       robust_on = 0:  theta_hat_dot  = Gamma_theta * s
%                       Lambda_hat_dot = Gamma_lambda .* u_h .* s
%       robust_on = 1:  theta_hat_dot  = Gamma_theta * (s - sigma_e*||s||*theta_hat)
%                       Lambda_hat_dot = Proj(Lambda_hat, Gamma_lambda .* u_h .* s)
%
%   Proj is the smooth projection onto lambda in [0.35, 1.10] defined by
%   lambda_proj_center / r0 / r1 (MATLAB Function in 'Lambda_dynamic').
%   sigma_lambda is NOT used by the live model (its leak term is commented out).
%
%   NOTE ON INITIAL CONDITIONS: theta_hat(0) = [0;0;0] and Lambda_hat(0) =
%   [1;1;1] are hard-coded in the integrators. The predictor state is reset to
%   the MEASURED [p; v] every time the adaptive subsystem is enabled.

arguments
    UAV   struct
    par_g (1,1) double
   opts.L_p = 150;
opts.L_v = 150;
opts.sigma_e (1,1) double = 0.02

% One gain set for Task 3, Task 4 and the Monte Carlo (the Monte Carlo
% values). The gamma_lambda z entry was raised x10 (0.15 -> 1.5): with the
% old value lambda_hat_z converged too slowly and theta_hat absorbed the
% package weight instead; x10 cuts the Monte Carlo median RMS by 43 % and
% brings lambda_hat_z to the mass/effectiveness value.
% Previous (Task 3 only): Gamma_theta diag(144,144,115.2), Gamma_lambda [0.9;0.9;0.18].
opts.Gamma_theta  = diag([144 144 57.6]);
opts.Gamma_lambda = [0.75; 0.75; 1.5];

opts.sigma_lambda = 0.01; % was 0.01

opts.adapt_mass = UAV_nominal.m;   % 1.2 kg                 % controller's believed mass
end

% The controller's belief about the mass. It deliberately stays at the NOMINAL
% vehicle mass: the extra package mass is precisely the uncertainty the
% adaptation is supposed to absorb (Task 3.1 says the package is unknown).
if isnan(opts.adapt_mass)
    ad.m = UAV.m;
else
    ad.m = opts.adapt_mass;
end

ad.g            = par_g;
ad.L_p          = opts.L_p;
ad.L_v          = opts.L_v;
ad.Gamma_theta  = opts.Gamma_theta * eye(3);
ad.Gamma_lambda = opts.Gamma_lambda;
ad.sigma_lambda = opts.sigma_lambda;
ad.sigma_e = opts.sigma_e;
ad.lambda_proj_center = 0.725;
ad.lambda_proj_r0     = 0.275;
ad.lambda_proj_r1     = 0.375;
end
