function dydtheta = computeFRFsensitivity( omega_vec, G, theta_0 )

%% COMPUTEFRFSENSITIVITY Compute sensitivity of the model FRF with respect to parameters.
% 
% input
%   omega_vec   frequency grid
%   G           model
%   theta_0     parameter value
% output
%   dydtheta    sensitivity of FRF (real and imag parts) wrt variations of
%               theta

% Author: 
%   Simone Panza <simone.panza@polimi.it>
% 
% changelog
%   v1
%       first version

% perturb theta_0 by a small fraction of it 
% PERT = 1e-5;     % perturb percentage
% 
% % K = length(omega_vec);
 nth = length(theta_0);
% 
 % evaluate FRF at theta_0
 Gth = G(theta_0);
 g = squeeze(freqresp(Gth, omega_vec));
% 
 Y0 = stackFRFData(g); % 4*k 
% 
 dydtheta = struct([]);
% 
% for ii = 1:nth
% 
%     % perturb theta
%     dtheta = PERT*abs(theta_0(ii));
%     pert_th = zeros(size(theta_0));
%     pert_th(ii) = dtheta;
%     theta = theta_0 + pert_th;   % perturbed parameter
% 
%     % evaluate FRF at theta
%     Gth = G(theta);
%     g = squeeze(freqresp(Gth, omega_vec));
%     Y = stackFRFData(g);
% 
%     % compute sensitivity
%     dydtheta(ii).dy = (Y - Y0)/dtheta;
% 
% end
PERT = 1e-3;          % start bigger than 1e-5
ABS_MIN = 1e-4;       % minimum absolute perturbation

for ii = 1:nth
    dtheta = max(PERT*abs(theta_0(ii)), ABS_MIN);

    theta_p = theta_0; 
    theta_p(ii) = theta_p(ii) + dtheta;
    theta_m = theta_0; 
    theta_m(ii) = theta_m(ii) - dtheta;

    gp = squeeze(freqresp(G(theta_p), omega_vec));
    gm = squeeze(freqresp(G(theta_m), omega_vec));

    Yp = stackFRFData(gp);
    Ym = stackFRFData(gm);

    dydtheta(ii).dy = (Yp - Ym) / (2*dtheta);   % central diff
end

end

