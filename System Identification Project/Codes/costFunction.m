function [Rinv_out, J, e ] = costFunction( g_m, omega_vec, G, theta,Rinv)
%% COST FUNCTION Evaluate cost function and error
% 
% input
%   g_m         measured FRF data
%   omega_vec   frequency grid [rad/s]
%   G           tf model parametrized in theta
%   theta       parameter vector
%   R measurement covariance matrix
% output
%   J           cost function value
%   e           estimation error as a function of frequency

% Author: 
%   Simone Panza <simone.panza@polimi.it>
% 
% changelog
%   v1
%       first version

% Stack complex FRF points into a vector of real and imaginary parts
Y_m = stackFRFData(g_m);

K = length(omega_vec);
nth = length(theta);

Gth = G(theta);                          % estimaned TF
g = squeeze(freqresp(Gth, omega_vec));   % estimaned FRF

% estimated output
Y = stackFRFData(g);

% error
e = Y_m - Y;

% cost function computation
J = 0;
R=zeros(4,4);
e_m=g_m-g;
for kk = 1:K
R = R + (1/(K-nth)) * abs((e(:,kk) * e(:,kk).'));
end
%Rdiag=diag([R(1,1) R(2,2)]);
%R_new = [Rdiag/2 zeros(2,2); zeros(2,2) Rdiag/2];
R_new = diag(diag(R));
%R_new=eye(4);
%R_new=diag([R/2 R/2]);
Rinv=inv(R_new);

for kk = 1:K
    J = J + 1/2*(e(:,kk).')*eye(4)*e(:, kk);
end
Rinv_out=Rinv;
end
