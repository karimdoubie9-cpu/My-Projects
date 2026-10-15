function y = stackFRFData(g)
%% STACKFRFDATA Stack complex FRF points into a vector of real and imaginary part.
% 
% input
%   g: [Kx2] vector of non-parametric frequency response data (complex)
% output
%   y: [4x(K)] matrix of FRF data stacked in column (real part, imag part)

% Author: 
%   Simone Panza <simone.panza@polimi.it>
% 
% changelog
%   v1
%       first version
g=g.'; % g is originally 2*k
K = size(g, 1);

y = zeros(4, K);
y(1,:) = real(g(:,1));
y(2,:) = real(g(:,2));
y(3,:) = imag(g(:,1));
y(4,:) = imag(g(:,2));

% y=zeros(2,K);
% y(1,:) = real(g(:,1));
% y(2,:) = imag(g(:,1));
end

