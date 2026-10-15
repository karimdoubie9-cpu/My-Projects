
function [A,B,C,D] = uav_longitudinal_grey(theta, Ts)
% Defines the continuous-time state-space model for the UAV's longitudinal
% dynamics, parameterized by the stability and control derivatives.
% This function is called by the 'idgrey' object during the estimation process.
%
% Parameter Vector Definition:
% theta = [Xu, Xq, Mu, Mq, X_delta, M_delta]'
%
% State Vector:       x = [u, q, theta]'
%   u:     Longitudinal velocity (body frame) [m/s]
%   q:     Pitch rate [rad/s]
%   theta: Pitch angle [rad]
%
% Input Vector:       u = [delta_lon]'
%   delta_lon: Total normalized pitching moment [-]
%
% Output Vector:      y = [q, ax]'
%   q:     Pitch rate [rad/s]
%   ax:    Longitudinal acceleration (body frame) [m/s^2]
%
% Argument 'Ts' is the sampling time, which is ignored here because we are
% defining a continuous-time ('c') model.

% -------------------------------------------------------------------------
% 1. Unpack Parameters from the Input Vector 'theta'
% -------------------------------------------------------------------------
Xu      = theta(1);
Xq      = theta(2);
Mu      = theta(3);
Mq      = theta(4);
X_delta = theta(5);
M_delta = theta(6);

% Define physical constants.
g = 9.81; % Acceleration due to gravity [m/s^2]

% -------------------------------------------------------------------------
% 2. Define the Continuous-Time State-Space Matrices (A, B, C, D)
% -------------------------------------------------------------------------

% The A matrix describes the system's internal dynamics (dx/dt = Ax + Bu).
% Rows correspond to d(u)/dt, d(q)/dt, d(theta)/dt.
A = [ Xu,  Xq, -g;  % u_dot = Xu*u + Xq*q - g*theta + ...
      Mu,  Mq,  0;   % q_dot = Mu*u + Mq*q + ...
       0,   1,  0 ];  % theta_dot = q

% The B matrix describes how the input affects the states.
B = [ X_delta;       % Input effect on u_dot
      M_delta;       % Input effect on q_dot
      0 ];           % Input has no direct effect on theta_dot


C = [ 0,   1,   0;   % y1 = q
      Xu,  Xq, 0 ]; % 
D = [ 0;             % No direct feedthrough to output q
      X_delta ];     % Feedthrough from delta_lon to ax

end