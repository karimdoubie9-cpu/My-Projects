function ref = multisine3d(T, opts)
%MULTISINE3D  Gentle 3-axis trajectory for PBMRAC validation.
%
%   ref = traj.multisine3d(180)
%
%   The default reference is
%
%       x_d = A_x sin(w_x t + phi_x)
%       y_d = A_y sin(w_y t + phi_y)
%       z_d = z_0 + A_z sin(w_z t + phi_z)
%       psi_d = constant
%
%   with different frequencies on the three axes.  The purpose is to vary the
%   commanded translational force persistently, while keeping the required
%   attitude excursions small.  This makes it useful for checking the PBMRAC
%   translational uncertainty model without deliberately exciting large
%   translation/rotation coupling.
%
%   Default values:
%       A   = [2.00; 1.50; 0.75] m
%       w   = [0.25; 0.31; 0.18] rad/s
%       z0  = 2 m
%       phase = [0; pi/3; 0] rad
%       yaw   = 0 rad
%
%   The default maximum component accelerations are approximately
%       [0.125, 0.144, 0.024] m/s^2.
%
%   Output follows the standard traj.pack convention:
%       [p(3), v(3), a(3), yaw].
%
%   See also traj.pack, traj.helix.

arguments
    T                 (1,1) double {mustBePositive}
    opts.amplitude    (3,1) double = [2.00; 1.50; 0.75]   % [m]
    opts.rate         (3,1) double = [0.25; 0.31; 0.18]   % [rad/s]
    opts.phase        (3,1) double = [0; pi/3; 0]            % [rad]
    opts.center       (3,1) double = [0; 0; 2]            % [m]
    opts.yaw          (1,1) double = 0                    % [rad]
    opts.dt           (1,1) double {mustBePositive} = 0.01
end

t = 0:opts.dt:T;

A   = opts.amplitude;
w   = opts.rate;
phi = opts.phase;
c   = opts.center;

arg = w*t + phi;

% Analytic kinematics: use exact derivatives rather than finite differences.
p = c + A .* sin(arg);
v =     A .* w .* cos(arg);
a =    -A .* (w.^2) .* sin(arg);

% Keep heading fixed.  The validation is aimed at the translational adaptive
% model, so unnecessary yaw excitation is intentionally avoided.
yaw = opts.yaw * ones(1, numel(t));

% Sanity check against the project trajectory-speed requirement.
speed = vecnorm(v, 2, 1);
if max(speed) > 5
    warning('traj:multisine3d:speedLimit', ...
        'Reference reaches %.2f m/s; reduce amplitude/rate.', max(speed));
end

ref = traj.pack(t(:), p, v, a, yaw);
end
