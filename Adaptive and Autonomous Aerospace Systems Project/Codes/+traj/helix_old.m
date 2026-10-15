function ref = helix_old(T, opts)
%HELIX_OLD  Original Task-2/Task-3 helix, with no approach manoeuvre.
%
%   This is a frozen copy of the original traj.helix used before the
%   rest-to-helix approach was introduced.  It is kept under a separate name
%   so the Monte Carlo can be run on the old trajectory without changing the
%   current +traj/helix.m.
%
%   Default speed: sqrt((R*w)^2 + vz^2) = 0.71 m/s.

arguments
    T           (1,1) double {mustBePositive}
    opts.radius (1,1) double = 5      % [m]
    opts.rate   (1,1) double = 0.1    % [rad/s]
    opts.climb  (1,1) double = 0.5    % [m/s]
    opts.z0     (1,1) double = 0      % [m]
    opts.dt     (1,1) double = 0.01   % [s]
end

t   = (0:opts.dt:T);
R   = opts.radius;
w   = opts.rate;
phi = w*t;

p = [ R*cos(phi);
      R*sin(phi);
      opts.z0 + opts.climb*t ];

v = [-R*w*sin(phi);
      R*w*cos(phi);
      opts.climb*ones(size(t)) ];

a = [-R*w^2*cos(phi);
     -R*w^2*sin(phi);
      zeros(size(t)) ];

yaw = unwrap(atan2(v(2,:), v(1,:)));

ref = traj.pack(t(:), p, v, a, yaw);
end
