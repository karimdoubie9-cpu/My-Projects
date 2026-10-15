function ref = slalom(T, opts)
%SLALOM  Climbing sinusoidal / slalom trajectory.
%
%   ref = traj.slalom(150)
%
%   ref = traj.slalom(150, ...
%       forward_speed = 1.5, ...
%       amplitude     = 3.0, ...
%       rate          = 0.5, ...
%       climb         = 0.2);
%
%   The trajectory is
%
%       x_d(t) = x0 + V*t
%       y_d(t) = y0 + A*sin(w*t)
%       z_d(t) = z0 + vz*t
%
%   where:
%       V   = constant forward speed
%       A   = lateral oscillation amplitude
%       w   = lateral angular frequency [rad/s]
%       vz  = constant climb rate
%
%   Unlike the helix, the lateral acceleration repeatedly changes sign:
%
%       a_y = -A*w^2*sin(w*t)
%
%   so the UAV must continuously reverse its lateral thrust direction.
%   For a coplanar multirotor this requires repeated roll/pitch manoeuvres,
%   making the trajectory useful for testing:
%
%       - position-loop tracking bandwidth,
%       - position/attitude cascade interaction,
%       - lateral control authority,
%       - payload-induced mass/inertia changes,
%       - adaptive augmentation under persistent excitation.
%
%   The heading is aligned with the horizontal velocity:
%
%       psi_d = atan2(v_y, v_x)
%
%   so yaw is exercised naturally as the UAV changes direction through
%   the slalom rather than being artificially held fixed.
%
%   DEFAULT TRAJECTORY
%   ------------------
%       V  = 1.5 m/s
%       A  = 3.0 m
%       w  = 0.5 rad/s
%       vz = 0.2 m/s
%
%   This gives
%
%       max lateral velocity     = A*w     = 1.50 m/s
%       max lateral acceleration = A*w^2   = 0.75 m/s^2
%       max total speed          ~= 2.13 m/s
%
%   which stays comfortably below the 4-5 m/s project limit while
%   producing substantially more lateral excitation than the old helix.
%
%   The initial phase is chosen as zero:
%
%       y(0)   = y0
%       a_y(0) = 0
%
%   so the trajectory does not begin with an instantaneous lateral
%   acceleration demand. It does, however, begin with non-zero lateral
%   velocity. harness.input initializes the real vehicle consistently
%   from the first reference velocity, so this does not create an
%   artificial initial-condition transient.
%
%   See also traj.helix, traj.setpoint, traj.pack.

arguments
    T                  (1,1) double {mustBePositive}
    opts.forward_speed (1,1) double {mustBePositive}    = 1.5   % [m/s]
    opts.amplitude     (1,1) double {mustBeNonnegative} = 3.0   % [m]
    opts.rate          (1,1) double {mustBePositive}    = 0.5   % [rad/s]
    opts.climb         (1,1) double                     = 0.2   % [m/s]
    opts.x0            (1,1) double                     = 0.0   % [m]
    opts.y0            (1,1) double                     = 0.0   % [m]
    opts.z0            (1,1) double                     = 1.0   % [m]
    opts.dt            (1,1) double {mustBePositive}    = 0.01  % [s]
end

%% Time vector
t = 0:opts.dt:T;

V  = opts.forward_speed;
A  = opts.amplitude;
w  = opts.rate;
vz = opts.climb;

%% Position reference
%
% Forward progression in x, sinusoidal lateral motion in y,
% constant climb in z.
p = [ ...
    opts.x0 + V*t;
    opts.y0 + A*sin(w*t);
    opts.z0 + vz*t
    ];

%% Velocity reference
%
% Analytic derivative of p(t). Do not obtain this numerically: the
% controller should receive a smooth and internally consistent reference.
v = [ ...
    V*ones(size(t));
    A*w*cos(w*t);
    vz*ones(size(t))
    ];

%% Acceleration reference
%
% This is the important excitation:
%
%   a_y = -A*w^2*sin(w*t)
%
% Its sign reverses every half-period, requiring alternating lateral
% force and therefore alternating vehicle tilt.
a = [ ...
    zeros(size(t));
    -A*w^2*sin(w*t);
    zeros(size(t))
    ];

%% Desired yaw
%
% Point the vehicle along the instantaneous horizontal direction of travel.
% Since V > 0, the heading remains well defined throughout the trajectory.
yaw = unwrap(atan2(v(2,:), v(1,:)));

%% Sanity check against the project speed limit
speed = vecnorm(v, 2, 1);
vmax  = max(speed);

if vmax > 5
    warning('traj:slalom:speedLimit', ...
        ['Slalom reference reaches %.2f m/s, above the 5 m/s project ' ...
         'limit. Reduce forward_speed, amplitude or rate.'], vmax);
end

%% Pack into the format consumed by Final.slx
ref = traj.pack(t(:), p, v, a, yaw);

end