function ref = setpoint(p0, pf, T, opts)
%SETPOINT  Step reference from p0 to pf, held for the rest of the run.
%
%   ref = traj.setpoint([0;0;1], [1;0;1], 35)
%   ref = traj.setpoint(p0, pf, T, t_step=0.001, yaw=0)
%
%   Used for the tuning cases of Task 2.2/2.3.1 (lateral and vertical steps).
%   The step is expressed with only four samples -- the From-Workspace block
%   interpolates linearly between them, and a near-zero rise time gives the
%   controller a genuine step rather than a fast ramp.
%
%   Velocity and acceleration references are held at zero throughout: the step
%   is a SETPOINT command, so there is no feedforward to hand the controller.

arguments
    p0          (3,1) double
    pf          (3,1) double
    T           (1,1) double {mustBePositive}
    opts.t_step (1,1) double = 0.001   % [s] time at which the setpoint jumps
    opts.yaw    (1,1) double = 0       % [rad] constant heading command
end

t = [0; opts.t_step; T];
p = [p0, pf, pf];
v = zeros(3,3);
a = zeros(3,3);
y = opts.yaw*ones(1,3);

ref = traj.pack(t, p, v, a, y);
end
