function pos = position(UAV, opts)
%POSITION  Outer-loop position gains from an LQR on the per-axis design model.
%
%   pos = design.position(UAV)
%   pos = design.position(UAV, max_pos_err=..., max_vel=..., ...)
%
%   Design model (one translational axis, thrust treated as a direct force --
%   valid because the attitude loop is faster: 5.16 vs 1.66 rad/s, ~3x):
%
%       xdot = [0 1 0; 0 0 0; 1 0 0] x + [0; 1/m; 0] f,   x = [p; v; INT(p-p_d)]
%
%   The integral state removes the steady-state offset caused by an unmodelled
%   constant force (added package weight, steady wind), which is what the
%   baseline needs before the adaptive augmentation is switched on.
%
%   Q and R come from Bryson's rule, i.e. each weight is the reciprocal of the
%   square of the largest acceptable excursion of that quantity.  Stating the
%   design that way makes the tuning defensible in the report rather than a set
%   of magic numbers.
%
%   Returns K_p, K_v, K_i (the ctrl.* gains the model reads) plus the design
%   model matrices, which run_task2_frequency.m reuses.

arguments
    UAV struct
    opts.max_pos_err (1,1) double = 0.5     % [m]    tolerable position error
    opts.max_vel     (1,1) double = 4     % [m/s]  speed cap (exam: 4-5 m/s)
    opts.max_int_err (1,1) double = 2     % [m s]  tolerable integrated error
    opts.max_force   (1,1) double = 4    % [N]    control force budget per axis
end

% A = [0 1 0;
%      0 0 0;
%      1 0 0];
% B = [0; 1/UAV.m; 0];
% 
% Q = diag([1/opts.max_pos_err^2, 1/opts.max_vel^2, 1/opts.max_int_err^2]);
% R = 1/opts.max_force^2;
% 
% K = lqr(A, B, Q, R);
% 
% pos.K_p = K(1);
% pos.K_v = K(2);
% pos.K_i = K(3);
% control system deisgn is used instead.
scale=1;
pos.K_p = 8.638*scale;
pos.K_v = 5.2135*scale;
pos.K_i = 1.5747*scale;

pos.K_p_Z = 10.274;
pos.K_v_Z = 4.4458;
pos.K_i_Z = 1.4;


%pos.A = A;  pos.B = B;  pos.Q = Q;  pos.R = R;  pos.K = K;
pos.limits = opts;
end
