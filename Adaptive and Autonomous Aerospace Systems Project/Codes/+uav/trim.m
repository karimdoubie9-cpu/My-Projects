function T = trim(UAV, opts)
%TRIM  Rotor commands that hold the UAV in equilibrium (Task 1.2 / 1.3).
%
%   T = uav.trim(UAV)                       % hover                  [1.2]
%   T = uav.trim(UAV, v_i=[3;0;0])          % constant inertial vel  [1.3]
%   T = uav.trim(UAV, v_i=[3;0;0], yaw=pi/4)
%
%   EQUILIBRIUM CONDITION
%   Steady flight means the body-frame generalized velocity is constant, so
%   nu_dot = 0 and the wrenches must cancel:
%
%       0 = w_gravity + w_aero + w_control - cross_six(nu, M*nu)
%
%   with, as implemented in the model,
%       w_gravity = g*[ m*R'*e_3 ; S*R'*e_3 ]
%       w_aero    = -D*nu
%       w_control = F1 * [T_tot; tau_x; tau_y; tau_z]
%
%   HOVER (1.2). R = I, nu = 0, so both the aerodynamic and gyroscopic terms
%   vanish and only weight has to be carried:  F1*u = -w_gravity.
%
%   CRUISE (1.3). At constant INERTIAL velocity v_i the body-frame velocity is
%   v_b = R'*v_i, which is constant only if the attitude is constant -- so
%   omega = 0 and the gyroscopic term drops, but the drag D_fv*v_b does not.
%   Level flight cannot cancel a horizontal drag force, because a coplanar rotor
%   set produces no side force: the vehicle must TILT so that a component of its
%   thrust opposes the drag. Attitude is therefore an unknown of the trim
%   problem, not an input, and the equations are solved for it.
%
%   The 4-vector u = [T_tot; tau] is mapped to individual rotor thrusts with the
%   pseudo-inverse and then to commanded speeds through T_r = k_f*Omega_r^2.
%
%   SIGN CONVENTION. Whether positive T_tot points along +z_body or -z_body
%   depends on the model's axis convention; rather than assume, the solved rotor
%   thrusts are checked and the collective sign flipped if they all come out
%   negative. T.thrust_sign records which convention was detected, and
%   T.residual is the norm of the wrench imbalance -- if that is not ~1e-12 the
%   trim is not a trim, and the caller should say so rather than plot it.

arguments
    UAV      struct
    opts.v_i (3,1) double = [0;0;0]   % desired constant inertial velocity [m/s]
    opts.yaw (1,1) double = 0         % desired heading [rad]
end

g   = 9.81;
e_3 = [0;0;1];

%% ---- Solve for the attitude that balances the horizontal drag -------
% Unknowns: roll and pitch (yaw is commanded). At hover the drag is zero and
% the solution collapses to level flight, so the same code covers both tasks.
solve = @(rp) residual_force(rp, UAV, opts, g, e_3);

rp0 = [0; 0];
if any(abs(opts.v_i(1:2)) > 0)
    % Small-angle first guess: tilt ~ drag / weight
    Fd  = UAV.D_fv * opts.v_i;
    rp0 = [ Fd(2)/(UAV.m*g); -Fd(1)/(UAV.m*g) ];
end

rp = fsolve_local(solve, rp0);

roll = rp(1); pitch = rp(2); yaw = opts.yaw;
R    = eul2rotm([yaw pitch roll], 'ZYX');

%% ---- Wrenches at the trimmed attitude -------------------------------
v_b   = R.' * opts.v_i;
nu    = [v_b; zeros(3,1)];
w_g   = g*[UAV.m*(R.'*e_3); UAV.S*(R.'*e_3)];
w_aer = -UAV.D*nu;
w_gyr = cross_six(nu, UAV.M*nu);

% Required control wrench, and the 4 components the rotors can actually produce
w_req = w_gyr - w_g - w_aer;
u     = UAV.F1 \ w_req;                  % [T_tot; tau_x; tau_y; tau_z]

%% ---- Rotor thrusts and commanded speeds -----------------------------
T_r = UAV.F2_pinv * u;

thrust_sign = 1;
if all(T_r < 0)
    % Opposite collective convention: thrust is reported along -z_body.
    thrust_sign = -1;
    T_r = -T_r;
    u   = -u;
end

Omega = sqrt(max(T_r, 0) / UAV.k_f);

%% ---- Verify ---------------------------------------------------------
w_ctrl   = UAV.F1 * (UAV.F2 * (thrust_sign*T_r));
residual = norm(w_g + w_aer + w_ctrl - w_gyr);

T = struct( ...
    'Omega',       Omega, ...            % [rad/s] commanded rotor speeds
    'Omega_rpm',   Omega*60/(2*pi), ...
    'T_r',         T_r, ...              % [N] individual rotor thrusts
    'u',           u, ...                % [T_tot; tau] collective + moments
    'R',           R, ...
    'euler',       [roll pitch yaw], ... % [rad]
    'v_i',         opts.v_i, ...
    'v_b',         v_b, ...
    'nu',          nu, ...
    'thrust_sign', thrust_sign, ...
    'residual',    residual, ...
    'feasible',    all(Omega >= UAV.Omega_min & Omega <= UAV.Omega_max));
end

% ---------------------------------------------------------------------
function r = residual_force(rp, UAV, opts, g, e_3)
%RESIDUAL_FORCE  Horizontal force imbalance for a candidate roll/pitch.
%   A coplanar rotor set produces force only along body z, so the two
%   horizontal components of the total force are what the attitude must zero.
R    = eul2rotm([opts.yaw rp(2) rp(1)], 'ZYX');
v_b  = R.'*opts.v_i;
nu   = [v_b; zeros(3,1)];
w_g  = g*[UAV.m*(R.'*e_3); UAV.S*(R.'*e_3)];
w_a  = -UAV.D*nu;
w_gy = cross_six(nu, UAV.M*nu);
w    = w_gy - w_g - w_a;                 % force the rotors must supply
r    = w(1:2);                           % ... but they can only supply w(3)
end

function c = cross_six(a, b)
%CROSS_SIX  The 6-D cross-product operator used by the plant model.
S = @(x)[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
c = [S(a(4:6)), zeros(3); S(a(1:3)), S(a(4:6))] * b;
end

function x = fsolve_local(f, x0)
%FSOLVE_LOCAL  Damped Newton on a 2-D residual, so no Optimization Toolbox.
x = x0;
for it = 1:100
    r = f(x);
    if norm(r) < 1e-14, break, end
    J = zeros(2,2);
    h = 1e-7;
    for j = 1:2
        dx = zeros(2,1); dx(j) = h;
        J(:,j) = (f(x+dx) - r)/h;
    end
    if rcond(J) < 1e-12, break, end
    x = x - J\r;
end
end
