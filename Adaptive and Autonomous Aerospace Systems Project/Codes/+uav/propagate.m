function S = propagate(UAV, T, tspan, opts)
%PROPAGATE  Open-loop 6-DOF integration with a FIXED rotor command.
%
%   S = uav.propagate(UAV, T, [0 20])
%
%   Integrates the same equations the Simulink plant implements, driven by the
%   constant rotor thrusts in the trim struct T. With no feedback whatsoever,
%   the state can only stay put if the commanded wrench genuinely balances the
%   model -- which is exactly the property Task 1.2/1.3 asks to demonstrate.
%   Any residual shows up as a drift that grows with time, so a flat trace at
%   the 1e-12 level is a real verification rather than a controller hiding an
%   error.
%
%   State: x = [p_i(3); q(4); nu(6)],  nu = [v_b; omega_b]
%
%       p_i_dot = R(q) * v_b
%       q_dot   = 0.5 * W(q) * omega
%       nu_dot  = M \ ( w_g + w_aero + w_ctrl - cross_six(nu, M*nu) )

arguments
    UAV       struct
    T         struct
    tspan     (1,2) double
    opts.p0   (3,1) double = [0;0;0]
    opts.reltol (1,1) double = 1e-12
    opts.abstol (1,1) double = 1e-14
end

g   = 9.81;
e_3 = [0;0;1];

q0 = rotm2quat(T.R).';
x0 = [opts.p0; q0; T.nu];

w_ctrl = UAV.F1 * (UAV.F2 * (T.thrust_sign * T.T_r));   % constant by definition

odeopt = odeset('RelTol', opts.reltol, 'AbsTol', opts.abstol);
[t, X] = ode45(@(~,x) deriv(x), tspan, x0, odeopt);

S.t     = t;
S.p     = X(:,1:3);
S.q     = X(:,4:7);
S.v_b   = X(:,8:10);
S.omega = X(:,11:13);

% Inertial velocity, for the cruise (Task 1.3) check
S.v_i = zeros(size(S.v_b));
for k = 1:numel(t)
    S.v_i(k,:) = (quat2rotm(S.q(k,:)) * S.v_b(k,:).').';
end

    function dx = deriv(x)
        q  = x(4:7)/norm(x(4:7));
        nu = x(8:13);
        R  = quat2rotm(q.');
        v  = nu(1:3);
        om = nu(4:6);

        w_g = g*[UAV.m*(R.'*e_3); UAV.S*(R.'*e_3)];
        w_a = -UAV.D*nu;
        w_c = cross_six(nu, UAV.M*nu);

        nud = UAV.M \ (w_g + w_a + w_ctrl - w_c);

        W  = [-q(2:4).'; skew(q(2:4)) + q(1)*eye(3)];
        dx = [R*v; 0.5*W*om; nud];
    end
end

% ---------------------------------------------------------------------
function c = cross_six(a, b)
c = [skew(a(4:6)), zeros(3); skew(a(1:3)), skew(a(4:6))] * b;
end

function S = skew(x)
S = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
end
