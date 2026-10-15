%% ========================================================================
%  TASK 1.2 / 1.3 - TRIM OF THE HEXACOPTER
%  ------------------------------------------------------------------------
%     M*nu_dot + cross_six(nu, M*nu) = w_gravity + w_aero + w_control     (1)
%
%     w_gravity = g * [ m*R'*e_3 ; S*R'*e_3 ]    weight, resolved in body axes
%     w_aero    = -D * nu                        linear damping
%     w_control = F1 * u,   u = [T; tau_x; tau_y; tau_z]
%
%  F1 = [0_(2x4); I_4] says the six rotors, being COPLANAR, generate force
%  along the body z axis only. 
%  Setting nu_dot = 0 in (1) gives the trim condition
%
%     0 = w_gravity + w_aero + w_control - cross_six(nu, M*nu)            (2)
%
%     Level 1  ANALYTIC
%     Level 2  NUMERIC  
%     Level 3  DYNAMIC  - open loop

clear; close all; clc;
viz.defaults();   % Times New Roman, heavier lines (see +viz/defaults.m)

UAV = uav.params();

m    = UAV.m;             % 1.2 kg
g    = 9.81;              % m/s^2
k_f  = UAV.k_f;           % 3.65e-6 N/(rad/s)^2
Nr   = UAV.Nr;            % 6 rotors
d_xy = UAV.D_fv(1,1);     % 0.26 N/(m/s) horizontal drag coefficient

T_END    = 20;            % [s] open-loop propagation horizon
V_CRUISE = 2;             % [m/s] cruise speed along inertial x (within 4-5)

fprintf('Hexacopter: m = %.2f kg, weight = %.4f N, k_f = %.3g, %d rotors\n\n', ...
        m, m*g, k_f, Nr);


%% ========================================================================
%  TASK 1.2  -  HOVER TRIM   (LEVEL 1: ANALYTIC)
%  ========================================================================
%  Conditions: p = const  =>  v_i = 0, and we choose the level attitude R = I.
%  Then nu = 0, so BOTH the damping term (-D*nu) and the gyroscopic term
%  cross_six(nu, M*nu) vanish. The centre of mass is at the body origin
%  (r_bg = 0 => S = 0), so gravity also produces no moment.
%
%  The trim condition (2) collapses to one scalar equation:
%
%        w_gravity + F1*u = 0
%     => [0; 0; m*g; 0; 0; 0] + [0; 0; T; tau] = 0
%     => T   = -m*g       (thrust opposes weight; the sign says it acts along
%     => tau = 0           -z_body, i.e. "up" in this NED frame)
%
%  ALLOCATION. F2 maps rotor thrusts to the wrench and its first row is all
%  ones, so the collective is just the sum:  sum_r T_r = T.
%
%        T_r = T / N_r = -m*g/6      for every r
%
%  Finally the propeller characteristic T_r = k_f*Omega_r^2 inverts to
%
%        Omega_r = sqrt( |T_r| / k_f )

T_hover      = -m*g;                        % [N] collective thrust
Tr_hover     = T_hover / Nr;                % [N] per rotor (signed)
Omega_hover  = sqrt(abs(Tr_hover) / k_f);   % [rad/s] commanded speed

fprintf('=== TASK 1.2  HOVER TRIM (analytic) ===\n');
fprintf('  collective   T   = -m*g              = %+.4f N\n', T_hover);
fprintf('  per rotor    T_r = T/6               = %+.4f N\n', Tr_hover);
fprintf('  command      W_r = sqrt(|T_r|/k_f)   = %.2f rad/s = %.0f rpm\n', ...
        Omega_hover, Omega_hover*60/(2*pi));
fprintf('  moments      tau = 0  (symmetric, level, CoM at body origin)\n');
fprintf('  actuator envelope [%.0f, %.0f] rpm -> %s\n', ...
        UAV.Omega_min*60/(2*pi), UAV.Omega_max*60/(2*pi), ...
        verdict(Omega_hover >= UAV.Omega_min && Omega_hover <= UAV.Omega_max));

% ---- LEVEL 2: independent numeric solve of the same condition ----------
Th = uav.trim(UAV);
fprintf('  cross-check vs numeric solver: |dOmega| = %.2e rad/s, residual = %.2e\n', ...
        max(abs(Th.Omega - Omega_hover)), Th.residual);

% ---- LEVEL 3: open-loop propagation -----------------------------------
%  WHAT "OPEN LOOP" MEANS HERE, CONCRETELY
%  uav.propagate takes the constant rotor thrusts, converts them ONCE into a
%  constant body wrench  w_control = F1*(F2*T_r), and integrates
%
%       p_i_dot = R(q)*v_b
%       q_dot   = 0.5*W(q)*omega
%       nu_dot  = M \ ( w_g + w_aero + w_control - cross_six(nu, M*nu) )
%

Sh    = uav.propagate(UAV, Th, [0 T_END], p0 = [0;0;0]);
err_p = Sh.p;                              % commanded position is the origin

figure('Name','Hover trim position error','Color','w','Position',[80 80 780 470]);
semilogy(Sh.t, abs(err_p) + eps, 'LineWidth', 2, 'Marker','^', ...
         'MarkerIndices', 1:round(numel(Sh.t)/25):numel(Sh.t), 'MarkerSize', 5);
grid on; xlabel('t  [s]'); ylabel('|p - p_d|  [m]');
title({'Task 1.2 - hover trim held OPEN LOOP', ...
       'residual is ode45 round-off, not a control error'});
legend('x','y','z','Location','best');


%% ========================================================================
%  TASK 1.3  -  CRUISE TRIM   (LEVEL 1: ANALYTIC)
%  ========================================================================
%  Conditions: v_i = [V; 0; 0] constant. A constant inertial velocity with a
%  constant attitude implies omega = 0, so the gyroscopic term still vanishes
%  - but the body velocity v_b = R'*v_i is now NON-zero, so the drag does not.
%
%  Take level yaw and a pure pitch theta about the body y axis, R = Ry(theta):
%
%       R'*e_3 = [-sin(theta); 0; cos(theta)]
%       v_b    = R'*v_i = [ V*cos(theta); 0; -V*sin(theta) ]
%
%  Because D_fv = diag(d, d, 0) has NO vertical entry, the drag force in body
%  axes is purely along body x:
%
%       f_drag = -D_fv*v_b = [ -d*V*cos(theta); 0; 0 ]
%
%  Now write the two non-trivial force balances; the rotors contribute only
%  the term T along body z.
%
%    body x:   m*g*(-sin(theta))  -  d*V*cos(theta)   = 0                (3)
%    body z:   m*g* cos(theta)    +  T                = 0                (4)
%
%  Equation (3) contains no T at all - the HORIZONTAL balance is what fixes
%  the attitude, and it is the geometric statement "tilt until the horizontal
%  component of weight cancels the drag":
%
%       tan(theta) = -d*V/(m*g)     =>   theta = -atan( d*V/(m*g) )
%
%  Substituting into (4) gives the collective:
%
%       T = -m*g*cos(theta)
%


drag_x       = d_xy * V_CRUISE;                   % [N] drag at this speed
theta_trim   = -atan( drag_x / (m*g) );           % [rad] required pitch
T_cruise     = -m*g*cos(theta_trim);              % [N] collective
Tr_cruise    = T_cruise / Nr;
Omega_cruise = sqrt(abs(Tr_cruise) / k_f);

fprintf('\n=== TASK 1.3  CRUISE TRIM at v_i = [%.1f 0 0] m/s (analytic) ===\n', V_CRUISE);
fprintf('  drag         d*V                    = %.4f N\n', drag_x);
fprintf('  pitch        theta = -atan(dV/mg)   = %+.4f deg\n', rad2deg(theta_trim));
fprintf('               check: tan(theta) = %.5f  vs  -dV/mg = %.5f\n', ...
        tan(theta_trim), -drag_x/(m*g));
fprintf('  collective   T = -m*g*cos(theta)    = %+.4f N   (hover was %+.4f)\n', ...
        T_cruise, T_hover);
fprintf('  per rotor    T_r                    = %+.4f N\n', Tr_cruise);
fprintf('  command      W_r                    = %.2f rad/s = %.0f rpm\n', ...
        Omega_cruise, Omega_cruise*60/(2*pi));
fprintf('  actuator envelope -> %s\n', ...
        verdict(Omega_cruise >= UAV.Omega_min && Omega_cruise <= UAV.Omega_max));

% ---- LEVEL 2 -----------------------------------------------------------
Tc = uav.trim(UAV, v_i = [V_CRUISE;0;0]);
fprintf('  cross-check vs numeric solver: |dtheta| = %.2e rad, |dOmega| = %.2e rad/s\n', ...
        abs(Tc.euler(2) - theta_trim), max(abs(Tc.Omega - Omega_cruise)));

% ---- LEVEL 3 -----------------------------------------------------------
Sc    = uav.propagate(UAV, Tc, [0 T_END], p0 = [0;0;0]);
err_v = Sc.v_i - repmat([V_CRUISE 0 0], numel(Sc.t), 1);

figure('Name','Cruise trim velocity error','Color','w','Position',[120 120 780 470]);
semilogy(Sc.t, abs(err_v) + eps, 'LineWidth', 2);
grid on; xlabel('t  [s]'); ylabel('|v_i - v_{i,d}|  [m/s]');
title({sprintf('Task 1.3 - cruise trim at %.1f m/s held OPEN LOOP', V_CRUISE), ...
       sprintf('pitched %.2f deg to balance %.2f N of drag', ...
               rad2deg(theta_trim), drag_x)});
legend('v_x','v_y','v_z','Location','best');


%% ========================================================================
%  VERDICT
%  ========================================================================
fprintf('\n=== Open-loop verification over %.0f s (no feedback of any kind) ===\n', T_END);
fprintf('  hover  : max |p - p_d|     = %.3e m\n',   max(abs(err_p(:))));
fprintf('  cruise : max |v_i - v_i,d| = %.3e m/s\n', max(abs(err_v(:))));
if max(abs(err_p(:))) < 1e-6 && max(abs(err_v(:))) < 1e-6
    fprintf(['  Both sit at ode45 round-off. With no controller present, the only\n' ...
             '  explanation is that the commanded wrench cancels gravity and drag\n' ...
             '  exactly - i.e. these really are equilibria of the model.\n']);
else
    warning('Trim did not hold open loop - the algebra or the sign convention is wrong.');
end

%% ---- Save the figures ------------------------------------------------
viz.export('Task1_trim');


%% ======================= LOCAL FUNCTIONS ============================
function s = verdict(tf)
if tf, s = 'FEASIBLE'; else, s = 'INFEASIBLE'; end
end
