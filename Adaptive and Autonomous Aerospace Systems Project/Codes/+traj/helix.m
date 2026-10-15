function ref = helix_min_snap(T, opts)
%HELIX_MIN_SNAP  Free-entry minimum-snap approach followed by a helix.
%
%   ref = traj.helix_min_snap(T, ...)
%
% The vehicle starts at opts.start_p at rest and reaches the helix after
% opts.approach_time. The entry phase on the helix is NOT prescribed.
% Instead, the function searches over the helix phase and selects the
% admissible entry that gives the smallest translational snap cost
%
%       J = integral_0^Ta || p''''(t) ||^2 dt
%
% for a single 7th-order polynomial approach.
%
% The approach matches position, velocity, acceleration AND jerk at the
% helix join. Therefore p, v, a and j are continuous at t = Ta.
%
% By default, phases that would reverse the horizontal heading direction
% during the approach are rejected. This avoids introducing an artificial
% yaw reversal into a Task-3 controller comparison.
%
% The returned reference has the same interface as traj.helix:
% time, position, velocity, acceleration and yaw are packed with traj.pack.

arguments
    T                       (1,1) double {mustBePositive}
    opts.radius             (1,1) double {mustBePositive} = 5
    opts.rate               (1,1) double = 0.1
    opts.climb              (1,1) double = 2.0
    opts.z0                 (1,1) double = 10
    opts.dt                 (1,1) double {mustBePositive} = 0.01

    opts.center             (2,1) double = [2; 2]
    opts.start_p            (3,1) double = [0; 0; 0]
    opts.approach_time      (1,1) double {mustBePositive} = 8

    % Free-entry search.
    opts.phase_samples      (1,1) double {mustBePositive, mustBeInteger} = 1440

    % Keep the velocity-defined yaw from changing rotation direction during
    % the approach. This is useful here because Task 3 is intended to test
    % the plant/controller mismatch, not an arbitrary heading reversal.
    opts.enforce_monotonic_yaw (1,1) logical = true
    opts.yaw_speed_epsilon  (1,1) double {mustBeNonnegative} = 1e-3
    opts.yaw_rate_tolerance (1,1) double {mustBeNonnegative} = 1e-6

    % Diagnostic only. The exam trajectory guideline is 4-5 m/s.
    opts.max_speed          (1,1) double {mustBePositive} = 5
    opts.verbose            (1,1) logical = false
end

dt = opts.dt;
Ta = opts.approach_time;

if Ta >= T
    error('approach_time must be smaller than total duration T.');
end

if opts.phase_samples < 72
    error('phase_samples should be at least 72 for a meaningful free-entry search.');
end

R  = opts.radius;
w  = opts.rate;
vz = opts.climb;
c  = opts.center;
p0 = opts.start_p;

%% ========================================================================
%  1. FREE HELIX-ENTRY PHASE
%  ========================================================================

Nphi = opts.phase_samples;
phi_grid = -pi + (0:Nphi-1) * (2*pi/Nphi);

J = inf(size(phi_grid));
feasible = false(size(phi_grid));

for i = 1:Nphi
    phi = phi_grid(i);

    [C_i, J_i] = approach_coefficients( ...
        phi, Ta, p0, c, R, w, vz, opts.z0);

    J(i) = J_i;

    if opts.enforce_monotonic_yaw
        feasible(i) = no_yaw_reversal( ...
            C_i, Ta, w, ...
            opts.yaw_speed_epsilon, ...
            opts.yaw_rate_tolerance);
    else
        feasible(i) = true;
    end
end

if any(feasible)
    J_select = J;
    J_select(~feasible) = inf;
    [~, idx] = min(J_select);
else
    warning(['No phase satisfied the monotonic-yaw screening. ' ...
             'Falling back to the unconstrained minimum-snap entry phase.']);
    [~, idx] = min(J);
end

phiE = phi_grid(idx);

[C, Jstar] = approach_coefficients( ...
    phiE, Ta, p0, c, R, w, vz, opts.z0);

%% ========================================================================
%  2. MINIMUM-SNAP APPROACH
%
%  Single 7th-order polynomial:
%
%      p(t) = c0 + c1 t + ... + c7 t^7
%
%  Boundary conditions:
%      t = 0  : p = start_p, v = 0, a = 0, j = 0
%      t = Ta : p,v,a,j exactly match the helix at phiE
%
%  For these eight endpoint constraints the degree-7 polynomial is unique,
%  and it is the single-segment minimum-snap trajectory.
%  ========================================================================

ta = 0:dt:Ta;

[pA, vA, aA, jA] = evaluate_poly7(C, ta);

% Exact boundary residuals at the join. This guards against constructing the
% polynomial in absolute rather than displacement coordinates.
pE_check = helix_position(phiE, c, R, opts.z0);
vE_check = helix_velocity(phiE, R, w, vz);
aE_check = helix_acceleration(phiE, R, w);
jE_check = helix_jerk(phiE, R, w);

join_residual = [ ...
    norm(pA(:,end) - pE_check), ...
    norm(vA(:,end) - vE_check), ...
    norm(aA(:,end) - aE_check), ...
    norm(jA(:,end) - jE_check)];

if any(join_residual > 1e-8)
    error('helix_min_snap:joinMismatch', ...
        ['Minimum-snap boundary mismatch: |Dp|=%.3e, |Dv|=%.3e, ' ...
         '|Da|=%.3e, |Dj|=%.3e'], join_residual);
end

%% ========================================================================
%  3. HELIX AFTER THE JOIN
%  ========================================================================

th = (Ta + dt):dt:T;
trel = th - Ta;

phi = phiE + w*trel;

pH = [c(1) + R*cos(phi);
      c(2) + R*sin(phi);
      opts.z0 + vz*trel];

vH = [-R*w*sin(phi);
       R*w*cos(phi);
       vz*ones(size(th))];

aH = [-R*w^2*cos(phi);
      -R*w^2*sin(phi);
       zeros(size(th))];

%% ========================================================================
%  4. MERGE
%  ========================================================================

t = [ta, th];
p = [pA, pH];
v = [vA, vH];
a = [aA, aH];

%% ========================================================================
%  5. HEADING
%
%  Keep heading aligned with horizontal velocity, as in the original helix.
%  The free-entry search rejects phases that introduce a yaw-rate reversal
%  during the approach (unless enforce_monotonic_yaw=false).
%  ========================================================================

yaw = zeros(size(t));

speed_xy = hypot(v(1,:), v(2,:));
moving = speed_xy > opts.yaw_speed_epsilon;

yaw(moving) = atan2(v(2,moving), v(1,moving));

firstMoving = find(moving, 1, 'first');
if ~isempty(firstMoving)
    yaw(1:firstMoving-1) = yaw(firstMoving);
end

yaw = unwrap(yaw);

%% ========================================================================
%  6. SANITY CHECKS
%  ========================================================================

speed = vecnorm(v, 2, 1);
vmax = max(speed);

if vmax > opts.max_speed + 1e-9
    warning('helix_min_snap:maxSpeed', ...
        'Reference reaches %.3f m/s, above the requested %.3f m/s limit.', ...
        vmax, opts.max_speed);
end

if opts.verbose
    pE = helix_position(phiE, c, R, opts.z0);
    vE = helix_velocity(phiE, R, w, vz);

    fprintf('\n=== FREE-ENTRY MINIMUM-SNAP HELIX ===\n');
    fprintf('approach time         : %.3f s\n', Ta);
    fprintf('selected entry phase  : %.3f deg\n', rad2deg(phiE));
    fprintf('entry position [m]    : [%+.3f %+.3f %+.3f]\n', pE);
    fprintf('entry velocity [m/s]  : [%+.3f %+.3f %+.3f]\n', vE);
    fprintf('snap cost             : %.6g\n', Jstar);
    fprintf('max reference speed   : %.3f m/s\n', vmax);
    fprintf('join residuals        : |Dp| %.3e, |Dv| %.3e, |Da| %.3e, |Dj| %.3e\n', ...
        join_residual(1), join_residual(2), join_residual(3), join_residual(4));
    fprintf('yaw reversal screening: %d\n', opts.enforce_monotonic_yaw);
end

ref = traj.pack(t(:), p, v, a, yaw);

end


% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function [C, J] = approach_coefficients(phi, T, p0, c, R, w, vz, z0)
% 7th-order polynomial coefficients satisfying p,v,a,j at both endpoints.

pE = helix_position(phi, c, R, z0);
vE = helix_velocity(phi, R, w, vz);
aE = helix_acceleration(phi, R, w);
jE = helix_jerk(phi, R, w);

% C(:,k+1) multiplies t^k.
C = zeros(3,8);

% Start from rest with zero acceleration and zero jerk.
C(:,1) = p0;
C(:,2) = 0;
C(:,3) = 0;
C(:,4) = 0;

A = [ ...
    T^4,       T^5,        T^6,         T^7;
    4*T^3,     5*T^4,      6*T^5,       7*T^6;
    12*T^2,    20*T^3,     30*T^4,      42*T^5;
    24*T,      60*T^2,     120*T^3,     210*T^4];

for k = 1:3
    rhs = [pE(k) - p0(k); vE(k); aE(k); jE(k)];
    C(k,5:8) = (A \ rhs).';
end

% Exact integral of squared snap.
% snap(t) = q0 + q1*t + q2*t^2 + q3*t^3
Q = [ ...
    24*C(:,5), ...
    120*C(:,6), ...
    360*C(:,7), ...
    840*C(:,8)];

G = zeros(4,4);
for i = 0:3
    for j = 0:3
        G(i+1,j+1) = T^(i+j+1)/(i+j+1);
    end
end

J = trace(Q * G * Q.');

end


function ok = no_yaw_reversal(C, T, w, speed_eps, yaw_rate_tol)
% Reject an approach whose horizontal velocity heading rotates opposite to
% the final helix rotation direction.

if abs(w) < eps
    ok = true;
    return
end

tt = linspace(0, T, 301);
[~, v, a, ~] = evaluate_poly7(C, tt);

vx = v(1,:);
vy = v(2,:);
ax = a(1,:);
ay = a(2,:);

s2 = vx.^2 + vy.^2;
moving = sqrt(s2) > speed_eps;

if ~any(moving)
    ok = false;
    return
end

yaw_rate = (vx(moving).*ay(moving) - vy(moving).*ax(moving)) ./ s2(moving);

ok = all(sign(w) * yaw_rate >= -yaw_rate_tol);

end


function [p, v, a, j] = evaluate_poly7(C, t)

t = reshape(t,1,[]);

p = zeros(3,numel(t));
v = zeros(3,numel(t));
a = zeros(3,numel(t));
j = zeros(3,numel(t));

for n = 0:7
    p = p + C(:,n+1) .* t.^n;
end

for n = 1:7
    v = v + n*C(:,n+1) .* t.^(n-1);
end

for n = 2:7
    a = a + n*(n-1)*C(:,n+1) .* t.^(n-2);
end

for n = 3:7
    j = j + n*(n-1)*(n-2)*C(:,n+1) .* t.^(n-3);
end

end


function p = helix_position(phi, c, R, z0)
p = [c(1) + R*cos(phi);
     c(2) + R*sin(phi);
     z0];
end


function v = helix_velocity(phi, R, w, vz)
v = [-R*w*sin(phi);
      R*w*cos(phi);
      vz];
end


function a = helix_acceleration(phi, R, w)
a = [-R*w^2*cos(phi);
     -R*w^2*sin(phi);
      0];
end


function j = helix_jerk(phi, R, w)
j = [ R*w^3*sin(phi);
     -R*w^3*cos(phi);
      0];
end
