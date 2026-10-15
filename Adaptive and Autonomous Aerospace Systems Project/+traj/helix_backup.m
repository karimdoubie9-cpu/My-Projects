% function ref = helix(T, opts)
% %HELIX  Helicoidal reference: constant-radius turn with a constant climb rate.
% %
% %   ref = traj.helix(150)
% %   ref = traj.helix(150, radius=5, rate=0.1, climb=0.5, dt=0.01)
% %
% %   This is the standard trajectory-tracking scenario for Tasks 2.3 and 3: it
% %   excites all three translational axes simultaneously and keeps the heading
% %   turning, so yaw control is exercised too, at a speed well inside the 4-5 m/s
% %   limit set by the exam text.
% %
% %   Default speed check:  |v| = sqrt((R*w)^2 + vz^2)
% %                             = sqrt(0.5^2 + 0.5^2) = 0.71 m/s.
% %
% %   Heading is aligned with the horizontal velocity and unwrapped, so the
% %   attitude planner sees a continuous yaw command rather than +-pi jumps.
% 
% arguments
%     T           (1,1) double {mustBePositive}
%     opts.radius (1,1) double = 5      % [m]     helix radius
%     opts.rate   (1,1) double = 0.1    % [rad/s] angular rate about the vertical
%     opts.climb  (1,1) double = 1    % [m/s]   constant climb rate
%     opts.z0     (1,1) double = 0      % [m]     initial altitude
%     opts.dt     (1,1) double = 0.01   % [s]     sample interval
% end
% 
% t   = (0:opts.dt:T);
% R   = opts.radius;
% w   = opts.rate;
% phi = w*t;
% 
% p = [ R*cos(phi);
%       R*sin(phi);
%       opts.z0 + opts.climb*t ];
% 
% v = [-R*w*sin(phi);
%       R*w*cos(phi);
%       opts.climb*ones(size(t)) ];
% 
% a = [-R*w^2*cos(phi);
%      -R*w^2*sin(phi);
%       zeros(size(t)) ];
% 
% yaw = unwrap(atan2(v(2,:), v(1,:)));
% 
% ref = traj.pack(t(:), p, v, a, yaw);
% end
function ref = helix(T, opts)
%HELIX  Smooth approach from rest followed by a helical trajectory.
%
% T is the TOTAL simulation duration.
%
% The vehicle starts at opts.start_p with zero velocity/acceleration,
% approaches the helix smoothly during opts.approach_time, and then joins
% the helix with continuous position, velocity and acceleration.

arguments
    T                  (1,1) double {mustBePositive}
    opts.radius        (1,1) double = 5
    opts.rate          (1,1) double = 0.1
    opts.climb         (1,1) double = 2.0
    opts.z0            (1,1) double = 10
    opts.dt            (1,1) double = 0.01

    % NEW
    opts.center        (2,1) double = [2; 2]    % helix axis in xy [m]
    opts.phase0        (1,1) double = pi        % entry angle [rad]
    opts.start_p       (3,1) double = [0;0;0]   % vehicle starts here
    opts.approach_time (1,1) double = 8         % [s]
end

dt = opts.dt;
Ta = opts.approach_time;

if Ta >= T
    error('approach_time must be smaller than total duration T.');
end

R   = opts.radius;
w   = opts.rate;
vz  = opts.climb;
phi0 = opts.phase0;
c   = opts.center;

%% =====================================================================
%  HELIX STATE AT THE JOIN
% ======================================================================

% Position at start of helix
pH0 = [c(1) + R*cos(phi0);
       c(2) + R*sin(phi0);
       opts.z0];

% Velocity at start of helix
vH0 = [-R*w*sin(phi0);
        R*w*cos(phi0);
        vz];

% Acceleration at start of helix
aH0 = [-R*w^2*cos(phi0);
       -R*w^2*sin(phi0);
        0];

%% =====================================================================
%  1. SMOOTH APPROACH
% ======================================================================

ta = (0:dt:Ta);
tau = ta / Ta;

% Quintic Hermite polynomial:
% matches p,v,a at BOTH endpoints.
%
% p(0)  = start_p
% v(0)  = 0
% a(0)  = 0
%
% p(Ta) = pH0
% v(Ta) = vH0
% a(Ta) = aH0

p0 = opts.start_p;
v0 = zeros(3,1);
a0 = zeros(3,1);

% polynomial coefficients p = c0+c1*t+c2*t^2+...+c5*t^5
C = zeros(3,6);

C(:,1) = p0;
C(:,2) = v0;
C(:,3) = a0/2;

A = [Ta^3,      Ta^4,       Ta^5;
     3*Ta^2,    4*Ta^3,     5*Ta^4;
     6*Ta,     12*Ta^2,    20*Ta^3];

for k = 1:3
    rhs = [pH0(k) - C(k,1) - C(k,2)*Ta - C(k,3)*Ta^2;
           vH0(k) - C(k,2) - 2*C(k,3)*Ta;
           aH0(k) - 2*C(k,3)];

    C(k,4:6) = A \ rhs;
end

pA = zeros(3,numel(ta));
vA = zeros(3,numel(ta));
aA = zeros(3,numel(ta));

for k = 1:3
    pA(k,:) = C(k,1) ...
            + C(k,2)*ta ...
            + C(k,3)*ta.^2 ...
            + C(k,4)*ta.^3 ...
            + C(k,5)*ta.^4 ...
            + C(k,6)*ta.^5;

    vA(k,:) = C(k,2) ...
            + 2*C(k,3)*ta ...
            + 3*C(k,4)*ta.^2 ...
            + 4*C(k,5)*ta.^3 ...
            + 5*C(k,6)*ta.^4;

    aA(k,:) = 2*C(k,3) ...
            + 6*C(k,4)*ta ...
            + 12*C(k,5)*ta.^2 ...
            + 20*C(k,6)*ta.^3;
end

%% =====================================================================
%  2. HELIX
% ======================================================================

th = (Ta+dt : dt : T);
trel = th - Ta;

phi = phi0 + w*trel;

pH = [c(1) + R*cos(phi);
      c(2) + R*sin(phi);
      opts.z0 + vz*trel];

vH = [-R*w*sin(phi);
       R*w*cos(phi);
       vz*ones(size(th))];

aH = [-R*w^2*cos(phi);
      -R*w^2*sin(phi);
       zeros(size(th))];

%% =====================================================================
%  MERGE
% ======================================================================

t = [ta, th];
p = [pA, pH];
v = [vA, vH];
a = [aA, aH];

% Heading from horizontal velocity
yaw = zeros(size(t));

speed_xy = hypot(v(1,:),v(2,:));
moving = speed_xy > 1e-8;

yaw(moving) = atan2(v(2,moving),v(1,moving));

% At exactly t=0 velocity is zero: use first meaningful direction.
firstMoving = find(moving,1,'first');
if ~isempty(firstMoving)
    yaw(1:firstMoving-1) = yaw(firstMoving);
end

yaw = unwrap(yaw);

ref = traj.pack(t(:), p, v, a, yaw);
end