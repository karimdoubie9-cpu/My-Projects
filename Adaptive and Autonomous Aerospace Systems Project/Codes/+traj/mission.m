function [ref, info] = mission(opts)
%MISSION  Urban delivery mission: map -> RRT* path -> minimum-snap trajectory.
%
%   [ref, info] = traj.mission()                     % 2-D mission (default)
%   [ref, info] = traj.mission(mode="3d")            % 3-D mission, see traj.mission3d
%   [ref, info] = traj.mission(seed=1, restarts=200, patience=60)
%
%   Pipeline:
%     1. an occupancyMap3D of the buildings, inflated by the vehicle radius so
%        the planner can treat the vehicle as a point;
%     2. RRT* between the point above takeoff and the point above the pad,
%        repeated with new seeds until the best route stops improving, each run
%        followed by line-of-sight shortcutting;
%     3. minsnappolytraj over takeoff -> start -> route -> goal -> pad with
%        TimeAllocation enabled, then a uniform time stretch to the speed cap.
%        The finished trajectory is checked for building clearance and for
%        the ground, and a warning is raised if either is violated.
%
%   MODE
%     "2d"  The route is searched AT the cruise altitude (z locked in the state
%           space), so every collision check is made where the vehicle flies.
%           Every default building is taller than cruise_alt + uav_radius, so the
%           route must go around them: climb once, cruise level, descend once.
%     "3d"  z is free between min_alt and the top of the map. The planner may
%           climb over a low building or go around a tall one, whichever is
%           shorter. The route is 3-D all the way; nothing is flattened.
%
%   SEARCH EFFORT
%     MATLAB's plannerRRTStar does not improve its returned path once the goal
%     is reached (with ContinueAfterGoalReached the cost found at iteration ~70
%     is still the cost after 6000 iterations). The extra search is therefore
%     spent on RESTARTS: RRT* is run with seeds seed, seed+1, ... and the
%     shortest route after shortcutting is kept. It stops after `restarts` runs,
%     or earlier when `patience` consecutive runs have not improved the best
%     route by more than 0.1 m. Deterministic for a given seed.
%
%   info returns the map, buildings, waypoints, phase times, clearance and the
%   planning internals (info.plan) used by viz.planning.
%
%   Requires: Navigation Toolbox.
%
%   See also traj.mission3d, viz.planning.

arguments
    opts.mode       (1,:) char {mustBeMember(opts.mode, {'2d','3d'})} = '2d'
    opts.map_size   (1,3) double = [170 70 40]
    opts.resolution (1,1) double = 1.0
    opts.uav_radius (1,1) double = 2.0   % [m] required distance from the finished trajectory to any building
    opts.plan_margin (1,1) double = 0.5  % [m] extra inflation used only for PLANNING, so the
                                         %     smoothed curve can bow without dropping below uav_radius
    % Default layout (2-D mission). Six blocks on a 16 x 14 m footprint, centred
    % at (30,35) (60,55) (60,20) (90,35) (110,20) (110,55), plus B7, the delivery
    % target, whose west facade faces the landing pad 8 m away. All are taller
    % than cruise_alt + uav_radius, so the 2-D route has to go around them.
    % An empty value selects the default of the chosen mode (see traj.mission3d
    % for the 3-D layout).
    %       xmin xmax  ymin ymax  zmin zmax
    opts.buildings  (:,6) double = zeros(0, 6)
    opts.takeoff_xy (1,2) double = [10 35]
    opts.landing_xy (1,2) double = [140 35]
    opts.takeoff_alt (1,1) double = 0    % [m] height the planned trajectory starts from (a hover point)
    opts.cruise_alt (1,1) double = 20    % [m] 2-D: the flight level.
                                         %     3-D: altitude of the route's start
                                         %     and end, above takeoff and pad
    opts.min_alt    (1,1) double = 4     % [m] 3-D only: lowest altitude the route may use
    opts.speed      (1,1) double = 3.0   % [m/s] initial timing guess
    opts.max_speed  (1,1) double = 4.5   % [m/s] hard cap; see time dilation below
    opts.max_wp_spacing (1,1) double = 8 % [m] densify the route polyline to this
    opts.dt         (1,1) double = 0.01
    opts.seed       (1,1) double = 1
    opts.restarts   (1,1) double {mustBeInteger, mustBePositive} = 200
    opts.patience   (1,1) double {mustBeInteger, mustBePositive} = 60
    opts.verbose    (1,1) logical = true
end

say = @(varargin) fprintf(varargin{:});
if ~opts.verbose, say = @(varargin) []; end
is3d = strcmp(opts.mode, '3d');

if isempty(opts.buildings)
    if is3d
        opts.buildings = default_buildings_3d();
    else
        opts.buildings = default_buildings_2d();
    end
end

%% ---- 1. Occupancy map ----------------------------------------------
say('[mission 1/3] building occupancy map (%s)...\n', upper(opts.mode));
map3D = occupancyMap3D(opts.resolution);

[xg, yg, zg] = meshgrid(0:opts.resolution:opts.map_size(1), ...
                        0:opts.resolution:opts.map_size(2), ...
                        0:opts.resolution:opts.map_size(3));
pts = [xg(:) yg(:) zg(:)];
occ = false(size(pts,1), 1);

% The map is inflated by uav_radius + plan_margin. The planner therefore keeps
% the route plan_margin further out than the requirement, and the smoothed curve
% (which bows slightly between waypoints) still clears uav_radius.
r = opts.uav_radius + opts.plan_margin;
for b = 1:size(opts.buildings,1)
    B = opts.buildings(b,:);
    occ = occ | ( pts(:,1) >= B(1)-r & pts(:,1) <= B(2)+r & ...
                  pts(:,2) >= B(3)-r & pts(:,2) <= B(4)+r & ...
                  pts(:,3) >= B(5)   & pts(:,3) <= B(6)+r );
end
% Observing the whole volume (not just the obstacles) marks the free space as
% genuinely free; otherwise unobserved cells stay at the 0.5 prior and the
% state validator rejects them.
setOccupancy(map3D, pts, double(occ));

%% ---- 2. RRT*, best of several runs ---------------------------------
say('[mission 2/3] planning collision-free route (RRT*)...\n');
p_takeoff = [opts.takeoff_xy, opts.takeoff_alt];
p_land    = [opts.landing_xy, 0];
p_climb   = [opts.takeoff_xy, opts.cruise_alt];   % start of the searched route
p_descend = [opts.landing_xy, opts.cruise_alt];   % end of the searched route

if is3d
    % Free altitude between min_alt and the top of the map (minus the radius,
    % so the vehicle never grazes the ceiling of the map).
    z_bounds = [opts.min_alt, opts.map_size(3) - r];
else
    % Locked to the flight level: every sample and every collision check is
    % made at the altitude actually flown, so no flattening step is needed.
    z_bounds = [opts.cruise_alt, opts.cruise_alt];
end
ss = stateSpaceSE3([0 opts.map_size(1); 0 opts.map_size(2); z_bounds; ...
                    -1 1; -1 1; -1 1; -1 1]);
ss.WeightQuaternion = 0;   % path cost = length only; the orientation part of SE(3) is unused
sv = validatorOccupancyMap3D(ss, Map = map3D);
% Edges are checked every 0.25 m (a quarter of the map resolution). At 1 m a
% route that hugs an inflated corner can have its samples straddle an occupied
% cell -- the edge passes, and a sub-segment of it later fails.
sv.ValidationDistance = opts.resolution / 4;

start_state = [p_climb,   1 0 0 0];
goal_state  = [p_descend, 1 0 0 0];
if ~isStateValid(sv, start_state) || ~isStateValid(sv, goal_state)
    error('traj:mission:invalidEndpoint', ...
          'Start or goal lies inside the inflated obstacle set (or outside the altitude band).');
end

% MaxConnectionDistance = 5 m: short edges give the tree many vertices, including
% some near the building corners, so line-of-sight shortcutting (which can only
% connect existing vertices) has good points to pull the route taut against.
% With 15 m edges the shortcut route kept visible detours, e.g. swinging wide of
% B1 because no vertex lay near its corner. Short edges cost more iterations to
% reach the goal; the 10000-iteration budget is ample (no failed runs in either
% mission).
planner = plannerRRTStar(ss, sv, MaxConnectionDistance = 5, ...
                                 MaxIterations = 10000, GoalBias = 0.12);

best  = struct('len', inf);
lens  = nan(1, opts.restarts);
stale = 0;
for k = 1:opts.restarts
    rng(opts.seed + k - 1);
    [pth_k, info_k] = plan(planner, start_state, goal_state);
    if pth_k.NumStates == 0, continue, end
    wp_k    = shortcut(pth_k.States(:,1:3), sv);
    lens(k) = sum(vecnorm(diff(wp_k, 1, 1), 2, 2));
    if lens(k) < best.len - 0.1
        best  = struct('len', lens(k), 'seed', opts.seed + k - 1, ...
                       'pth', pth_k, 'info', info_k, 'wp', wp_k);
        stale = 0;
    else
        stale = stale + 1;
        if stale >= opts.patience, break, end
    end
end
n_runs = k;
if ~isfinite(best.len)
    error('traj:mission:noPath', 'RRT* found no path in %d runs; check geometry/clearance.', n_runs);
end
pth = best.pth;  solnInfo = best.info;  route = best.wp;
lens = lens(1:n_runs);
say('  RRT*: %d runs, best is seed %d at %.1f m after shortcutting (runs ranged %.1f-%.1f m, straight line %.1f m)\n', ...
    n_runs, best.seed, best.len, min(lens), max(lens), norm(p_descend - p_climb));
% 3-D only: climb vertically to 2 m below the route's first waypoint (and start
% the vertical descent 2 m below its last one), instead of to the fixed search
% start/end altitude. The route then leaves the vertical leg almost level, so
% the curve does not have to turn a vertical climb into a steep diagonal one
% right after takeoff. Both modified segments are re-validated; if either is
% not collision-free the original endpoints are kept.
if is3d && size(route,1) >= 3
    q0      = [1 0 0 0];
    new_up  = [opts.takeoff_xy, max(opts.min_alt, route(2,3)     - 2)];
    new_dn  = [opts.landing_xy, max(opts.min_alt, route(end-1,3) - 2)];
    ok = isStateValid(sv, [new_up q0]) && isStateValid(sv, [new_dn q0]) && ...
         isMotionValid(sv, [new_up q0], [route(2,:) q0]) && ...
         isMotionValid(sv, [route(end-1,:) q0], [new_dn q0]);
    if ok
        say('  vertical legs: climb to %.1f m, descend from %.1f m (2 m below the first/last route waypoint)\n', ...
            new_up(3), new_dn(3));
        route(1,:)   = new_up;   p_climb   = new_up;
        route(end,:) = new_dn;   p_descend = new_dn;
    else
        say('  vertical legs: raised endpoints not collision-free, keeping %.1f m\n', opts.cruise_alt);
    end
end
shortcut_wp = route;                     % kept for info.plan (before densification)

%% ---- 3. Minimum-snap trajectory ------------------------------------
say('[mission 3/3] computing minimum-snap trajectory...\n');

% minsnappolytraj only constrains the curve AT the waypoints; between two distant
% ones it is free to bow sideways, which is how a collision-free waypoint list
% can still give a trajectory that shaves an inflated footprint. Densifying pins
% the curve to the route.
dense = densify(route, opts.max_wp_spacing);
spacing = opts.max_wp_spacing;
% Every segment of the densified route must be collision-free: they are
% sub-segments of shortcut edges that were checked, so this can only fail if
% something upstream is wrong -- in which case stop rather than fly it.
for i = 1:size(dense,1)-1
    if ~isMotionValid(sv, [dense(i,:) 1 0 0 0], [dense(i+1,:) 1 0 0 0])
        error('traj:mission:invalidRoute', ...
              'Route segment %d (%s -> %s) is not collision-free.', ...
              i, mat2str(dense(i,:),3), mat2str(dense(i+1,:),3));
    end
end
T = solve_trajectory([p_takeoff; p_climb; dense; p_descend; p_land], opts);
if T.k > 1
    say('  speed cap: peak was %.2f m/s, dilating time by %.2fx -> %.2f m/s\n', ...
        T.peak, T.k, max(vecnorm(T.v,2,1)));
end

% Checks on the FINISHED trajectory (reported, not repaired).
clr   = clearance(T.p, opts.buildings);
z_min = min(T.p(3,:));
if clr < opts.uav_radius
    warning('traj:mission:clearance', ...
        ['Planned trajectory passes %.2f m from a building, inside the %.1f m ' ...
         'inflation. Reduce max_wp_spacing or raise uav_radius.'], clr, opts.uav_radius);
end
if z_min < -0.01
    warning('traj:mission:ground', 'Planned trajectory goes %.2f m below the ground.', -z_min);
end

t = T.t; p = T.p; v = T.v; a = T.a;

% Heading follows the horizontal velocity, holding the last valid value while
% nearly stationary (atan2 is meaningless at takeoff/landing where |v_xy| ~ 0).
yaw = zeros(size(t)); last = 0;
for i = 1:numel(t)
    if norm(v(1:2,i)) > 0.12, last = atan2(v(2,i), v(1,i)); end
    yaw(i) = last;
end
yaw = unwrap(yaw);

ref = traj.pack(t(:), p, v, a, yaw);

%% ---- info ----------------------------------------------------------
[~, i_climb]   = min(vecnorm(p - p_climb.',   2, 1));
[~, i_descend] = min(vecnorm(p - p_descend.', 2, 1));

info.mode       = opts.mode;
info.map        = map3D;
info.buildings  = opts.buildings;
info.map_size   = opts.map_size;
info.uav_radius = opts.uav_radius;   % so viz.map2d draws the true inflation
info.cruise_alt = opts.cruise_alt;
info.waypoints  = T.wp;
info.cruise_wp  = dense;
info.p_takeoff  = p_takeoff;
info.p_land     = p_land;
info.p_climb    = p_climb;
info.p_descend  = p_descend;
info.t_climb    = t(i_climb);
info.t_descend  = t(i_descend);
info.phases     = [t(i_climb) t(i_descend)];
info.plan.pre         = T.pre;                        % before the time stretch
info.plan.tree        = solnInfo.TreeData(:, 1:3);   % NaN-separated edges (best run)
info.plan.raw         = pth.States(:, 1:3);          % RRT* path, before shortcutting
info.plan.raw_quat    = pth.States(:, 4:7);          % orientation part of SE(3): unused
info.plan.shortcut    = shortcut_wp;                 % after shortcutting
info.plan.dense       = dense;                       % after densification
info.plan.max_spacing = spacing;                     % spacing actually used
info.plan.wp          = T.wp;                        % list handed to minsnappolytraj
info.plan.t_wp_init   = T.timepoints;                % distance / speed guess
info.plan.t_wp_opt    = T.tWp;                       % after time allocation
info.plan.speed       = opts.speed;
info.plan.k           = max(T.k, 1);
info.plan.runs        = n_runs;                      % RRT* restarts performed
info.plan.run_lengths = lens;                        % route length of each run [m]
info.plan.best_seed   = best.seed;
info.plan.route_length = best.len;
info.T_end      = t(end);
info.max_speed  = opts.max_speed;         % the cap the profile was scaled to
info.peak_speed = max(vecnorm(v, 2, 1));  % what it actually reaches

% Clearance of the FINAL trajectory, not of the waypoints: the smoothing happens
% after planning, and only this describes what the vehicle is asked to fly.
info.clearance = clr;
info.z_min     = z_min;               % lowest point of the trajectory (ground check)
info.peak_acc  = max(vecnorm(a, 2, 1));
say('  min clearance to any building: %.2f m (inflation %.1f m)\n', clr, opts.uav_radius);
say(['mission: takeoff 0-%.1f s | transit %.1f-%.1f s | landing %.1f-%.1f s ' ...
     '(total %.1f s)\n'], info.t_climb, info.t_climb, info.t_descend, ...
     info.t_descend, info.T_end, info.T_end);
end

% ---------------------------------------------------------------------
function T = solve_trajectory(wp, opts)
%SOLVE_TRAJECTORY  Minimum snap with time allocation, then the speed-cap stretch.

% Drop near-coincident waypoints: they add no geometry but create a segment
% whose duration is a small fraction of a second.
keep      = [true; vecnorm(diff(wp,1,1),2,2) > 0.5];
keep(end) = true;                       % never drop the landing point
wp        = wp(keep, :);

seg        = vecnorm(diff(wp,1,1), 2, 2) / opts.speed;
timepoints = [0; cumsum(seg)].';
nSamples   = max(200, round(timepoints(end)/opts.dt));

% The segment-time bounds must BRACKET every initial guess, or minsnappolytraj
% rejects the problem; deriving them from the actual spread keeps it feasible
% for any geometry.
minSeg = max(1e-3, 0.5*min(seg));
maxSeg = 1.5*max(seg) + 5;

[p, v, a, ~, ~, ~, tWp, tS] = minsnappolytraj(wp.', timepoints, nSamples, ...
    MinSegmentTime = minSeg, ...
    MaxSegmentTime = maxSeg, ...
    TimeAllocation = true, TimeWeight = 80);
t = tS(:).';
pre = struct('t', t, 'p', p, 'v', v, 'a', a);

% Uniform time dilation: t -> k*t leaves the path unchanged and gives
% v -> v/k, a -> a/k^2, so one factor brings the whole profile under the cap.
peak = max(vecnorm(v, 2, 1));
k    = peak / opts.max_speed;
if k > 1
    t = t*k;  v = v/k;  a = a/k^2;
end
T = struct('t', t, 'p', p, 'v', v, 'a', a, 'k', k, 'peak', peak, 'pre', pre, ...
           'wp', wp, 'timepoints', timepoints, 'tWp', tWp);
end

function B = default_buildings_2d()
%       xmin xmax  ymin ymax  zmin zmax
B = [    22   38    28   42     0   28     % B1 (30,35)
         52   68    48   62     0   26     % B2 (60,55)
         52   68    13   27     0   32     % B3 (60,20)
         76   92    32   46     0   32     % B4 (90,35)
        102  118    13   27     0   24     % B5 (110,20)
        102  118    48   62     0   30     % B6 (110,55)
        148  164    27   43     0   30 ];  % B7 (156,35) delivery target
end

function B = default_buildings_3d()
% Same footprints, heights varied on purpose: B1 and B4 sit on the direct line
% from takeoff to the pad and are LOW (overflying them is shorter than going
% around), B3 is tall enough that nothing can pass over it, and the rest are
% in between. The planner has to trade climbing against detouring.
%       xmin xmax  ymin ymax  zmin zmax
B = [    22   38    28   42     0   12     % B1 (30,35)  low - on the direct line
         52   68    48   62     0   26     % B2 (60,55)
         52   68    13   27     0   36     % B3 (60,20)  tallest - cannot be overflown
         82   98    28   42     0    40     % B4 (90,35)  very low - on the direct line
        102  118    13   27     0   18     % B5 (110,20)
        102  118    48   62     0   30     % B6 (110,55)
        148  164    27   43     0   30 ];  % B7 (156,35) delivery target
end

function wp = densify(wp, maxStep)
%DENSIFY  Insert collinear waypoints so no segment exceeds maxStep metres.
out = wp(1,:);
for i = 1:size(wp,1)-1
    a = wp(i,:); b = wp(i+1,:);
    n = max(1, ceil(norm(b-a)/maxStep));
    for k = 1:n
        out(end+1,:) = a + (b-a)*(k/n);   %#ok<AGROW>
    end
end
wp = out;
end

function d = clearance(p, buildings)
%CLEARANCE  Smallest distance from the sampled path to any building box.
%   Zero inside a box, otherwise the Euclidean distance to its surface.
d = inf;
for b = 1:size(buildings,1)
    B  = buildings(b,:);
    dx = max([B(1)-p(1,:); p(1,:)-B(2); zeros(1,size(p,2))]);
    dy = max([B(3)-p(2,:); p(2,:)-B(4); zeros(1,size(p,2))]);
    dz = max([B(5)-p(3,:); p(3,:)-B(6); zeros(1,size(p,2))]);
    d  = min(d, min(sqrt(dx.^2 + dy.^2 + dz.^2)));
end
end

function wp = shortcut(raw, sv)
%SHORTCUT  Line-of-sight path shortening using the validator's motion check.
%   Greedily connects each retained waypoint to the furthest later waypoint it
%   can reach through free space. Every accepted segment is checked, including
%   the fallback to the very next point, so a shortcut can never introduce a
%   segment that clips an obstacle.
if size(raw,1) <= 2, wp = raw; return; end
wp   = raw(1,:);
curr = 1;
while curr < size(raw,1)
    next = size(raw,1);
    while next > curr + 1 && ~isMotionValid(sv, [raw(curr,:) 1 0 0 0], [raw(next,:) 1 0 0 0])
        next = next - 1;
    end
    wp(end+1,:) = raw(next,:);   %#ok<AGROW>
    curr = next;
end
end
