function [ref, info] = mission3d(varargin)
%MISSION3D  3-D delivery mission: buildings of varying heights, free altitude.
%
%   [ref, info] = traj.mission3d()
%   [ref, info] = traj.mission3d(seed=2, restarts=80, cruise_alt=10)
%
%   Same pipeline as traj.mission, with mode="3d": the RRT* search is not
%   locked to a flight level. The ALTITUDE PROFILE IS CHOSEN BY THE PLANNER:
%   the route starts and ends cruise_alt (10 m) above the takeoff point and
%   the pad -- traj.mission then raises these endpoints to 2 m below the first
%   and last route waypoints -- and in between may use any altitude from
%   min_alt (4 m) to the top of the map. It can therefore climb over the low
%   buildings (B1, roof 12 m) or go around the tall ones (B4, roof 40 m),
%   whichever gives the shorter route. Nothing is flattened. With the default
%   map it flies OVER B1 (~15.5 m) and passes north of B4 (~18 m).
%
%   Why the route starts at 10 m and not lower: with a 4 m vertical leg
%   followed directly by a steep climb, the minimum-snap curve dips about 6 m
%   below the ground right after takeoff (traj.mission then warns).
%
%   Any traj.mission option can be passed through (seed, restarts, patience,
%   buildings, max_speed, ...).
%
%   See also traj.mission.

[ref, info] = traj.mission(varargin{:}, mode = '3d', cruise_alt = alt_default(varargin));
end

function z = alt_default(args)
% Respect a cruise_alt passed by the caller; otherwise start and end the
% route 10 m above takeoff and pad (see the note in the help text).
z = 10;
for i = 1:2:numel(args)-1
    if (ischar(args{i}) || isstring(args{i})) && strcmpi(args{i}, 'cruise_alt')
        z = args{i+1};
    end
end
end
