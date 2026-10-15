function fig = map2d(M, opts)
%MAP2D  Top-down mission view: footprints, inflation margin, waypoints, paths.
%
%   viz.map2d(M)                                  % plan only
%   viz.map2d(M, runs=runs, tags=tags)            % plan + flown paths
%   viz.map2d(M, inflation=2.0)
%
%   The dashed red rectangles are the obstacle set the PLANNER actually saw:
%   each building grown by the vehicle clearance radius. Drawing only the true
%   footprints makes the reference look needlessly wide of the buildings; the
%   inflated boundary is what it is really hugging, and showing it is what makes
%   the clearance argument legible.
%
%   M is the info struct from traj.mission.

arguments
    M              struct
    opts.runs      cell = {}
    opts.tags      cell = {}
    opts.ref       struct = struct([])
    opts.inflation (1,1) double = NaN
    opts.name      (1,:) char = 'Top-down flight path'
end

infl = opts.inflation;
if isnan(infl)
    if isfield(M,'uav_radius'), infl = M.uav_radius; else, infl = 2.0; end
end

s   = viz.style();
fig = figure('Name', opts.name, 'Color','w', 'Position',[100 100 980 580]);
hold on; grid on; axis equal

for b = 1:size(M.buildings,1)
    B = M.buildings(b,:);
    rectangle('Position',[B(1) B(3) B(2)-B(1) B(4)-B(3)], ...
              'FaceColor',[0.60 0.60 0.65], 'EdgeColor','k');
    rectangle('Position',[B(1)-infl B(3)-infl ...
                          (B(2)-B(1))+2*infl (B(4)-B(3))+2*infl], ...
              'EdgeColor',[0.85 0.20 0.20], 'LineStyle','--');
end
% Legend proxies for the two rectangle styles (rectangle() takes no DisplayName)
plot(NaN, NaN, 's', 'MarkerFaceColor',[0.6 0.6 0.65], 'MarkerEdgeColor','k', ...
     'DisplayName','Buildings');
plot(NaN, NaN, '--', 'Color',[0.85 0.2 0.2], ...
     'DisplayName', sprintf('Inflated by %.1f m', infl));

if isfield(M,'cruise_wp') && ~isempty(M.cruise_wp)
    plot(M.cruise_wp(:,1), M.cruise_wp(:,2), 'bo', 'MarkerFaceColor','b', ...
         'MarkerSize',6, 'DisplayName','Waypoints');
end

if ~isempty(fieldnames(opts.ref))
    p = traj.pos(opts.ref);
    plot(p(1,:), p(2,:), '--k', 'LineWidth',1.8, 'DisplayName','Planned reference');
elseif ~isempty(opts.runs)
    plot(opts.runs{1}.p_d(:,1), opts.runs{1}.p_d(:,2), '--k', ...
         'LineWidth',1.8, 'DisplayName','Planned reference');
end

for k = 1:numel(opts.runs)
    plot(opts.runs{k}.p(:,1), opts.runs{k}.p(:,2), 'LineWidth',2.1, ...
         'Color', s.colors{1+mod(k-1,numel(s.colors))}, 'DisplayName', opts.tags{k});
end

plot(M.p_takeoff(1), M.p_takeoff(2), 'r^', 'MarkerSize',12, ...
     'MarkerFaceColor','r', 'DisplayName','Takeoff');
plot(M.p_land(1), M.p_land(2), 'gv', 'MarkerSize',12, ...
     'MarkerFaceColor','g', 'DisplayName','Landing');

xlabel('x  [m]'); ylabel('y  [m]');
xlim([0 M.map_size(1)]); ylim([0 M.map_size(2)]);
title(opts.name); legend('Location','eastoutside');
end
