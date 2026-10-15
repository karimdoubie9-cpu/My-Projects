function planning(ref, M)
%PLANNING  Figures explaining the Task 4 planning pipeline, step by step.
%
%   viz.planning(ref, M)        % ref, M = traj.mission(...)
%
%   Uses the internals traj.mission keeps in M.plan (tree, raw and shortcut
%   paths, waypoint lists, initial / optimised segment times, pre-dilation
%   profile). One figure per step:
%     1  real vs inflated footprints, route of the point robot
%     2  RRT*: growing tree / rejected edge / accepted route
%     3  raw RRT* path vs shortcut path, climb and descent legs
%     4  densification of a long segment, full waypoint list in altitude
%     5  angular polyline vs minimum-snap curve, acceleration of each
%     6  initial vs optimised segment times, speed not bounded by 3 m/s
%     7  uniform time dilation: speed curve rescaled, route unchanged
%     8  10-channel reference feeding the controller; unused quaternion

P  = M.plan;
C  = struct('real', [0.30 0.30 0.30], 'infl', [0.35 0.35 0.35], 'route', [0.00 0.45 0.74], ...
            'raw', [0.85 0.33 0.10], 'rej', [0.80 0.00 0.00], 'tree', [0.70 0.70 0.70], 'ms', [0.47 0.67 0.19]);
pr = traj.pos(ref); vr = traj.vel(ref); ar = traj.acc(ref); tr = ref.time(:).';
cruise = tr >= M.t_climb & tr <= M.t_descend;

%% 1 ---- real vs inflated footprint ----------------------------------------
figure('Name', 'Planning map: real vs inflated footprints', 'Color', 'w', 'Position', [100 100 900 420]);
ax = axes; drawMap(ax, M, C, true);
plot(ax, pr(1,:), pr(2,:), '-', 'Color', C.route, 'LineWidth', 2.2, 'DisplayName', 'route of the point robot');
plot(ax, P.dense(:,1), P.dense(:,2), '.', 'Color', C.route, 'MarkerSize', 12, 'HandleVisibility', 'off');
endpoints(ax, M);
legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal');
title(ax, sprintf('Buildings grown by the vehicle radius (%.0f m): the vehicle becomes a point', M.uav_radius));

%% 2 ---- RRT*: growing tree, rejected edge, accepted route ------------------
E = treeEdges(P.tree);                                  % one row per edge [x1 y1 x2 y2]
nE = size(E, 1); part = E(1:max(1, round(0.25*nE)), :);
[e0, e1] = rejectedEdge(part, M);
figure('Name', 'RRT star: growing tree, rejected edge, accepted route', 'Color', 'w', 'Position', [100 100 1500 420]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact');
ax = nexttile(tl); drawMap(ax, M, C, false); drawEdges(ax, part, C.route, 1.0);
plot(ax, M.p_climb(1), M.p_climb(2), 'ko', 'MarkerFaceColor', 'k'); title(ax, '(a) tree grows from the start');
ax = nexttile(tl); drawMap(ax, M, C, false); drawEdges(ax, part, C.tree, 0.8);
plot(ax, [e0(1) e1(1)], [e0(2) e1(2)], '--', 'Color', C.rej, 'LineWidth', 2.4);
plot(ax, e1(1), e1(2), 'x', 'Color', C.rej, 'MarkerSize', 14, 'LineWidth', 2.5);
plot(ax, e0(1), e0(2), 'o', 'Color', C.route, 'MarkerFaceColor', C.route);
title(ax, '(b) edge through an inflated footprint: rejected');
ax = nexttile(tl); drawMap(ax, M, C, false); drawEdges(ax, E, [0.85 0.85 0.85], 0.5);
plot(ax, P.raw(:,1), P.raw(:,2), '-', 'Color', C.raw, 'LineWidth', 2.8);
endpoints(ax, M); title(ax, '(c) accepted route');

%% 3 ---- raw RRT* vs shortcut, with the vertical legs ----------------------
figure('Name', 'Raw RRT star path vs shortcut path', 'Color', 'w', 'Position', [100 100 1000 560]);
ax = axes; hold(ax, 'on'); grid(ax, 'on'); drawBoxes(ax, M);
plot3(ax, P.raw(:,1), P.raw(:,2), P.raw(:,3), '.-', 'Color', C.raw, 'LineWidth', 1.4, 'MarkerSize', 10, 'DisplayName', 'raw RRT* path');
plot3(ax, P.shortcut(:,1), P.shortcut(:,2), P.shortcut(:,3), 'o-', 'Color', C.route, 'LineWidth', 2.6, ...
      'MarkerFaceColor', C.route, 'DisplayName', 'shortcut path (cruise route)');
plot3(ax, M.p_takeoff(1)*[1 1], M.p_takeoff(2)*[1 1], [0 M.cruise_alt], 'k--', 'LineWidth', 2.2, 'DisplayName', 'climb / descent legs');
plot3(ax, M.p_land(1)*[1 1], M.p_land(2)*[1 1], [M.cruise_alt 0], 'k--', 'LineWidth', 2.2, 'HandleVisibility', 'off');
plot3(ax, M.p_takeoff(1), M.p_takeoff(2), 0, 'k^', 'MarkerSize', 11, 'MarkerFaceColor', [0.2 0.7 0.2], 'DisplayName', 'start (ground)');
plot3(ax, M.p_land(1), M.p_land(2), 0, 'kv', 'MarkerSize', 11, 'MarkerFaceColor', [0.85 0.2 0.2], 'DisplayName', 'destination (ground)');
text(ax, M.p_takeoff(1), M.p_takeoff(2), M.cruise_alt/2, '  climb', 'FontSize', 12);
text(ax, M.p_land(1), M.p_land(2), M.cruise_alt/2, '  descent', 'FontSize', 12);
axis(ax, 'equal'); xlim(ax, [0 M.map_size(1)]); ylim(ax, [0 M.map_size(2)]); zlim(ax, [0 M.map_size(3)]);
view(ax, [-25 35]); xlabel(ax, 'x  [m]'); ylabel(ax, 'y  [m]'); zlabel(ax, 'z  [m]');
legend(ax, 'Location', 'northeastoutside');
title(ax, sprintf('Shortcutting: %d raw vertices %s %d', size(P.raw,1), char(8594), size(P.shortcut,1)));

%% 4 ---- densification and the full waypoint list --------------------------
seg = vecnorm(diff(P.shortcut(:,1:2)), 2, 2); [~, j] = max(seg);
a0 = P.shortcut(j,1:2); a1 = P.shortcut(j+1,1:2);
on = abs(cross2(P.dense(:,1:2) - a0, a1 - a0)) < 1e-6 & ...
     dot2(P.dense(:,1:2) - a0, a1 - a0) >= -1e-9 & dot2(P.dense(:,1:2) - a1, a0 - a1) >= -1e-9;
D = P.dense(on, 1:2); gap = max(vecnorm(diff(sortrows(D)), 2, 2));
figure('Name', 'Waypoint densification and assembled list', 'Color', 'w', 'Position', [100 100 1300 430]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); drawMap(ax, M, C, false);
plot(ax, [a0(1) a1(1)], [a0(2) a1(2)], '-', 'Color', C.route, 'LineWidth', 2.2);
plot(ax, D(:,1), D(:,2), '.', 'Color', C.raw, 'MarkerSize', 22);
plot(ax, [a0(1) a1(1)], [a0(2) a1(2)], 'o', 'Color', C.route, 'MarkerSize', 10, 'MarkerFaceColor', C.route);
pad = 12; xlim(ax, sort([a0(1) a1(1)]) + [-pad pad]); ylim(ax, sort([a0(2) a1(2)]) + [-pad pad]);
title(ax, sprintf('Longest segment %.1f m: points every %.1f m (%s %.0f m)', seg(j), gap, char(8804), P.max_spacing));
full = P.wp; if norm(full(end,:) - M.p_land) > 1e-6, full = [full; M.p_land]; end
s = [0; cumsum(vecnorm(diff(full), 2, 2))];
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
plot(ax, s, full(:,3), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
plot(ax, s, full(:,3), 'o', 'Color', C.route, 'MarkerFaceColor', C.route, 'MarkerSize', 5);
plot(ax, s([1 end]), full([1 end],3), 's', 'Color', 'k', 'MarkerFaceColor', [0.85 0.2 0.2], 'MarkerSize', 10);
yline(ax, 0, 'k-', 'ground', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'center');
xlabel(ax, 'distance along the waypoint list  [m]'); ylabel(ax, 'altitude  [m]'); ylim(ax, [-3 M.cruise_alt + 5]);
title(ax, sprintf('Full list: %d waypoints (first at the hover point, last on the pad)', size(full,1)));

%% 5 ---- angular polyline vs smooth curve ----------------------------------
wpc = P.wp(P.wp(:,3) > M.cruise_alt - 1e-6, :);          % cruise waypoints
tc  = tr(cruise) - M.t_climb; pc = pr(:, cruise); ac = vecnorm(ar(:, cruise), 2, 1);
[tp, app, vmax] = polylineAccel(wpc, tc(end));
figure('Name', 'Angular polyline vs minimum snap curve', 'Color', 'w', 'Position', [100 100 1300 430]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); drawMap(ax, M, C, false);
plot(ax, wpc(:,1), wpc(:,2), 'o-', 'Color', C.raw, 'LineWidth', 1.6, 'MarkerSize', 5, 'MarkerFaceColor', C.raw, 'DisplayName', 'waypoint polyline');
plot(ax, pc(1,:), pc(2,:), '-', 'Color', C.route, 'LineWidth', 2.4, 'DisplayName', 'minimum-snap curve');
legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal'); title(ax, 'Same waypoints');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
top = 3*max(ac); app(app < 0.02*top) = 0;
stem(ax, tp(app > 0), min(app(app > 0), top), 'Color', C.raw, 'LineWidth', 1.6, 'Marker', '^', 'MarkerFaceColor', C.raw, 'DisplayName', 'polyline: impulses at corners');
plot(ax, tc, ac, '-', 'Color', C.route, 'LineWidth', 2.2, 'DisplayName', 'minimum snap: continuous');
ylim(ax, [0 top*1.05]); xlim(ax, [0 tc(end)]);
xlabel(ax, 'cruise time  [s]'); ylabel(ax, '|a|  [m/s^2]'); legend(ax, 'Location', 'northeast');
title(ax, sprintf('Acceleration (polyline flown at %.1f m/s; spikes clipped)', vmax));

%% 6 ---- initial vs optimised segment times --------------------------------
figure('Name', 'Time allocation: initial vs optimised times', 'Color', 'w', 'Position', [100 100 1300 460]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on');
ti = P.t_wp_init(:).'; to = P.t_wp_opt(:).';
for i = 1:numel(ti), plot(ax, [ti(i) to(i)], [1 0], '-', 'Color', [0.8 0.8 0.8]); end
plot(ax, ti, ones(size(ti)), '|', 'Color', C.raw, 'MarkerSize', 16, 'LineWidth', 2);
plot(ax, to, zeros(size(to)), '|', 'Color', C.route, 'MarkerSize', 16, 'LineWidth', 2);
plot(ax, [0 max([ti to])], [1 1], '-', 'Color', C.raw); plot(ax, [0 max([ti to])], [0 0], '-', 'Color', C.route);
yticks(ax, [0 1]); yticklabels(ax, {sprintf('optimised (%.1f s)', to(end)), sprintf('initial: length / %.0f m/s (%.1f s)', P.speed, ti(end))});
ylim(ax, [-0.5 1.5]); xlabel(ax, 'waypoint arrival time  [s]'); title(ax, 'Same route, two timelines');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
sp = vecnorm(P.pre.v, 2, 1);
plot(ax, P.pre.t, sp, '-', 'Color', C.route, 'LineWidth', 2.2);
yline(ax, P.speed, '--', sprintf('%.0f m/s: initial guess only, not a bound', P.speed), 'Color', C.raw, 'LineWidth', 1.6);
[pk, ip] = max(sp); plot(ax, P.pre.t(ip), pk, 'v', 'Color', 'k', 'MarkerFaceColor', C.rej, 'MarkerSize', 10);
text(ax, P.pre.t(ip), pk, sprintf('  peak %.1f m/s', pk), 'VerticalAlignment', 'bottom');
xlabel(ax, 't  [s]'); ylabel(ax, 'speed  [m/s]'); xlim(ax, [0 P.pre.t(end)]);

%% 7 ---- uniform time dilation ---------------------------------------------
k = P.k;
figure('Name', 'Time dilation: speed rescaled, route unchanged', 'Color', 'w', 'Position', [100 100 1400 430]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact');
ax = nexttile(tl, [1 2]); hold(ax, 'on'); grid(ax, 'on');
plot(ax, P.pre.t, sp, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 2, 'DisplayName', 'optimised');
plot(ax, P.pre.t*k, sp/k, '-', 'Color', C.route, 'LineWidth', 2.4, 'DisplayName', sprintf('dilated: t %s %.2f t,  v %s v / %.2f', char(8594), k, char(8594), k));
yline(ax, M.max_speed, '--', sprintf('cap %.1f m/s', M.max_speed), 'Color', C.rej, 'LineWidth', 1.6, 'HandleVisibility', 'off');
xlabel(ax, 't  [s]'); ylabel(ax, 'speed  [m/s]'); xlim(ax, [0 P.pre.t(end)*k]); legend(ax, 'Location', 'northwest');
title(ax, 'Time stretched, speed reduced by the same factor');
in = nexttile(tl); drawMap(in, M, C, false);
plot(in, P.pre.p(1,:), P.pre.p(2,:), '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 3);
plot(in, pr(1, tr <= M.t_descend), pr(2, tr <= M.t_descend), '--', 'Color', C.route, 'LineWidth', 1.6);
set(in, 'XTick', [], 'YTick', []); title(in, 'route: unchanged', 'FontWeight', 'normal');

%% 7b ---- landing: separate rest-to-rest descent ---------------------------
if isfield(M, 'descent')
    Ld = tr >= M.t_descend; tt = tr(Ld) - M.t_descend;
    figure('Name', 'Landing descent planned from an acceleration budget', 'Color', 'w', 'Position', [100 100 1300 400]);
    tl = tiledlayout(1, 3, 'TileSpacing', 'compact');
    ax = nexttile(tl); plot(ax, tt, pr(3, Ld), 'Color', C.route, 'LineWidth', 2.2); grid(ax, 'on');
    xlabel(ax, 't - t_{descend}  [s]'); ylabel(ax, 'z  [m]'); title(ax, sprintf('%.0f m in %.1f s', M.descent.height, M.descent.T));
    ax = nexttile(tl); plot(ax, tt, vr(3, Ld), 'Color', C.route, 'LineWidth', 2.2); grid(ax, 'on');
    xlabel(ax, 't - t_{descend}  [s]'); ylabel(ax, 'v_z  [m/s]'); title(ax, sprintf('peak %.2f m/s', M.descent.v_peak));
    ax = nexttile(tl); plot(ax, tt, ar(3, Ld), 'Color', C.route, 'LineWidth', 2.2); grid(ax, 'on'); hold(ax, 'on');
    yline(ax, [-1 1]*M.land_acc, '--', 'Color', C.rej, 'LineWidth', 1.4);
    xlabel(ax, 't - t_{descend}  [s]'); ylabel(ax, 'a_z  [m/s^2]'); title(ax, sprintf('|a_z| %s %.1f m/s^2', char(8804), M.land_acc));
end

%% 8 ---- reference interface -------------------------------------------------
referenceDiagram(C);
end

% =========================================================================
function drawMap(ax, M, C, withLegend)
hold(ax, 'on'); axis(ax, 'equal'); box(ax, 'on');
r = M.uav_radius;
for b = 1:size(M.buildings, 1)
    B = M.buildings(b, :);
    h1 = rectangle(ax, 'Position', [B(1) B(3) B(2)-B(1) B(4)-B(3)], 'FaceColor', C.real, 'EdgeColor', C.real);
    h2 = rectangle(ax, 'Position', [B(1)-r B(3)-r B(2)-B(1)+2*r B(4)-B(3)+2*r], 'EdgeColor', C.infl, ...
                   'LineStyle', '--', 'LineWidth', 1.3);
    if b == 1 && withLegend      % rectangles have no legend entry: add proxies
        plot(ax, NaN, NaN, 's', 'MarkerFaceColor', C.real, 'MarkerEdgeColor', C.real, 'MarkerSize', 12, 'DisplayName', 'real footprint');
        plot(ax, NaN, NaN, '--', 'Color', C.infl, 'LineWidth', 1.3, 'DisplayName', sprintf('inflated by %.0f m', r));
    end
    %#ok<*NASGU>
end
xlim(ax, [0 M.map_size(1)]); ylim(ax, [0 M.map_size(2)]);
xlabel(ax, 'x  [m]'); ylabel(ax, 'y  [m]');
end

function drawBoxes(ax, M)
F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
for b = 1:size(M.buildings, 1)
    B = M.buildings(b, :);
    V = [B(1) B(3) B(5); B(2) B(3) B(5); B(2) B(4) B(5); B(1) B(4) B(5);
         B(1) B(3) B(6); B(2) B(3) B(6); B(2) B(4) B(6); B(1) B(4) B(6)];
    patch(ax, 'Vertices', V, 'Faces', F, 'FaceColor', [0.45 0.45 0.45], 'FaceAlpha', 0.35, ...
          'EdgeColor', [0.3 0.3 0.3], 'HandleVisibility', 'off');
end
end

function endpoints(ax, M)
plot(ax, M.p_takeoff(1), M.p_takeoff(2), 'k^', 'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.7 0.2], 'DisplayName', 'start');
plot(ax, M.p_land(1), M.p_land(2), 'kv', 'MarkerSize', 10, 'MarkerFaceColor', [0.85 0.2 0.2], 'DisplayName', 'destination');
end

function E = treeEdges(T)
% TreeData lists edges as [parent; child; NaN] triples in insertion order.
E = zeros(0, 4); i = 1;
while i < size(T, 1)
    if all(isfinite(T(i,:))) && all(isfinite(T(i+1,:)))
        E(end+1, :) = [T(i,1:2) T(i+1,1:2)]; i = i + 3; %#ok<AGROW>
    else
        i = i + 1;
    end
end
end

function drawEdges(ax, E, col, lw)
X = [E(:,1) E(:,3) nan(size(E,1),1)].'; Y = [E(:,2) E(:,4) nan(size(E,1),1)].';
plot(ax, X(:), Y(:), '-', 'Color', col, 'LineWidth', lw);
end

function [e0, e1] = rejectedEdge(E, M)
% A tree node near a building, extended 15 m towards the building's centre:
% the motion check finds the inflated footprint on the way and rejects it.
nodes = unique([E(:,1:2); E(:,3:4)], 'rows');
best = inf; e0 = nodes(1,:); e1 = e0;
r = M.uav_radius;
for b = 1:size(M.buildings, 1)
    B = M.buildings(b, :); c = [mean(B(1:2)) mean(B(3:4))];
    for n = 1:size(nodes, 1)
        q = nodes(n,:); u = (c - q) / norm(c - q); t = q + 15*u;
        if inBox(q, B, r), continue, end
        crosses = any(arrayfun(@(s) inBox(q + s*(t - q), B, r), linspace(0, 1, 60)));
        if crosses && norm(c - q) < best, best = norm(c - q); e0 = q; e1 = t; end
    end
end
end

function tf = inBox(q, B, r)
tf = q(1) >= B(1)-r && q(1) <= B(2)+r && q(2) >= B(3)-r && q(2) <= B(4)+r;
end

function [t, amag, v] = polylineAccel(W, T)
% The waypoint polyline flown at constant speed in the same total time: the
% velocity jumps at every corner, i.e. an impulse of acceleration.
L = vecnorm(diff(W(:,1:2)), 2, 2); v = sum(L) / T; dt = 0.01;
tk = [0; cumsum(L)/v];
t = 0:dt:T; amag = zeros(size(t));
for i = 2:numel(tk)-1
    d1 = (W(i,1:2) - W(i-1,1:2)) / L(i-1); d2 = (W(i+1,1:2) - W(i,1:2)) / L(i);
    [~, j] = min(abs(t - tk(i))); amag(j) = norm(v*(d2 - d1)) / dt;
end
end

function z = cross2(a, b), z = a(:,1).*b(2) - a(:,2).*b(1); end
function z = dot2(a, b),   z = a(:,1).*b(1) + a(:,2).*b(2); end

function referenceDiagram(C)
figure('Name', 'Reference interface: 10 channels into the controller', 'Color', 'w', 'Position', [100 100 1400 560]);
ax = axes('Position', [0 0 1 1]); axis(ax, [0 140 0 56]); axis(ax, 'off'); hold(ax, 'on');
blk = @(x, y, w, h, fc, ec) rectangle(ax, 'Position', [x y w h], 'Curvature', 0.08, 'FaceColor', fc, 'EdgeColor', ec, 'LineWidth', 1.6);
arw = @(x0, y0, x1, y1, col) annotation('arrow', [x0 x1]/140, [y0 y1]/56, 'Color', col, 'LineWidth', 1.6, 'HeadLength', 9, 'HeadWidth', 9);

% RRT* state: position used, quaternion unused
text(ax, 13, 51, 'RRT* state (SE(3))', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 13);
lab = {'x','y','z','q_w','q_x','q_y','q_z'};
for i = 1:7
    used = i <= 3;
    blk(2 + (i-1)*3.3, 42, 3.1, 5, ternary(used, [0.86 0.92 0.98], [0.93 0.93 0.93]), ternary(used, C.route, [0.7 0.7 0.7]));
    text(ax, 3.55 + (i-1)*3.3, 44.5, lab{i}, 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', ternary(used, [0 0 0], [0.6 0.6 0.6]));
end
plot(ax, [11.9 25.1], [42 47], '-', 'Color', C.rej, 'LineWidth', 1.8);
text(ax, 18.5, 39.5, 'orientation: identity, never used', 'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', C.rej);
arw(8, 41.5, 8, 33, C.route);

% trajectory generation
blk(1, 20, 26, 12.5, [0.96 0.96 0.96], [0.4 0.4 0.4]);
text(ax, 14, 29.5, 'waypoints  \rightarrow  minimum snap', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(ax, 14, 26, 'p(t), v(t), a(t)', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(ax, 14, 22.3, '\psi(t) = atan2(v_y, v_x)', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', C.raw, 'FontWeight', 'bold');
arw(27, 26, 34, 26, [0.3 0.3 0.3]);

% 10-channel reference
text(ax, 58, 16, 'traj\_ref: 10-channel From Workspace signal', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 13);
ch = {'p_x','p_y','p_z','v_x','v_y','v_z','a_x','a_y','a_z','\psi'};
for i = 1:10
    isyaw = i == 10;
    blk(35 + (i-1)*4.6, 23, 4.4, 6, ternary(isyaw, [0.99 0.90 0.80], [0.86 0.92 0.98]), ternary(isyaw, C.raw, C.route));
    text(ax, 37.2 + (i-1)*4.6, 26, ch{i}, 'HorizontalAlignment', 'center', 'FontSize', 12);
end
text(ax, 55.7, 20.5, 'p_d, v_d, a_d', 'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', C.route);
text(ax, 78.7, 20.5, '\psi_d', 'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', C.raw);

% controller (Task 3)
blk(92, 36, 22, 10, [0.86 0.92 0.98], C.route); text(ax, 103, 42.5, 'Position controller', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(ax, 103, 39, '+ PBMRAC', 'HorizontalAlignment', 'center', 'FontSize', 12);
blk(92, 21, 22, 9, [0.99 0.90 0.80], C.raw);  text(ax, 103, 25.5, 'Attitude planner', 'HorizontalAlignment', 'center', 'FontSize', 12);
blk(118, 21, 20, 9, [0.96 0.96 0.96], [0.4 0.4 0.4]); text(ax, 128, 25.5, 'Attitude ctrl', 'HorizontalAlignment', 'center', 'FontSize', 12);
blk(118, 6, 20, 9, [0.96 0.96 0.96], [0.4 0.4 0.4]);  text(ax, 128, 10.5, 'Allocation', 'HorizontalAlignment', 'center', 'FontSize', 12);
arw(58, 29.5, 92, 41, C.route);
arw(80, 23, 92, 25.5, C.raw);
arw(103, 36, 103, 30, [0.3 0.3 0.3]); text(ax, 104, 33, 'f_c', 'FontSize', 11);
arw(114, 25.5, 118, 25.5, [0.3 0.3 0.3]);
arw(128, 21, 128, 15, [0.3 0.3 0.3]);
text(ax, 125, 50, 'Task 3 controller', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 13);
end

function out = ternary(c, a, b), if c, out = a; else, out = b; end, end
