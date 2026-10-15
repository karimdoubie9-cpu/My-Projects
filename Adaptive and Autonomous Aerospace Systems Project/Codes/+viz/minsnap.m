function fig = minsnap(ref, opts)
%MINSNAP  Minimum-snap reference profiles: x(t), y(t), |v|(t), |a|(t).
%
%   viz.minsnap(ref)
%   viz.minsnap(ref, phases=M.phases, limit=M.max_speed)
%
%   Four stacked panels. x and y are shown as signed coordinates because the
%   horizontal route is what the planner solved for; speed and acceleration are
%   shown as magnitudes because those are the quantities that carry limits --
%   the speed cap from the exam, and the acceleration the attitude loop has to
%   deliver. Per-axis velocity would just repeat viz.kinematics.

arguments
    ref         struct
    opts.phases double = []
    opts.limit  (1,1) double = NaN
    opts.name   (1,:) char = 'Minimum-snap reference profiles'
end

t = ref.time;
p = traj.pos(ref);
v = traj.vel(ref);
a = traj.acc(ref);

fig = figure('Name', opts.name, 'Color','w', 'Position',[140 140 850 740]);

series = {p(1,:), p(2,:), vecnorm(v,2,1), vecnorm(a,2,1)};
labels = {'x  [m]','y  [m]','|v|  [m/s]','|a|  [m/s^2]'};

for i = 1:4
    subplot(4,1,i); hold on; grid on
    plot(t, series{i}, 'LineWidth', 2.1, 'HandleVisibility','off');
    if i == 3 && ~isnan(opts.limit)
        yline(opts.limit, '--', sprintf('%.1f m/s limit', opts.limit), ...
              'Color',[0.8 0 0], 'LineWidth',1.7, ...
              'LabelHorizontalAlignment','left', 'HandleVisibility','off');
    end
    for k = 1:numel(opts.phases)
        xline(opts.phases(k), '--k', 'HandleVisibility','off');
    end
    ylabel(labels{i}); xlim([t(1) t(end)]);
    if i == 1, title(opts.name); end
    if i == 4, xlabel('t  [s]'); end
end
end
