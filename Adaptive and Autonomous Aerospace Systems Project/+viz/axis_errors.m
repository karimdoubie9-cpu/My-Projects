function fig = axis_errors(runs, tags, opts)
%AXIS_ERRORS  Per-axis tracking errors on log scale, with a 3D companion panel.
%
%   viz.axis_errors({R}, {'Baseline'})
%   viz.axis_errors({R0,R1}, {'No package','With package'}, show3d=true)
%
%   Three panels, the layout used throughout the deck:
%     1. |eps_p| per axis (x, y, z), semilog-y
%     2. |eps_v| per axis, semilog-y
%     3. the 3D trajectory, reference dashed plus one line per run
%
%   NOTATION: per-axis components are eps (epsilon), reserved here; the 2-norms
%   are e and live in viz.errors. Keeping the two apart matters because a large
%   |eps_z| with small |eps_x|,|eps_y| is a very different failure from an
%   isotropic error of the same norm.
%
%   COMPARISON STYLE: with two runs the first is drawn dotted and the second
%   solid, sharing the axis colour, so x/y/z stay identifiable across both
%   cases -- six solid lines in three colours would be unreadable. The colour
%   order is reset between the two so axis k has the same colour in both.

arguments
    runs         cell
    tags         cell
    opts.show3d  (1,1) logical = true
    opts.n       (1,1) double = 6000
    opts.events  double = zeros(0,2)
    opts.name    (1,:) char = 'Per-axis tracking error'
end

nr = numel(runs);
M  = cell(1,nr);
for k = 1:nr
    M{k} = harness.metrics(runs{k}, n = opts.n);
end

nPanels = 2 + double(opts.show3d);
fig = figure('Name', opts.name, 'Color','w', ...
             'Position',[80 80 400*nPanels 430]);

styles = {':','-','--','-.'};
labels = {{'\epsilon_x','\epsilon_y','\epsilon_z'}, ...
          {'\epsilon_{v_x}','\epsilon_{v_y}','\epsilon_{v_z}'}};
ylabs  = {'|\epsilon_p|  [m]', '|\epsilon_v|  [m/s]'};
titles = {'Position error per axis', 'Velocity error per axis'};

for panel = 1:2
    subplot(1, nPanels, panel); hold on; grid on
    for k = 1:nr
        if panel == 1, E = M{k}.ep; else, E = M{k}.ev; end
        set(gca, 'ColorOrderIndex', 1);        % same colour per axis in every run
        semilogy(M{k}.t, abs(E), styles{1+mod(k-1,numel(styles))}, 'LineWidth', 1.8);
    end
    set(gca, 'YScale', 'log');
    markEvents(opts.events);
    xlabel('t  [s]'); ylabel(ylabs{panel});
    if nr > 1
        title(sprintf('%s\n(%s: %s,  %s: %s)', titles{panel}, ...
              styles{1}, tags{1}, styles{2}, tags{2}));
    else
        title(titles{panel});
    end
    legend(labels{panel}, 'Location','best');
end

if opts.show3d
    subplot(1, nPanels, 3); hold on; grid on
    s = viz.style();
    plot3(runs{1}.p_d(:,1), runs{1}.p_d(:,2), runs{1}.p_d(:,3), ...
          s.ref_style, 'Color', s.ref_color, 'LineWidth', s.lw_ref, ...
          'DisplayName','Reference');
    for k = 1:nr
        plot3(runs{k}.p(:,1), runs{k}.p(:,2), runs{k}.p(:,3), ...
              styles{1+mod(k-1,numel(styles))}, 'LineWidth', 2, ...
              'Color', s.colors{1+mod(k-1,numel(s.colors))}, 'DisplayName', tags{k});
    end
    axis equal; view(40,25);
    xlabel('x  [m]'); ylabel('y  [m]'); zlabel('z  [m]');
    title('3D trajectory'); legend('Location','best');
end
end

function markEvents(events)
for e = 1:size(events,1)
    xline(events(e,1), '--r', 'HandleVisibility','off');
end
end
