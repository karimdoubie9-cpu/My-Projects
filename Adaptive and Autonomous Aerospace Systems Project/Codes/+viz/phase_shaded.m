function fig = phase_shaded(ref, phases, opts)
%PHASE_SHADED  Planned x, y, z with the mission phases shaded behind them.
%
%   viz.phase_shaded(ref, M.phases)
%   viz.phase_shaded(ref, [t1 t2], names={'Ascent','Cruise','Descent'})
%
%   `phases` holds the interior boundary times, so N boundaries give N+1 bands
%   spanning [0, T]. Shading the phase behind the curve rather than marking the
%   boundaries with vertical lines makes it obvious which segment of x(t) belongs
%   to the climb and which to the cruise -- with three stacked axes, vertical
%   lines alone force the reader to track the boundary across all three.

arguments
    ref         struct
    phases      double
    opts.names  cell = {}
    opts.name   (1,:) char = 'Planned trajectory with phase shading'
end

t = ref.time;
p = traj.pos(ref);

edges = [t(1), reshape(phases,1,[]), t(end)];
nph   = numel(edges) - 1;

names = opts.names;
if isempty(names)
    default = {'Ascent','Cruise','Descent','Phase 4','Phase 5'};
    names = default(1:min(nph, numel(default)));
    while numel(names) < nph, names{end+1} = sprintf('Phase %d', numel(names)+1); end %#ok<AGROW>
end

shades = {[0.85 0.95 0.85], [0.85 0.90 1.00], [1.00 0.90 0.85], ...
          [0.95 0.90 1.00], [1.00 0.98 0.85]};

fig = figure('Name', opts.name, 'Color','w', 'Position',[120 120 850 680]);
labels = {'x  [m]','y  [m]','z  [m]'};

for ax = 1:3
    subplot(3,1,ax); hold on
    yl = [min(p(ax,:)) max(p(ax,:))];
    if diff(yl) < 1e-3, yl = yl + [-1 1]; end
    yl = yl + 0.05*diff(yl)*[-1 1];

    for ph = 1:nph
        patch([edges(ph) edges(ph+1) edges(ph+1) edges(ph)], ...
              [yl(1) yl(1) yl(2) yl(2)], shades{1+mod(ph-1,numel(shades))}, ...
              'EdgeColor','none', 'HandleVisibility','off');
    end
    plot(t, p(ax,:), 'k', 'LineWidth', 2.3, 'HandleVisibility','off');
    grid on; ylabel(labels{ax}); ylim(yl); xlim([t(1) t(end)]);

    if ax == 1
        % Zero-area patches purely to populate the legend with the phase names.
        for ph = 1:nph
            patch(NaN, NaN, shades{1+mod(ph-1,numel(shades))}, 'DisplayName', names{ph});
        end
        legend('Location','eastoutside');
        title(opts.name);
    end
end
xlabel('t  [s]');
end
