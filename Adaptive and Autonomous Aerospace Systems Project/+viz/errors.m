function fig = errors(runs, tags, opts)
%ERRORS  Position and velocity error norms, log scale, with optional events.
%
%   viz.errors({R1,R2}, {'Baseline','PBMRAC'})
%   viz.errors(runs, tags, events=[60 80; 90 150], eventNames={'Gust','Fault'})
%
%   Log scale is the right default here: the interesting comparison spans two
%   or three decades, and on a linear axis the adaptive trace collapses onto
%   the x-axis and looks like it does nothing.
%
%   'events' shades [start stop] windows so the reader can attribute each error
%   excursion to the disturbance that caused it instead of guessing from the
%   time axis.

arguments
    runs             cell
    tags             cell
    opts.events      double = zeros(0,2)
    opts.eventNames  cell   = {}
    opts.name        (1,:) char = 'Error norms'
end

s = viz.style();
M = cellfun(@(R) harness.metrics(R), runs, UniformOutput=false);

fig = figure('Name', opts.name, 'Color', 'w', 'Position', [100 100 900 620]);

for row = 1:2
    subplot(2,1,row); hold on; grid on
    if row == 1, series = @(m) m.nep; ylab = '|e_p|  [m]';   ttl = 'Position error norm';
    else,        series = @(m) m.nev; ylab = '|e_v|  [m/s]'; ttl = 'Velocity error norm';
    end

    allvals = cell2mat(cellfun(@(m) series(m), M, UniformOutput=false));
    ylim_hi = 1.3*max(allvals(:));
    if ~(ylim_hi > 0), ylim_hi = 1; end        % degenerate: all errors zero
    pos     = allvals(allvals > 0);
    if isempty(pos), ylim_lo = ylim_hi*1e-4;
    else,            ylim_lo = max(min(pos), ylim_hi*1e-4);
    end

    shadeEvents(opts.events, opts.eventNames, ylim_lo, ylim_hi, s, row == 1);

    for k = 1:numel(M)
        semilogy(M{k}.t, series(M{k}), 'LineWidth', s.lw, ...
                 'Color', s.colors{1+mod(k-1,numel(s.colors))}, 'DisplayName', tags{k});
    end
    set(gca, 'YScale', 'log'); ylim([ylim_lo ylim_hi])
    ylabel(ylab); title(ttl)
    if row == 1, legend('Location','best'); else, xlabel('t  [s]'); end
end
end

% ---------------------------------------------------------------------
function shadeEvents(events, names, lo, hi, s, showLegend)
for e = 1:size(events,1)
    if e <= numel(names), nm = names{e}; else, nm = sprintf('event %d', e); end
    if showLegend, vis = 'on'; else, vis = 'off'; end
    patch([events(e,1) events(e,2) events(e,2) events(e,1)], [lo lo hi hi], ...
          s.colors{3+mod(e-1,2)}, 'FaceAlpha', s.event_alpha, 'EdgeColor','none', ...
          'DisplayName', nm, 'HandleVisibility', vis);
end
end
