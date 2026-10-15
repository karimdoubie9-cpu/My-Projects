function fig = tracking(runs, tags, opts)
%TRACKING  Per-axis position (and optionally velocity) against the reference.
%
%   viz.tracking({R1,R2}, {'Baseline','PBMRAC'})
%   viz.tracking(runs, tags, show='both', events=[60 80])
%   viz.tracking(runs, tags, overlay=wind, overlayLabel='F_{gust}  [N]')
%
%   The reference is drawn once, dashed and black; the flown traces overlay it.
%   Plotting the reference per-run draws it two or three times on top of itself,
%   which makes the dashed line look solid.
%
%   OVERLAY: pass a From-Workspace style signal (fields time and
%   signals.values) to draw the disturbance on a right-hand y-axis behind the
%   tracking traces. Sharing the time axis with the response is what lets a
%   reader attribute an excursion to the disturbance that caused it, rather
%   than inferring it from a shaded band. Channel k of the overlay is drawn on
%   axis k when the widths match, otherwise its norm is drawn on every axis.

arguments
    runs             cell
    tags             cell
    opts.show        (1,:) char {mustBeMember(opts.show,{'position','both'})} = 'position'
    opts.events      double = zeros(0,2)
    opts.overlay     struct = struct([])
    opts.overlayLabel (1,:) char = 'disturbance'
    opts.name        (1,:) char = 'Tracking'
end

s          = viz.style();
both       = strcmp(opts.show, 'both');
rows       = 3 * (1 + both);
hasOverlay = ~isempty(fieldnames(opts.overlay));

if hasOverlay
    ov_t = opts.overlay.time;
    ov_v = opts.overlay.signals.values;
end

fig = figure('Name', opts.name, 'Color', 'w', ...
             'Position', [100 100 900 190*rows]);

for ax = 1:3
    subplot(rows, 1, ax);
    if hasOverlay
        yyaxis right
        plot(ov_t, overlayChannel(ov_v, ax), '-', 'Color', [0 0.45 0.74 0.55], ...
             'LineWidth', 1.4, 'DisplayName', opts.overlayLabel);
        ylabel(opts.overlayLabel);
        yyaxis left
        set(gca, 'YColor', 'k');
    end
    hold on; grid on
    plot(runs{1}.t, runs{1}.p_d(:,ax), s.ref_style, 'Color', s.ref_color, ...
         'LineWidth', s.lw_ref, 'DisplayName', 'Reference');
    for k = 1:numel(runs)
        plot(runs{k}.t, runs{k}.p(:,ax), 'LineWidth', s.lw, ...
             'Color', s.colors{1+mod(k-1,numel(s.colors))}, 'DisplayName', tags{k});
    end
    markEvents(opts.events);
    ylabel(sprintf('%s  [m]', s.axis_labels{ax}));
    if ax == 1, legend('Location','best'); title('Position tracking'); end
    if ~both && ax == 3, xlabel('t  [s]'); end
end

if both
    for ax = 1:3
        subplot(rows, 1, 3+ax); hold on; grid on
        plot(runs{1}.t, runs{1}.v_d(:,ax), s.ref_style, 'Color', s.ref_color, ...
             'LineWidth', s.lw_ref, 'DisplayName', 'Reference');
        for k = 1:numel(runs)
            plot(runs{k}.t, runs{k}.v(:,ax), 'LineWidth', s.lw, ...
                 'Color', s.colors{1+mod(k-1,numel(s.colors))}, 'DisplayName', tags{k});
        end
        markEvents(opts.events);
        ylabel(sprintf('v_%s  [m/s]', s.axis_labels{ax}));
        if ax == 1, title('Velocity tracking'); end
        if ax == 3, xlabel('t  [s]'); end
    end
end
end

% ---------------------------------------------------------------------
function y = overlayChannel(V, ax)
%OVERLAYCHANNEL  Channel ax of the overlay, or its norm if there is no channel.
if size(V,2) >= ax
    y = V(:,ax);
else
    y = vecnorm(V, 2, 2);
end
end

function markEvents(events)
for e = 1:size(events,1)
    xline(events(e,1), '--', 'Color', [0.8 0.2 0.2], 'HandleVisibility','off');
end
end
