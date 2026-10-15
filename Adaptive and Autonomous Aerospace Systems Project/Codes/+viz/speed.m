function fig = speed(runs, tags, opts)
%SPEED  Flown speed profile against the planned one.
%
%   viz.speed({R_bl, R_ad}, {'Baseline','PBMRAC'})
%   viz.speed(runs, tags, limit=4.5, events=E, eventNames={'Gust','Fault'})
%
%   Two panels:
%     1. speed magnitude |v|, which is the quantity the exam constrains
%        (no larger than 4-5 m/s) -- the limit is drawn when `limit` is given,
%        so a violation is visible rather than buried in the numbers;
%     2. vertical velocity v_z, which separates the mission phases: a positive
%        plateau during the climb, ~0 through the cruise, negative on descent.
%        The horizontal speed alone cannot show that structure.
%
%   The planned profile is taken from the logged reference (v_d) of the first
%   run, so it is drawn once rather than once per controller.

arguments
    runs            cell
    tags            cell
    opts.limit      (1,1) double = NaN
    opts.events     double = zeros(0,2)
    opts.eventNames cell = {}
    opts.phases     double = []
    opts.name       (1,:) char = 'Speed profile'
end

s   = viz.style();
fig = figure('Name', opts.name, 'Color', 'w', 'Position', [110 110 900 620]);

spd = @(V) vecnorm(V, 2, 2);

%% ---- Panel 1: speed magnitude --------------------------------------
subplot(2,1,1); hold on; grid on

ref_spd = spd(runs{1}.v_d);
allmax  = max(ref_spd);
for k = 1:numel(runs), allmax = max(allmax, max(spd(runs{k}.v))); end
ymax = 1.15*allmax;

shadeEvents(opts.events, opts.eventNames, ymax, s);

plot(runs{1}.t, ref_spd, s.ref_style, 'Color', s.ref_color, ...
     'LineWidth', s.lw_ref, 'DisplayName', 'Planned');
for k = 1:numel(runs)
    plot(runs{k}.t, spd(runs{k}.v), 'LineWidth', s.lw, ...
         'Color', s.colors{1+mod(k-1,numel(s.colors))}, 'DisplayName', tags{k});
end

if ~isnan(opts.limit)
    yline(opts.limit, '--', sprintf('%.1f m/s limit', opts.limit), ...
          'Color', [0.8 0 0], 'LineWidth', 1.8, ...
          'LabelHorizontalAlignment','left', 'HandleVisibility','off');
    ymax = max(ymax, 1.15*opts.limit);
end
markPhases(opts.phases);
ylim([0 ymax]);
ylabel('|v|  [m/s]'); title('Speed magnitude'); legend('Location','best');

%% ---- Panel 2: vertical velocity ------------------------------------
subplot(2,1,2); hold on; grid on
plot(runs{1}.t, runs{1}.v_d(:,3), s.ref_style, 'Color', s.ref_color, ...
     'LineWidth', s.lw_ref, 'DisplayName', 'Planned');
for k = 1:numel(runs)
    plot(runs{k}.t, runs{k}.v(:,3), 'LineWidth', s.lw, ...
         'Color', s.colors{1+mod(k-1,numel(s.colors))}, 'DisplayName', tags{k});
end
yline(0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
markPhases(opts.phases);
xlabel('t  [s]'); ylabel('v_z  [m/s]');
title('Vertical velocity (climb / cruise / descent)');

%% ---- Numbers --------------------------------------------------------
fprintf('\n%-14s %12s %12s %12s\n', 'case', 'peak |v|', 'mean |v|', 'peak |v_z|');
fprintf('%s\n', repmat('-', 1, 54));
fprintf('%-14s %12.3f %12.3f %12.3f\n', 'Planned', max(ref_spd), ...
        mean(ref_spd), max(abs(runs{1}.v_d(:,3))));
for k = 1:numel(runs)
    v = spd(runs{k}.v);
    fprintf('%-14s %12.3f %12.3f %12.3f\n', tags{k}, max(v), mean(v), ...
            max(abs(runs{k}.v(:,3))));
    if ~isnan(opts.limit) && max(v) > opts.limit
        fprintf('    NOTE: %s exceeds the %.1f m/s limit by %.2f m/s\n', ...
                tags{k}, opts.limit, max(v) - opts.limit);
    end
end
end

% ---------------------------------------------------------------------
function shadeEvents(events, names, ymax, s)
for e = 1:size(events,1)
    if e <= numel(names), nm = names{e}; else, nm = sprintf('event %d', e); end
    patch([events(e,1) events(e,2) events(e,2) events(e,1)], ...
          [0 0 ymax ymax], s.colors{3+mod(e-1,2)}, ...
          'FaceAlpha', s.event_alpha, 'EdgeColor','none', 'DisplayName', nm);
end
end

function markPhases(phases)
for k = 1:numel(phases)
    xline(phases(k), '--k', 'HandleVisibility','off');
end
end
