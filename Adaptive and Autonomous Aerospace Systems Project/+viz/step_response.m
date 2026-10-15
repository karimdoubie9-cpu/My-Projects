function [fig, S] = step_response(R, ax, opts)
%STEP_RESPONSE  Annotated setpoint response with the design objectives drawn on.
%
%   [fig, S] = viz.step_response(R, 1)
%   [fig, S] = viz.step_response(R, 1, os_target=35, ts_target=15)
%
%   Plots the commanded and actual position on the chosen axis and marks every
%   quantity the objectives are stated in, so the figure can be read without a
%   separate table:
%
%     * the +-2% settling band around the command (shaded)
%     * the peak, with the overshoot percentage
%     * the settling time, and the target if one is given
%     * the 10-90% rise window
%
%   Drawing the BAND rather than just quoting a number is the point: settling
%   time is meaningless without saying what tolerance it settles to, and a
%   reader can see immediately whether the response is creeping across the
%   boundary or genuinely converged.

arguments
    R              struct
    ax             (1,1) double {mustBeMember(ax,[1 2 3])} = 1
    opts.threshold (1,1) double = 0.02
    opts.os_target (1,1) double = NaN
    opts.ts_target (1,1) double = NaN
    opts.label     (1,:) char = ''
    opts.name      (1,:) char = 'Step response'
    opts.quiet     (1,1) logical = false
end

S = harness.step_metrics(R, ax, threshold = opts.threshold, ...
                         os_target = opts.os_target, ts_target = opts.ts_target, ...
                         label = opts.label, quiet = opts.quiet);

s      = viz.style();
axname = 'xyz';
t      = R.t;
y      = R.p(:,ax);
yd     = R.p_d(:,ax);

fig = figure('Name', opts.name, 'Color','w', 'Position',[100 100 860 520]);
hold on; grid on

% Settling band around the commanded value
band = abs(opts.threshold * S.step);
xl   = [t(1) t(end)];
patch([xl(1) xl(2) xl(2) xl(1)], ...
      [S.to-band S.to-band S.to+band S.to+band], [0.4 0.7 0.4], ...
      'FaceAlpha', 0.15, 'EdgeColor','none', ...
      'DisplayName', sprintf('\\pm%.0f%% band', 100*opts.threshold));

plot(t, yd, s.ref_style, 'Color', s.ref_color, 'LineWidth', s.lw_ref, ...
     'DisplayName', 'Command');
plot(t, y, 'LineWidth', s.lw, 'Color', s.colors{1}, 'DisplayName', 'Response');

% Peak
if isfinite(S.PeakTime)
    plot(S.t_step + S.PeakTime, S.PeakValue, 'v', 'MarkerSize', 9, ...
         'MarkerFaceColor', s.colors{2}, 'MarkerEdgeColor','k', ...
         'DisplayName', sprintf('peak: %.1f%% OS', S.Overshoot));
end

% Settling time, achieved and targeted
if isfinite(S.SettlingTime)
    xline(S.t_step + S.SettlingTime, '-', sprintf('T_s = %.2f s', S.SettlingTime), ...
          'Color', s.colors{2}, 'LineWidth', 2.1, ...
          'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
end
if ~isnan(opts.ts_target)
    xline(S.t_step + opts.ts_target, ':', sprintf('target %.0f s', opts.ts_target), ...
          'Color', [0 0.5 0], 'LineWidth', 1.8, 'HandleVisibility','off');
end

% Overshoot bound, if an objective was given
if ~isnan(opts.os_target)
    yline(S.to + opts.os_target/100*S.step, ':', ...
          sprintf('OS bound %.0f%%', opts.os_target), 'Color',[0.8 0 0], ...
          'LineWidth', 1.8, 'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
end

xlabel('t  [s]'); ylabel(sprintf('%s  [m]', axname(ax)));
title(sprintf('%s: T_s = %.2f s, OS = %.1f%%, e_{ss} = %.4f m', ...
              ternary(isempty(opts.label), sprintf('%s-axis step', axname(ax)), opts.label), ...
              S.SettlingTime, S.Overshoot, S.SSError));
% Legend at the right-middle and the metrics box anchored in the lower-right
% corner: after a step the response sits at the top of the axes, so both
% regions are empty. Placing them explicitly (not 'best') keeps them apart
% whatever the font size - 'best' is decided before the box exists.
legend('Location','east');

% Metrics box, so the figure stands alone on a slide. Kept monospaced: the
% values are aligned in columns with spaces, which a proportional font breaks.
txt = sprintf(['rise (10-90%%)  %.2f s\n' ...
               'peak           %.3f m\n' ...
               'overshoot      %.1f %%\n' ...
               'settling       %.2f s\n' ...
               'steady error   %.4f m'], ...
               S.RiseTime, S.PeakValue, S.Overshoot, S.SettlingTime, S.SSError);
% Left-aligned text keeps the columns straight; the box is then slid left by
% its own measured width so its right edge sits just inside the axes.
h = text(0.98, 0.04, txt, 'Units','normalized', 'FontName','FixedWidth', ...
         'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
         'BackgroundColor','w', 'EdgeColor',[0.4 0.4 0.4], 'FontSize', 10);
h.Position(1) = 0.98 - h.Extent(3);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
