function fig = histograms(MC, opts)
%HISTOGRAMS  Monte Carlo metric distributions, one panel per metric.
%
%   viz.histograms(MC)
%   viz.histograms(MC, metrics={'rms_p','peak_p'})
%
%   Overlays the configurations on shared bin edges -- comparing histograms
%   drawn with independently chosen bins is meaningless, and MATLAB's automatic
%   binning would pick different edges for each series.

arguments
    MC           struct
    opts.metrics cell = {'rms_p','peak_p','ss_p'}
    opts.name    (1,:) char = 'Monte Carlo distributions'
end

s  = viz.style();
nm = numel(opts.metrics);
labels = struct('rms_p','RMS |e_p|  [m]', 'peak_p','peak |e_p|  [m]', ...
                'ss_p','steady-state |e_p|  [m]', 'rms_v','RMS |e_v|  [m/s]');

fig = figure('Name', opts.name, 'Color','w', 'Position',[100 100 420*nm 420]);

for i = 1:nm
    f = opts.metrics{i};
    X = MC.(f);
    subplot(1, nm, i); hold on; grid on

    all = X(~isnan(X));
    edges = linspace(min(all), max(all), 30);

    for c = 1:size(X,2)
        x = X(~isnan(X(:,c)), c);
        histogram(x, edges, 'FaceColor', s.colors{1+mod(c-1,numel(s.colors))}, ...
                  'FaceAlpha', 0.55, 'EdgeColor','none', 'DisplayName', MC.tags{c});
        xline(median(x), '--', 'Color', s.colors{1+mod(c-1,numel(s.colors))}, ...
              'LineWidth', 2.1, 'HandleVisibility','off');
    end
    xlabel(labels.(f)); ylabel('count');
    title(sprintf('%s (dashed = median)', f), 'Interpreter','none');
    if i == 1, legend('Location','best'); end
end
end
