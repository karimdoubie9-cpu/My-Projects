function [fig, S] = sensitivity(MC, opts)
%SENSITIVITY  Attribute the Monte Carlo spread to each uncertainty.
%
%   [fig, S] = viz.sensitivity(MC)
%   [fig, S] = viz.sensitivity(MC, metric='rms_p')
%
%   Left panel  : share of variance per factor (from the standardized
%                 regression coefficients) -- the "who is to blame" chart.
%   Right panels: the metric against each factor, with a least-squares trend.
%
%   Read the scatter panels with the left panel in hand. In a joint campaign
%   every point is confounded -- when mass is high, lambda is whatever it
%   happened to be -- so a weak-looking scatter does NOT mean a weak effect. The
%   variance shares come from a multiple regression that holds the other factors
%   fixed, which the scatter cannot do.

arguments
    MC          struct
    opts.metric (1,:) char = 'rms_p'
    opts.name   (1,:) char = 'Monte Carlo sensitivity'
end

S = harness.sensitivity(MC, metric = opts.metric, quiet = true);

lam = MC.samples.lambda;
Xall  = [MC.samples.mass(:), mean(lam,2), min(lam,[],2), max(lam,[],2)-min(lam,[],2)];
nmAll = {'mass [kg]','mean(\lambda)','min(\lambda)','spread(\lambda)'};
keep  = std(Xall,0,1) > 1e-12;
X     = Xall(:,keep);
nm    = nmAll(keep);

nf  = size(X,2);
sty = viz.style();
Y   = MC.(opts.metric);

fig = figure('Name', opts.name, 'Color','w', ...
             'Position',[70 70 min(1500, 380+300*nf) 430]);

%% ---- Variance attribution ------------------------------------------
subplot(1, nf+1, 1);
V = zeros(numel(S), nf);
for c = 1:numel(S), V(c,:) = 100*S(c).var_frac; end
b = bar(V.', 'grouped'); grid on
for c = 1:numel(b), b(c).FaceColor = sty.colors{1+mod(c-1,numel(sty.colors))}; end
set(gca, 'XTickLabel', nm, 'XTickLabelRotation', 30);
ylabel('share of variance  [%]');
title('Attribution of the spread');
legend({S.tag}, 'Location','best');

%% ---- One scatter per factor ----------------------------------------
for j = 1:nf
    subplot(1, nf+1, j+1); hold on; grid on
    for c = 1:size(Y,2)
        y  = Y(:,c);
        ok = ~isnan(y);
        col = sty.colors{1+mod(c-1,numel(sty.colors))};
        scatter(X(ok,j), y(ok), 16, 'filled', 'MarkerFaceColor', col, ...
                'MarkerFaceAlpha', 0.45, 'DisplayName', S(c).tag);
        p  = polyfit(X(ok,j), y(ok), 1);
        xr = linspace(min(X(:,j)), max(X(:,j)), 50);
        plot(xr, polyval(p, xr), '-', 'Color', col, 'LineWidth', 2.6, ...
             'HandleVisibility','off');
    end
    xlabel(nm{j}); ylabel(opts.metric, 'Interpreter','none');
    title(sprintf('SRC %.2f / %.2f', S(1).src(j), S(min(2,numel(S))).src(j)));
    if j == 1, legend('Location','best'); end
end
end
