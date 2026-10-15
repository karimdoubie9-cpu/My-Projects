function [fig, info] = convergence(MC, opts)
%CONVERGENCE  Is the sample count large enough? Running-statistics diagnostic.
%
%   [fig, info] = viz.convergence(MC)
%   [fig, info] = viz.convergence(MC, metric='rms_p', tol=0.02)
%
%   Three panels, answering three different questions:
%
%   1. RUNNING MEAN with a 95% confidence band.
%      The estimate after k samples, mean(x(1:k)), with band +-1.96*s_k/sqrt(k).
%      The curve should flatten and the band should close around it. This is the
%      qualitative "has it settled" check.
%
%   2. RELATIVE STANDARD ERROR, the quantitative criterion.
%      RSE(k) = 1.96*s_k/(sqrt(k)*mean_k) -- the 95% half-width as a fraction of
%      the estimate. Convergence is declared when it drops below `tol`
%      (default 2%). Because RSE falls as 1/sqrt(k), the required sample size
%      can be extrapolated in closed form:
%
%          n_required = ( 1.96 * s / (tol * mean) )^2
%
%      which is reported in `info` and printed. That number, not a visual
%      impression of flatness, is what justifies N in the report.
%
%   3. RUNNING p95, the tail statistic.
%      Percentiles converge more slowly than means, so a campaign can be plenty
%      for the mean and still noisy at p95. Since the summary table quotes p95,
%      it has to be shown to have settled too.
%
%   ORDERING SENSITIVITY: the samples are i.i.d., so any ordering is a valid
%   history, but a single one can flatten early by luck. Faint grey curves show
%   `shuffles` random re-orderings; the spread between them is the honest
%   uncertainty in the convergence history itself.

arguments
    MC            struct
    opts.metric   (1,:) char = 'rms_p'
    opts.tol      (1,1) double = 0.02    % target relative half-width (2%)
    opts.shuffles (1,1) double = 30
    opts.seed     (1,1) double = 0
    opts.name     (1,:) char = 'Monte Carlo convergence'
end

X  = MC.(opts.metric);
nc = size(X, 2);
s  = viz.style();
z  = 1.96;

fig = figure('Name', opts.name, 'Color','w', 'Position',[80 80 1250 420]);
info = struct('config', {}, 'tag', {}, 'n', {}, 'mean', {}, 'std', {}, ...
              'rse', {}, 'n_required', {}, 'converged', {});

rng(opts.seed);

for c = 1:nc
    x = X(:,c);
    x = x(~isnan(x));
    n = numel(x);
    k = (1:n).';

    % Running mean and running standard deviation (unbiased, k>=2)
    rmean = cumsum(x)./k;
    rvar  = max((cumsum(x.^2) - k.*rmean.^2)./max(k-1,1), 0);
    rstd  = sqrt(rvar);
    sem   = rstd./sqrt(k);
    rse   = z*sem./abs(rmean);

    % Closed-form requirement from the final (best) variance estimate
    n_req = ceil((z*rstd(end)/(opts.tol*abs(rmean(end))))^2);

    info(c).config     = MC.configs(c);
    info(c).tag        = MC.tags{c};
    info(c).n          = n;
    info(c).mean       = rmean(end);
    info(c).std        = rstd(end);
    info(c).rse        = rse(end);
    info(c).n_required = n_req;
    info(c).converged  = rse(end) <= opts.tol;

    col = s.colors{1+mod(c-1,numel(s.colors))};

    %% Panel 1: running mean with CI band
    subplot(1,3,1); hold on; grid on
    hi = rmean + z*sem;  lo = rmean - z*sem;
    valid = k >= 2;
    fill([k(valid); flipud(k(valid))], [hi(valid); flipud(lo(valid))], col, ...
         'FaceAlpha', 0.15, 'EdgeColor','none', 'HandleVisibility','off');
    plot(k, rmean, 'LineWidth', 2.3, 'Color', col, 'DisplayName', MC.tags{c});
    yline(rmean(end), ':', 'Color', col, 'HandleVisibility','off');

    %% Panel 2: relative standard error
    subplot(1,3,2); hold on; grid on
    % Ordering sensitivity, drawn first so the main curve sits on top
    for r = 1:opts.shuffles
        xs = x(randperm(n));
        m2 = cumsum(xs)./k;
        v2 = max((cumsum(xs.^2) - k.*m2.^2)./max(k-1,1), 0);
        plot(k, z*sqrt(v2)./sqrt(k)./abs(m2), '-', 'Color', [col 0.10], ...
             'LineWidth', 0.7, 'HandleVisibility','off');
    end
    plot(k, rse, 'LineWidth', 2.3, 'Color', col, 'DisplayName', MC.tags{c});

    %% Panel 3: running p95
    subplot(1,3,3); hold on; grid on
    % Start once there are enough samples for a percentile to mean anything,
    % but never past the end: a short campaign (n < 10, e.g. a smoke test) would
    % otherwise produce an empty grid and index into nothing.
    step  = max(1, floor(n/200));
    kk    = unique([(min(10,n):step:n).'; n]);
    if ~isempty(kk)
        rp95 = arrayfun(@(j) harness.percentile(x(1:j), 95), kk);
        plot(kk, rp95, 'LineWidth', 2.3, 'Color', col, 'DisplayName', MC.tags{c});
        yline(rp95(end), ':', 'Color', col, 'HandleVisibility','off');
    end
end

%% ---- Panel cosmetics ------------------------------------------------
subplot(1,3,1)
xlabel('samples used'); ylabel(sprintf('running mean of %s', opts.metric));
title('Running mean \pm 95% CI', 'Interpreter','tex');
legend('Location','best'); set(gca,'XScale','log')

subplot(1,3,2)
yline(opts.tol, '--k', sprintf('%.0f%% target', 100*opts.tol), ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');
xlabel('samples used'); ylabel('relative 95% half-width');
title('Relative standard error (faint = re-orderings)');
set(gca,'XScale','log','YScale','log'); legend('Location','best');

subplot(1,3,3)
xlabel('samples used'); ylabel(sprintf('running p95 of %s', opts.metric));
title('Tail statistic (p95)'); legend('Location','best'); set(gca,'XScale','log')

%% ---- Verdict --------------------------------------------------------
fprintf('\n=== Convergence of %s (target: %.0f%% relative 95%% half-width) ===\n', ...
        opts.metric, 100*opts.tol);
fprintf('%-10s %7s %11s %11s %10s %12s %s\n', ...
        'case','n','mean','std','RSE','n required','verdict');
fprintf('%s\n', repmat('-',1,76));
for c = 1:numel(info)
    fprintf('%-10s %7d %11.4f %11.4f %9.2f%% %12d   %s\n', ...
        info(c).tag, info(c).n, info(c).mean, info(c).std, 100*info(c).rse, ...
        info(c).n_required, ternary(info(c).converged,'converged','MORE RUNS NEEDED'));
end

worst = max([info.n_required]);
if all([info.converged])
    fprintf(['\n%d samples is sufficient: every case is inside the %.0f%% band, ' ...
             'and the\nbinding requirement was %d.\n'], info(1).n, 100*opts.tol, worst);
else
    fprintf(['\nNot converged at %d samples. The binding case needs about %d ' ...
             '(n scales as\n1/tol^2, so halving the tolerance quadruples the ' ...
             'runs).\n'], info(1).n, worst);
end
end

% ---------------------------------------------------------------------
function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
