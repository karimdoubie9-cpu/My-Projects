function S = sensitivity(MC, opts)
%SENSITIVITY  Which uncertainty drives the spread? Costs no extra simulations.
%
%   S = harness.sensitivity(MC)
%   S = harness.sensitivity(MC, metric='rms_p', quiet=true)
%
%   A joint Monte Carlo varies every factor at once, so no single scatter plot
%   separates their effects -- each point is confounded by the other factor.
%   This recovers the individual contributions from the SAME samples by
%   regression, which is why it needs no extra runs.
%
%   INPUT FACTORS
%     mass        the payload mass draw
%     mean(lam)   average rotor effectiveness -- the collective thrust the
%                 vehicle can still produce
%     min(lam)    the worst single rotor -- what actually binds the allocation,
%                 since the pseudo-inverse has to redistribute around it
%     spread(lam) max-min, i.e. how ASYMMETRIC the degradation is; an even loss
%                 across all six rotors is a very different disturbance from one
%                 dead arm, and only this regressor can tell them apart
%
%   MEASURES
%     SRC   Standardized Regression Coefficient. Regress the standardized metric
%           on the standardized factors: beta_j is then the change in the metric
%           (in standard deviations) per one-standard-deviation change in factor
%           j, holding the others fixed. That "holding the others fixed" is
%           exactly what a raw scatter plot cannot do.
%           For a near-linear response, beta_j^2 is the FRACTION OF VARIANCE
%           explained by factor j, so the betas partition the spread.
%           Trust this only when the model R^2 is high (say > 0.8) -- it is
%           reported so the reader can judge.
%
%     rho   Spearman rank correlation: monotonic association, robust to
%           nonlinearity and to outliers. Where |rho| is large but SRC is small,
%           the relationship is monotonic but not linear.
%
%   Returns S(c) per configuration with .names .src .rho .r2 .var_frac.

arguments
    MC          struct
    opts.metric (1,:) char = 'rms_p'
    opts.quiet  (1,1) logical = false
end

lam = MC.samples.lambda;
X   = [ MC.samples.mass(:), ...
        mean(lam, 2), ...
        min(lam, [], 2), ...
        max(lam, [], 2) - min(lam, [], 2) ];
names = {'mass', 'mean(lambda)', 'min(lambda)', 'spread(lambda)'};

% Drop factors that were held constant (a one-at-a-time campaign) -- a
% zero-variance regressor makes the design matrix singular.
keep  = std(X, 0, 1) > 1e-12;
X     = X(:, keep);
names = names(keep);

Y  = MC.(opts.metric);
nc = size(Y, 2);
S  = struct('tag',{},'names',{},'src',{},'rho',{},'r2',{},'var_frac',{});

for c = 1:nc
    y  = Y(:,c);
    ok = ~isnan(y) & all(isfinite(X), 2);
    xc = X(ok,:);
    yc = y(ok);

    Z  = (xc - mean(xc)) ./ std(xc);
    zy = (yc - mean(yc)) / std(yc);

    beta = [ones(size(Z,1),1), Z] \ zy;    % intercept + standardized slopes
    src  = beta(2:end).';

    resid = zy - [ones(size(Z,1),1), Z]*beta;
    r2    = 1 - sum(resid.^2)/sum((zy - mean(zy)).^2);

    rho = zeros(1, size(xc,2));
    for j = 1:size(xc,2)
        rho(j) = corr_(rank_(xc(:,j)), rank_(yc));
    end

    S(c).tag      = MC.tags{c};
    S(c).names    = names;
    S(c).src      = src;
    S(c).rho      = rho;
    S(c).r2       = r2;
    S(c).var_frac = src.^2 / max(sum(src.^2), eps);   % normalised contribution
    S(c).n        = numel(yc);
end

if opts.quiet, return, end

%% ---- Report ---------------------------------------------------------
vary = 'both';
if isfield(MC.samples,'vary'), vary = MC.samples.vary; end
fprintf('\n=== Sensitivity of %s (campaign varied: %s) ===\n', opts.metric, vary);
for c = 1:numel(S)
    fprintf('\n%s   [n = %d, model R^2 = %.3f]\n', S(c).tag, S(c).n, S(c).r2);
    fprintf('  %-16s %9s %9s %14s\n','factor','SRC','rho','share of var');
    fprintf('  %s\n', repmat('-',1,52));
    for j = 1:numel(S(c).names)
        fprintf('  %-16s %9.3f %9.3f %12.1f%%\n', S(c).names{j}, ...
                S(c).src(j), S(c).rho(j), 100*S(c).var_frac(j));
    end
    if S(c).r2 < 0.8
        fprintf(['  NOTE: R^2 = %.2f is low, so the response is not well described\n' ...
                 '        by a linear fit -- read rho rather than SRC here.\n'], S(c).r2);
    end
end
end

% ---------------------------------------------------------------------
function r = rank_(x)
%RANK_  Ranks with ties averaged (Spearman needs midranks).
n        = numel(x);
[xs, ix] = sort(x(:));
r        = zeros(n,1);
r(ix)    = 1:n;
i = 1;
while i <= n
    j = i;
    while j < n && xs(j+1) == xs(i), j = j + 1; end
    if j > i, r(ix(i:j)) = mean(i:j); end
    i = j + 1;
end
end

function c = corr_(a, b)
a = a(:) - mean(a);
b = b(:) - mean(b);
c = (a.'*b) / max(sqrt((a.'*a)*(b.'*b)), eps);
end
