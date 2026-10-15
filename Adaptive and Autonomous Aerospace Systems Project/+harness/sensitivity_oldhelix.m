function S = sensitivity_oldhelix(MC, opts)
%SENSITIVITY  Which uncertainty drives the spread? Costs no extra simulations.
%
%   S = harness.sensitivity_oldhelix(MC)
%   S = harness.sensitivity_oldhelix(MC, metric='rms_p', quiet=true)
%
%   A joint Monte Carlo varies all uncertain factors in the same run.  The
%   standardized regression below recovers each factor's contribution while
%   holding the others fixed, without requiring a second campaign.
%
%   INPUT FACTORS
%     mass          payload mass
%     r_xy          horizontal package-CoM lever arm sqrt(r_cx^2+r_cy^2)
%     depth         distance of the package CoM below the body origin, -r_cz
%     mean(lam)     average rotor effectiveness
%     min(lam)      worst single-rotor effectiveness
%     spread(lam)   max(lam)-min(lam), measuring asymmetry
%
%   r_xy is used instead of signed r_cx/r_cy because a symmetric +/- lateral
%   offset can have a strong physical effect while having almost zero linear
%   correlation with a signed coordinate.  The actual sampled x/y/z coordinates
%   remain stored in MC.samples.r_c for exact reproduction of every run.
%
%   Older campaign files without MC.samples.r_c remain supported: in that case
%   the original four mass/lambda regressors are used.
%
%   MEASURES
%     SRC   Standardized Regression Coefficient.
%     rho   Spearman rank correlation.
%
%   Returns S(c) per configuration with .names .src .rho .r2 .var_frac.

arguments
    MC          struct
    opts.metric (1,:) char = 'rms_p'
    opts.quiet  (1,1) logical = false
end

lam = MC.samples.lambda;

X = MC.samples.mass(:);
names = {'mass'};

if isfield(MC.samples,'r_c')
    rc = MC.samples.r_c;
    if size(rc,2) ~= 3
        error('harness:sensitivity_oldhelix:badRc', ...
              'MC.samples.r_c must be n-by-3 [x_c y_c z_c].');
    end
    rxy   = hypot(rc(:,1), rc(:,2));
    depth = -rc(:,3);

    X = [X, rxy, depth]; %#ok<AGROW>
    names = [names, {'r_xy(CoM)', 'depth(CoM)'}]; %#ok<AGROW>
end

X = [X, ...
     mean(lam, 2), ...
     min(lam, [], 2), ...
     max(lam, [], 2) - min(lam, [], 2) ];

names = [names, {'mean(lambda)', 'min(lambda)', 'spread(lambda)'}];

% Drop factors that were held constant. A zero-variance regressor makes the
% standardized design matrix singular.
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

    if numel(yc) < 3 || std(yc) <= eps
        S(c).tag      = MC.tags{c};
        S(c).names    = names;
        S(c).src      = nan(1,size(xc,2));
        S(c).rho      = nan(1,size(xc,2));
        S(c).r2       = NaN;
        S(c).var_frac = nan(1,size(xc,2));
        S(c).n        = numel(yc);
        continue
    end

    Z  = (xc - mean(xc,1)) ./ std(xc,0,1);
    zy = (yc - mean(yc)) / std(yc);

    beta = [ones(size(Z,1),1), Z] \ zy;
    src  = beta(2:end).';

    resid = zy - [ones(size(Z,1),1), Z]*beta;
    den   = sum((zy - mean(zy)).^2);
    r2    = 1 - sum(resid.^2)/max(den, eps);

    rho = zeros(1, size(xc,2));
    for j = 1:size(xc,2)
        rho(j) = corr_(rank_(xc(:,j)), rank_(yc));
    end

    S(c).tag      = MC.tags{c};
    S(c).names    = names;
    S(c).src      = src;
    S(c).rho      = rho;
    S(c).r2       = r2;
    S(c).var_frac = src.^2 / max(sum(src.^2), eps);
    S(c).n        = numel(yc);
end

if opts.quiet, return, end

%% ---- Report ---------------------------------------------------------
vary = 'all';
if isfield(MC.samples,'vary'), vary = MC.samples.vary; end
fprintf('\n=== Sensitivity of %s (campaign varied: %s) ===\n', opts.metric, vary);

for c = 1:numel(S)
    fprintf('\n%s   [n = %d, model R^2 = %.3f]\n', S(c).tag, S(c).n, S(c).r2);
    fprintf('  %-18s %9s %9s %14s\n','factor','SRC','rho','share of var');
    fprintf('  %s\n', repmat('-',1,54));

    for j = 1:numel(S(c).names)
        fprintf('  %-18s %9.3f %9.3f %12.1f%%\n', S(c).names{j}, ...
                S(c).src(j), S(c).rho(j), 100*S(c).var_frac(j));
    end

    if isfinite(S(c).r2) && S(c).r2 < 0.8
        fprintf(['  NOTE: R^2 = %.2f is low, so the response is not well described\n' ...
                 '        by a linear fit -- read rho and the raw distributions\n' ...
                 '        rather than treating SRC as a complete attribution.\n'], S(c).r2);
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
