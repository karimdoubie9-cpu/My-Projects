function q = percentile(x, p)
%PERCENTILE  Linear-interpolation percentile without the Statistics Toolbox.
%
%   q = harness.percentile(x, 95)
%   q = harness.percentile(x, [5 50 95])
%
%   Uses the midpoint plotting positions (i-0.5)/n, which is the convention
%   prctile follows, so results match it where the toolbox is available.
%   NaNs are dropped rather than propagated: a diverged Monte Carlo run should
%   reduce the sample count, not poison the whole statistic.

arguments
    x (:,1) double
    p (1,:) double {mustBeInRange(p, 0, 100)}
end

x = sort(x(~isnan(x)));

if isempty(x)
    q = nan(size(p));
    return
end
if isscalar(x)
    q = repmat(x, size(p));
    return
end

pos = 100*((1:numel(x)) - 0.5)/numel(x);
q   = interp1(pos, x, p, 'linear', 'extrap');
q   = min(max(q, x(1)), x(end));    % never extrapolate past the sample range
end
