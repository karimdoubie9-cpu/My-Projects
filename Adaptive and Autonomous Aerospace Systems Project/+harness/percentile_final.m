function q = percentile_final(x,p)
%PERCENTILE_FINAL  Linear-interpolation percentile without Statistics Toolbox.
%
% NaNs are removed before calculation.

arguments
    x (:,1) double
    p (1,:) double {mustBeInRange(p,0,100)}
end

x = sort(x(~isnan(x)));

if isempty(x)
    q = nan(size(p));
    return
end

if isscalar(x)
    q = repmat(x,size(p));
    return
end

pos = 100*((1:numel(x))-0.5)/numel(x);
q   = interp1(pos,x,p,'linear','extrap');
q   = min(max(q,x(1)),x(end));
end
