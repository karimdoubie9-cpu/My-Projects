function sig = constant(t, value)
%CONSTANT  Constant vector signal in the standard From-Workspace format.
%
%   sig = disturb.constant(t, [Fx; Fy; Fz])
%
%   value may have any number of channels.  The same value is applied from
%   the first to the last sample in t.
%
%   This is useful when the uncertainty model explicitly assumes a constant
%   exogenous term, rather than a wind-velocity disturbance whose resulting
%   aerodynamic force changes with relative airspeed.

arguments
    t     (:,1) double
    value (:,1) double
end

if any(diff(t) <= 0)
    error('disturb:constant:nonMonotonicTime', ...
        'Time vector must be strictly increasing.');
end

n = numel(value);

sig.time               = t(:);
sig.signals.values     = repmat(value(:).', numel(t), 1);
sig.signals.dimensions = n;
end
