function sig = none(t, n, value)
%NONE  Neutral From-Workspace signal spanning the whole run.
%
%   sig = disturb.none(t, 3)      -> 3 channels of zeros (no wind)
%   sig = disturb.none(t, 6, 1)   -> 6 channels of ones  (nominal rotors)
%
%   Every scenario must define wind_vel, wind_force AND lambda_sig, otherwise
%   the From-Workspace blocks fall back to whatever is left in the workspace.
%   Building the neutral case explicitly is what makes "I forgot to clear the
%   disturbance from the previous run" impossible.

arguments
    t     (:,1) double
    n     (1,1) double {mustBeInteger, mustBePositive}
    value (1,1) double = 0
end

% Two samples spanning the run are enough: the block holds/interpolates between.
sig.time               = [t(1); t(end)];
sig.signals.values     = value * ones(2, n);
sig.signals.dimensions = n;
end
