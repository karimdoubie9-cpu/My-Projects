function sig = step(t, opts)
%STEP  Sustained force step, for disturbance-rejection testing.
%
%   sig = disturb.step(t, t0=60, value=[3;0;0])
%   sig = disturb.step(t, t0=60, value=[3;0;0], duration=20)
%
%   Zero before t0, then constant at VALUE until t0+duration (default: to the
%   end of the run).  Feed it to the Disturbance_i channel, which the model
%   reads as a force in the INERTIAL frame and rotates into the body before
%   adding it to the actuator wrench.
%
%   Why a step rather than the 1-cosine gust: a sustained force is the classic
%   disturbance-rejection test.  The proportional term alone leaves a standing
%   offset against a constant force, so what the step exposes is the INTEGRAL
%   action driving that offset back to zero, plus the tilt the vehicle must
%   hold to generate the opposing force.  disturb.gust is the right choice when
%   the transient response is the question instead.
%
%   The step is expressed with a duplicated sample at t0 (and at the end of the
%   window), so the From-Workspace block's interpolation cannot turn the edge
%   into a ramp spread over one sample interval.
%
%   See also disturb.gust, disturb.fault, disturb.none.

arguments
    t             (:,1) double
    opts.t0       (1,1) double = 60          % [s] onset
    opts.value    (3,1) double = [3; 0; 0]   % [N] force in the inertial frame
    opts.duration (1,1) double = inf         % [s] length; inf = to the end
end

t_end = t(end);
t1    = opts.t0 + opts.duration;

% Breakpoints: start, both sides of the onset edge, both sides of the release
% edge (when there is one), end.
eps_t = 1e-9;
tk = [t(1); opts.t0 - eps_t; opts.t0];
vk = [zeros(2,3); opts.value.'];

if isfinite(t1) && t1 < t_end
    tk = [tk; t1 - eps_t; t1];
    vk = [vk; opts.value.'; zeros(1,3)];
end

tk = [tk; t_end];
vk = [vk; vk(end,:)];

% Guard against an onset at or before the first sample, which would otherwise
% leave a non-monotonic time vector.
keep = [true; diff(tk) > 0];
tk   = tk(keep);
vk   = vk(keep,:);

sig.time               = tk;
sig.signals.values     = vk;
sig.signals.dimensions = 3;
end
