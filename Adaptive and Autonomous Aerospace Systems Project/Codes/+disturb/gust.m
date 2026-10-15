function sig = gust(t, opts)
%GUST  1-cosine ("1-minus-cosine") discrete gust, the classic airworthiness shape.
%
%   sig = disturb.gust(t, t0=60, duration=20, peak=[8;0;0])
%
%   Amplitude profile over the gust window tau = t - t0, 0 <= tau <= duration:
%
%       s(tau) = 0.5*(1 - cos(2*pi*tau/duration))
%
%   which rises smoothly from zero, peaks at mid-window and returns to zero, so
%   the vehicle is not hit with a discontinuity that would make the response an
%   artefact of the input rather than of the controller.
%
%   Use it for BOTH wind channels the model reads:
%     * wind_vel   -- a velocity gust [m/s], which the aerodynamic block turns
%                     into a force through D_fv (so the force depends on the
%                     relative airspeed, the physically meaningful version);
%     * wind_force -- a force gust [N] injected directly.
%   Pass whichever the scenario calls for; the other stays at disturb.none.

arguments
    t             (:,1) double
    opts.t0       (1,1) double = 60      % [s] gust start
    opts.duration (1,1) double = 20      % [s] gust length
    opts.peak     (3,1) double = [8;0;0] % peak amplitude, per axis
end

tau   = t - opts.t0;
inside = tau >= 0 & tau <= opts.duration;

shape         = zeros(numel(t), 1);
shape(inside) = 0.5*(1 - cos(2*pi*tau(inside)/opts.duration));

sig.time               = t;
sig.signals.values     = shape * opts.peak.';   % N-by-3
sig.signals.dimensions = 3;
end
