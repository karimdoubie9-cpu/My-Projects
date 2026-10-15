function sig = fault(t, opts)
%FAULT  Time-gated loss of rotor effectiveness (lambda_sig, 6 channels).
%
%   sig = disturb.fault(t, t0=90, lambda=[0.8 0.6 1 1 0.4 1])
%   sig = disturb.fault(t, t0=0,  lambda=0.7)        % all rotors, from t=0
%
%   lambda_r = 1 is a healthy rotor; lambda_r < 1 scales the thrust that rotor
%   actually produces for a given command.  Before t0 all rotors are nominal;
%   from t0 onward the given vector applies, held to the end of the run.
%
%   The transition is a genuine step: a degradation that ramps in slowly would
%   be trivially absorbed by the integral action and would not test the
%   adaptation at all.
%
%   With the default vector rotors 2 and 5 (opposing arms) lose the most
%   authority, which biases the allocation and forces a visible attitude
%   correction rather than a pure thrust deficit.

arguments
    t           (:,1) double
    opts.t0     (1,1) double = 90
    opts.lambda (:,1) double = [0.8; 0.6; 0.78; 0.82; 0.4; 0.92]
end

lam = opts.lambda;
if isscalar(lam), lam = lam*ones(6,1); end
if numel(lam) ~= 6
    error('disturb:fault:size', 'lambda must have 6 entries (one per rotor).');
end

degraded = t >= opts.t0;
values   = repmat(ones(1,6), numel(t), 1);
values(degraded, :) = repmat(lam.', nnz(degraded), 1);

sig.time               = t;
sig.signals.values     = values;
sig.signals.dimensions = 6;
end
