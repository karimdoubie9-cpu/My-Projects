function M = metrics(R, opts)
%METRICS  Tracking-error metrics for one run, on a uniform time grid.
%
%   M = harness.metrics(R)
%   M = harness.metrics(R, window=[60 80], ss_fraction=0.4)
%
%   Fields:
%     M.t        uniform evaluation grid
%     M.ep,M.ev  N-by-3 per-axis position / velocity errors (reference - actual)
%     M.nep,M.nev  their 2-norms
%     M.rms_p, M.peak_p, M.ss_p   scalar summaries
%     M.rms_v, M.peak_v
%     M.window_peak_p   peak within the requested window, if one is given
%
%   ss_p is the mean error norm over the last ss_fraction of the run, i.e. the
%   residual the controller settles to once transients have decayed -- the
%   figure of merit that actually separates the baseline from the adaptive
%   augmentation, since both eventually stop moving but not at the same offset.

arguments
    R                struct
    opts.n           (1,1) double = 4000
    opts.window      (1,2) double = [NaN NaN]
    opts.ss_fraction (1,1) double = 0.4
end

[tu, iu] = unique(R.t, 'stable');
t = linspace(tu(1), tu(end), opts.n).';

M.t   = t;
M.ep  = interp1(tu, R.p_d(iu,:) - R.p(iu,:), t);
M.ev  = interp1(tu, R.v_d(iu,:) - R.v(iu,:), t);
M.nep = vecnorm(M.ep, 2, 2);
M.nev = vecnorm(M.ev, 2, 2);

M.rms_p  = sqrt(mean(M.nep.^2));
M.peak_p = max(M.nep);
M.rms_v  = sqrt(mean(M.nev.^2));
M.peak_v = max(M.nev);

iss    = t > (1 - opts.ss_fraction)*t(end);
M.ss_p = mean(M.nep(iss));

if ~any(isnan(opts.window))
    iw = t >= opts.window(1) & t <= opts.window(2);
    M.window_peak_p = max(M.nep(iw));
    M.window_rms_p  = sqrt(mean(M.nep(iw).^2));
end

if isfield(R, 'label'), M.label = R.label; else, M.label = ''; end
end
