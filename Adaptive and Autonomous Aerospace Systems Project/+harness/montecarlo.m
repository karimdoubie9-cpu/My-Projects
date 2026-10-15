function MC = montecarlo(ref, opts)
%MONTECARLO  Parametric uncertainty campaign for Task 3.4.
%
%   MC = harness.montecarlo(ref, n=200)
%   MC = harness.montecarlo(ref, n=500, mrac=[0 1], parallel=true)
%
%   Samples the two uncertainties named in the exam text and runs the scenario
%   once per sample, for each requested controller configuration:
%
%     package mass      m_p ~ N(mu, sigma), truncated to [lo, hi]
%                       Default: mu = 0.4 kg, sigma = 0.08 kg, clipped to
%                       [0.15, 0.65] -- a +-60% spread, which is realistic for a
%                       delivery payload whose contents are not weighed.
%
%     rotor efficiency  lambda_r ~ N(mu, sigma) i.i.d. per rotor, truncated to
%                       [lo, 1].  Default: mu = 0.90, sigma = 0.07, clipped to
%                       [0.55, 1] -- wear and fouling degrade a rotor, they
%                       never make it produce more thrust than commanded, hence
%                       the hard upper bound at 1.
%
%   Truncation (rather than an unbounded normal) keeps every draw physically
%   admissible: a negative mass or an effectiveness above unity would be
%   meaningless, and a single such draw would poison the statistics.
%
%   Returns MC.results (nRuns-by-nCfg struct array of metrics), MC.samples
%   (the drawn parameters, so any run can be reproduced) and MC.summary.
%
%   parallel=true uses parsim, which is the reason harness.input builds a
%   self-contained SimulationInput per run: the workers have no shared base
%   workspace to inherit from.

arguments
    ref                struct
    opts.n             (1,1) double {mustBeInteger,mustBePositive} = 200
    opts.mrac          (1,:) double = [0 1]
    opts.package       (1,1) double = 1
    opts.mass_mu       (1,1) double = 0.4
    opts.mass_sigma    (1,1) double = 0.08
    opts.mass_range    (1,2) double = [0.15 0.65]
    opts.lambda_mu     (1,1) double = 0.90
    opts.lambda_sigma  (1,1) double = 0.07
    opts.lambda_range  (1,2) double = [0.55 1.0]
    opts.vary          (1,:) char {mustBeMember(opts.vary, ...
                         {'both','mass','lambda'})} = 'both'
    opts.wind          struct = struct([])
    opts.seed          (1,1) double = 42
    opts.parallel      (1,1) logical = false
    opts.fast_restart  (1,1) logical = true
    opts.progress      (1,:) char {mustBeMember(opts.progress, ...
                         {'auto','inplace','lines','none'})} = 'auto'
    opts.save_to       (1,:) char = 'auto'   % 'auto' = timestamped, '' = no save
    opts.model         (1,:) char = 'Final'
end

rng(opts.seed);
n    = opts.n;
cfgs = opts.mrac;
nc   = numel(cfgs);

%% ---- Draw the samples ----------------------------------------------
mass = harness.truncnorm(opts.mass_mu,   opts.mass_sigma,   opts.mass_range,   [n 1]);
lam  = harness.truncnorm(opts.lambda_mu, opts.lambda_sigma, opts.lambda_range, [n 6]);

% One-at-a-time campaigns: freeze the other factor at its DISTRIBUTION MEAN,
% not at a "perfect" value, so all three campaigns sit at the same nominal
% operating point and any difference in spread is attributable to the factor
% being varied rather than to a shift in the mean condition.
switch opts.vary
    case 'mass'
        lam(:)  = opts.lambda_mu;
    case 'lambda'
        mass(:) = opts.mass_mu;
end

MC.samples.mass   = mass;
MC.samples.lambda = lam;
MC.samples.seed   = opts.seed;
MC.samples.vary   = opts.vary;

%% ---- Build every SimulationInput up front --------------------------
nTotal = n*nc;
fprintf('Monte Carlo: %d samples x %d configuration(s) = %d runs\n', n, nc, nTotal);

buildStep = harness.progress(nTotal, label = 'building inputs', ...
                             min_interval = 0.5, mode = opts.progress);
inputs(nTotal) = Simulink.SimulationInput(opts.model);
idx = 0;
for c = 1:nc
    for i = 1:n
        idx = idx + 1;
        P   = harness.config(package_mass = mass(i));
        % A constant per-rotor effectiveness for the whole run: this campaign
        % is about parametric spread, not about fault timing.
        lamSig = disturb.fault(ref.time, t0 = 0, lambda = lam(i,:).');
        inputs(idx) = harness.input(ref, config = P, mrac = cfgs(c), ...
                                    package = opts.package, lambda = lamSig, ...
                                    wind = opts.wind, model = opts.model);
        buildStep(idx);
    end
end

%% ---- Run -----------------------------------------------------------
tic
if opts.parallel
    % No base-workspace transfer is needed or requested: harness.input already
    % put every variable the model resolves onto each SimulationInput, which is
    % exactly what makes these runs safe to farm out to workers.
    outs = num2cell(parsim(inputs, ShowProgress = 'on', StopOnError = 'off'));
else
    % Fast restart keeps the model compiled across the batch (~2.3x here). The
    % guard restores the previous setting even if a run throws.
    if opts.fast_restart
        % Deliberately unread: this is an onCleanup handle whose destructor
        % restores the FastRestart setting when the function exits or throws.
        frGuard = harness.fastrestart(opts.model); %#ok<NASGU>
    end
    runStep = harness.progress(nTotal, label = 'simulating', mode = opts.progress);
    outs    = cell(1, nTotal);
    nfail   = 0;
    for k = 1:nTotal
        try
            outs{k} = sim(inputs(k));
        catch ME
            outs{k} = ME;    % recorded as a failure in the reduction below
            nfail   = nfail + 1;
        end
        runStep(k);
    end
    if nfail > 0
        fprintf('  %d of %d runs threw during simulation\n', nfail, nTotal);
    end
end
MC.elapsed = toc;

%% ---- Reduce to metrics ---------------------------------------------
fields = {'rms_p','peak_p','ss_p','rms_v'};
for f = 1:numel(fields), MC.(fields{f}) = nan(n, nc); end
MC.failed = false(n, nc);

redStep = harness.progress(nTotal, label = 'reducing metrics', ...
                           min_interval = 0.5, mode = opts.progress);
idx = 0;
for c = 1:nc
    for i = 1:n
        idx = idx + 1;
        try
            m = harness.metrics(harness.unpack(outs{idx}));
            for f = 1:numel(fields), MC.(fields{f})(i,c) = m.(fields{f}); end
        catch
            MC.failed(i,c) = true;   % diverged or errored: kept as NaN
        end
        redStep(idx);
    end
end

MC.configs = cfgs;
MC.tags    = arrayfun(@(c) ternary(c,'PBMRAC','Baseline'), cfgs, UniformOutput=false);

%% ---- Persist --------------------------------------------------------
% Saved BEFORE the summary is printed and before any caller starts plotting:
% a campaign of several hundred runs is expensive enough that a downstream
% error must not be able to destroy it. Everything needed to redraw the figures
% or reproduce the campaign is in the struct; the raw SimulationOutputs are
% deliberately not kept (400 runs of logged states is gigabytes, and every
% metric derived from them is already here).
MC.ref  = ref;
MC.meta = struct('timestamp',   datetime('now'), ...
                 'model',       opts.model, ...
                 'matlab',      version('-release'), ...
                 'n',           n, ...
                 'configs',     cfgs, ...
                 'parallel',    opts.parallel, ...
                 'fast_restart',opts.fast_restart, ...
                 'seed',        opts.seed, ...
                 'vary',        opts.vary, ...
                 'mass_dist',   [opts.mass_mu   opts.mass_sigma   opts.mass_range], ...
                 'lambda_dist', [opts.lambda_mu opts.lambda_sigma opts.lambda_range], ...
                 'elapsed_s',   MC.elapsed);

if ~isempty(opts.save_to)
    file = opts.save_to;
    if strcmp(file, 'auto')
        % Timestamped by default: a later 20-run sanity check must not silently
        % overwrite the 400-run campaign that took three quarters of an hour.
        file = sprintf('montecarlo_%s.mat', ...
                       char(datetime('now', Format='yyyyMMdd_HHmmss')));
    end
    [~,~,ext] = fileparts(file);
    if isempty(ext), file = [file '.mat']; end
    if ~isAbsolute(file), file = fullfile(pwd, file); end
    MC.meta.file = file;
    try
        save(file, 'MC', '-v7.3');
        d = dir(file);
        fprintf('\nSaved results to %s (%.1f MB)\n', file, d.bytes/2^20);
    catch ME
        % The data is still in the returned struct, so a failed write is a
        % warning, not an error -- do not throw away a completed campaign.
        warning('harness:montecarlo:saveFailed', ...
                'Could not write %s (%s). Results are still returned in MC.', ...
                file, ME.message);
    end
end

%% ---- Summary -------------------------------------------------------
% Shared with the reload path in run_task3_montecarlo.m, so a campaign read back
% from disk reports identically to one just computed.
harness.mc_summary(MC);
end

% ---------------------------------------------------------------------
function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function tf = isAbsolute(p)
%ISABSOLUTE  True for a drive-qualified or UNC Windows path, or a POSIX root.
tf = ~isempty(regexp(p, '^([A-Za-z]:[\\/]|\\\\|/)', 'once'));
end
