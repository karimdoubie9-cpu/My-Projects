function MC = run_montecarlo_final(opts)
%RUN_MONTECARLO_FINAL  Final crash-safe Task 3.4 Monte Carlo campaign.
%
%   MC = run_montecarlo_final
%   MC = run_montecarlo_final(n=50)     % quick check
%   MC = run_montecarlo_final(resume='results/montecarlo_final_....mat')
%   MC = run_montecarlo_final(parallel=true)            % parsim on a local pool
%   MC = run_montecarlo_final(parallel=true, workers=6, batch=24)
%
% FINAL CAMPAIGN DEFINITION
% -------------------------
% Trajectory (same final helix used in Task 3):
%   radius        = 10 m
%   rate          = 0.2 rad/s
%   climb         = 2 m/s
%   center        = [10;10] m
%   entry phase   = searched by traj.helix (free-entry minimum snap)
%   start_p       = [0;0;0] m
%   approach_time = 20 s
%
% Uncertainties (jointly sampled every run):
%   package mass:
%       uniform on mass_range = [0.15, 0.70] kg
%       (mass_mu / mass_sigma are only used by the commented-out Gaussian option)
%   package CoM:
%       NOT sampled: fixed at the nominal r_c = [0.02; 0.02; -0.15] m
%       (rc_mu / rc_sigma / rc_range are only used by the commented-out option)
%   six independent rotor efficiencies:
%       N(0.90,0.07^2), truncated to [0.65,1.00]
%
% Final PBMRAC tuning (forced explicitly for reproducibility):
%   L_p = L_v = 150
%   Gamma_theta  = diag([144 144 57.6])
%   Gamma_lambda = [0.75; 0.75; 1.5]   (same as +design/adaptive.m)
%   sigma_e      = 0.02  (e-modification on theta)
%   Lambda projection:
%       center = 0.725, r0 = 0.275, r1 = 0.375
%
% Baseline and PBMRAC are run on EXACTLY the same sampled plants.
%
% EXTRA ROBUSTNESS METRICS
% ------------------------
% Besides RMS/peak/steady-state position error and RMS velocity error, each
% run is reduced immediately to:
%   max_tilt_deg        maximum |roll| or |pitch|
%   sat_upper_frac      fraction of samples with any rotor physically at
%                       the upper actuator stop
%   max_omega_eq        maximum post-limit equivalent rotor speed
%   min_lambda_hat      minimum adaptive effectiveness estimate
%   proj_region_frac    fraction of samples with any Lambda_hat < 0.45
%
% The actuator metric is intentionally based on T_out/lambda: it measures
% PHYSICAL upper-limit occupancy after the saturation block. It is not the
% same as reconstructing an unsaturated command from an applied wrench.
%
% The runner preserves the old campaign's crash-safety:
%   * one simulation in memory at a time (parallel: one batch at a time);
%   * scopes and unused loggers disabled only while the campaign runs;
%   * solver untouched;
%   * checkpoint after every run (parallel: after every batch);
%   * resume support.
%
% Recommended final call:
%   MC = run_montecarlo_final(n=500);
%

arguments
    opts.n             (1,1) double {mustBeInteger,mustBePositive} = 500
    opts.duration      (1,1) double {mustBePositive} = 100
    opts.seed          (1,1) double = 42
    opts.mrac          (1,:) double = [0 1]
    opts.package       (1,1) double = 1
    opts.robust_adapt  (1,1) double {mustBeMember(opts.robust_adapt,[0 1])} = 1

    % Uncertainty distributions
    % Package-mass distribution parameters
    opts.mass_mu       (1,1) double = 0.40
    opts.mass_sigma    (1,1) double = 0.12
    opts.mass_range    (1,2) double = [0.15 0.7]

    opts.lambda_mu     (1,1) double = 0.90
    opts.lambda_sigma  (1,1) double = 0.07
    opts.lambda_range  (1,2) double = [0.65 1.00]

    opts.rc_mu         (3,1) double = [0.02; 0.02; -0.15]
    opts.rc_sigma      (3,1) double = [0.02; 0.02; 0.05]
    opts.rc_range      (3,2) double = [-0.03  0.07; ...
                                       -0.03  0.07; ...
                                       -0.30 -0.05]

    % Final helix
    opts.helix_radius      (1,1) double = 10
    opts.helix_rate        (1,1) double = 0.2
    opts.helix_climb       (1,1) double = 2
    opts.helix_center      (2,1) double = [10;10]
    opts.helix_phase0      (1,1) double = 0   % unused: traj.helix now searches the entry phase
    opts.helix_start_p     (3,1) double = [0;0;0]
    opts.helix_approach    (1,1) double = 20

    % Final PBMRAC tuning
    opts.L_p               (1,1) double = 150
    opts.L_v               (1,1) double = 150
    opts.Gamma_theta       (3,3) double = diag([144 144 57.6])
    opts.Gamma_lambda      (3,1) double = [0.75;0.75;1.5]   % z raised x10 (was 0.15)
    opts.sigma_e           (1,1) double = 0.02
    opts.lambda_proj_center (1,1) double = 0.725
    opts.lambda_proj_r0     (1,1) double = 0.275
    opts.lambda_proj_r1     (1,1) double = 0.375

    % Campaign mechanics
    opts.timeout       (1,1) double {mustBePositive} = 300
    opts.min_free_gb   (1,1) double {mustBeNonnegative} = 2
    opts.resume        (1,:) char = ''
    opts.label         (1,:) char = 'montecarlo_final'
    opts.outdir        (1,:) char = 'results'
    opts.plots         (1,1) logical = true
    opts.export        (1,:) char = 'Task3_montecarlo_final'
    opts.model         (1,:) char = 'Final'

    % Parallel execution (Parallel Computing Toolbox). Same samples, same
    % metrics and same checkpoint file as the serial path; only the order in
    % which simulations finish differs.
    opts.parallel      (1,1) logical = false
    opts.workers       (1,1) double {mustBeInteger,mustBeNonnegative} = 0   % 0 = pool default
    opts.batch         (1,1) double {mustBeInteger,mustBeNonnegative} = 0   % runs per parsim call; 0 = 4 x workers
end

%% ========================================================================
%  Preconditions
%  ========================================================================
if bdIsLoaded(opts.model)
    error('run_montecarlo_final:modelOpen', ...
          ['''%s'' is already loaded. Close it (close_system(''%s'',0)) or run ' ...
           'the campaign in a fresh MATLAB process.'], opts.model, opts.model);
end

if ~isvarname(matlab.lang.makeValidName(opts.label))
    error('run_montecarlo_final:badLabel', 'label must be filename-safe.');
end

if ~isfolder(opts.outdir)
    mkdir(opts.outdir);
end

%% ========================================================================
%  Fresh campaign or resume
%  ========================================================================
if isempty(opts.resume)

    validateRanges(opts);

    rng(opts.seed);

    % Preserve historical draw order. With the same seed and distributions,
    % sample i stays directly comparable with the old-helix campaign.
    % ---- Package mass distribution: choose ONE --------------------------

    % OPTION A - truncated Gaussian
    % mass = harness.truncnorm( ...
    %     opts.mass_mu, ...
    %     opts.mass_sigma, ...
    %     opts.mass_range, ...
    %     [opts.n 1]);
    
    % OPTION B - uniform
    mass = opts.mass_range(1) + ...
           diff(opts.mass_range) * rand(opts.n,1);


% ---- Rotor efficiencies ---------------------------------------------

    lam = harness.truncnorm(opts.lambda_mu, opts.lambda_sigma, ...
                            opts.lambda_range, [opts.n 6]);

    % rc = zeros(opts.n,3);
    % for j = 1:3
    %     rc(:,j) = harness.truncnorm(opts.rc_mu(j), opts.rc_sigma(j), ...
    %                                 opts.rc_range(j,:), [opts.n 1]);
    % end
    % ==============================================================
    % Package CoM - fixed for Task 3.4
    % ==============================================================
    
    rc_nominal = [0.02; 0.02; -0.15];
    
    rc = repmat(rc_nominal.', opts.n, 1);

    % FINAL helix - this is intentionally NOT MC.ref from the old campaign.
    ref = traj.helix(opts.duration, ...
        radius        = opts.helix_radius, ...
        rate          = opts.helix_rate, ...
        climb         = opts.helix_climb, ...
        center        = opts.helix_center, ...
        start_p       = opts.helix_start_p, ...
        approach_time = opts.helix_approach);

    n  = opts.n;
    nc = numel(opts.mrac);

    MC.samples = struct( ...
        'mass', mass, ...
        'lambda', lam, ...
        'r_c', rc, ...
        'seed', opts.seed, ...
        'vary', 'mass + CoM + six independent rotor efficiencies');

    MC.configs = opts.mrac;
    MC.tags    = arrayfun(@tagFor, opts.mrac, 'UniformOutput', false);
    MC.ref     = ref;

    % Tracking metrics
    for f = {'rms_p','peak_p','ss_p','rms_v'}
        MC.(f{1}) = nan(n,nc);
    end

    % Robustness / actuator metrics
    MC.max_tilt_deg     = nan(n,nc);
    MC.sat_upper_frac   = nan(n,nc);
    MC.max_omega_eq     = nan(n,nc);
    MC.min_lambda_hat   = nan(n,nc);
    MC.proj_region_frac = nan(n,nc);

    MC.failed    = false(n,nc);
    MC.completed = false(n,nc);
    MC.errors    = cell(n,nc);
    MC.run_s     = nan(n,nc);
    MC.elapsed   = 0;

    file = fullfile(opts.outdir, sprintf('%s_%s.mat', opts.label, ...
                    char(datetime('now','Format','yyyyMMdd_HHmmss'))));

    % Store fixed plant/controller constants used by the report.
    Pmeta = harness.config(package_mass = mass(1), package_r_c = rc(1,:).');
    uavMeta = struct( ...
        'Omega_min', Pmeta.UAV.Omega_min, ...
        'Omega_max', Pmeta.UAV.Omega_max, ...
        'k_f',       Pmeta.UAV.k_f);

    controllerMeta = struct( ...
        'L_p', opts.L_p, ...
        'L_v', opts.L_v, ...
        'Gamma_theta', opts.Gamma_theta, ...
        'Gamma_lambda', opts.Gamma_lambda, ...
        'sigma_e', opts.sigma_e, ...
        'lambda_proj_center', opts.lambda_proj_center, ...
        'lambda_proj_r0', opts.lambda_proj_r0, ...
        'lambda_proj_r1', opts.lambda_proj_r1);

    helixMeta = struct( ...
        'duration', opts.duration, ...
        'radius', opts.helix_radius, ...
        'rate', opts.helix_rate, ...
        'climb', opts.helix_climb, ...
        'center', opts.helix_center, ...
        'entry', 'searched (free-entry minimum snap)', ...
        'start_p', opts.helix_start_p, ...
        'approach_time', opts.helix_approach);

    MC.meta = struct( ...
        'file', file, ...
        'timestamp', datetime('now'), ...
        'runner', 'run_montecarlo_final', ...
        'model', opts.model, ...
        'matlab', version('-release'), ...
        'n', n, ...
        'configs', opts.mrac, ...
        'parallel', opts.parallel, ...
        'seed', opts.seed, ...
        'vary', 'all', ...
        'mass_dist', [opts.mass_mu opts.mass_sigma opts.mass_range], ...
        'lambda_dist', [opts.lambda_mu opts.lambda_sigma opts.lambda_range], ...
        'rc_dist', struct('mu',opts.rc_mu,'sigma',opts.rc_sigma,'range',opts.rc_range), ...
        'helix', helixMeta, ...
        'controller', controllerMeta, ...
        'uav', uavMeta, ...
        'elapsed_s', 0, ...
        'settings', opts);

else
    loaded = load(opts.resume,'MC');
    MC = loaded.MC;

    if ~isfield(MC,'meta') || ~isfield(MC.meta,'runner') || ...
            ~strcmp(MC.meta.runner,'run_montecarlo_final')
        error('run_montecarlo_final:badResume', ...
              'Resume file is not a checkpoint produced by run_montecarlo_final.');
    end

    % A checkpoint owns the settings with which it was created.
    % How to run (not what to run) is taken from this call, so a campaign can
    % be resumed serially or in parallel regardless of how it was started.
    runMode = struct('parallel',opts.parallel,'workers',opts.workers, ...
                     'batch',opts.batch,'resume',opts.resume);
    opts = MC.meta.settings;
    for f = fieldnames(runMode).'
        opts.(f{1}) = runMode.(f{1});
    end

    file = MC.meta.file;
    ref  = MC.ref;
    mass = MC.samples.mass;
    lam  = MC.samples.lambda;
    rc   = MC.samples.r_c;

    n  = size(MC.rms_p,1);
    nc = size(MC.rms_p,2);

    % Defensive support for a checkpoint made by an early final-runner draft.
    MC = ensureMetricFields(MC,n,nc);

    fprintf('Resuming %s: %d of %d runs already done.\n', ...
        opts.resume, nnz(MC.completed), numel(MC.completed));
end

%% ========================================================================
%  Prepare model: keep only signals required by the scalar metrics
%  ========================================================================
keepLogs = {'p','p_d','v','v_d','T_out','euler','lambda_hat'};

% Serial: guard the model in this session. Parallel: parsim needs the
% model unmodified here; every worker applies the same guard itself
% (harness.mc_worker_setup).
if ~opts.parallel
    guard = harness.mc_guard_final(opts.model, keep=keepLogs); %#ok<NASGU>
end
closer = onCleanup(@() closeQuietly(opts.model));          %#ok<NASGU>

fprintf('\n============================================================\n');
fprintf(' FINAL TASK 3.4 MONTE CARLO\n');
fprintf('============================================================\n');
fprintf('%d samples x %d configuration(s) = %d simulations\n', ...
    n,nc,n*nc);
fprintf('duration     : %.1f s\n',ref.time(end));
fprintf('trajectory   : R=%.2f m, rate=%.3f rad/s, climb=%.2f m/s, approach=%.1f s\n', ...
    opts.helix_radius,opts.helix_rate,opts.helix_climb,opts.helix_approach);
fprintf('mass         : U(%.2f, %.2f) kg\n', ...
        opts.mass_range(1), opts.mass_range(2));
fprintf('lambda       : N(%.2f,%.2f^2), truncated [%.2f %.2f], six independent draws\n', ...
    opts.lambda_mu,opts.lambda_sigma,opts.lambda_range);
fprintf('L_p,L_v      : %.1f, %.1f\n',opts.L_p,opts.L_v);
fprintf('Gamma_theta  : diag([%.3f %.3f %.3f])\n',diag(opts.Gamma_theta));
fprintf('Gamma_lambda : [%.3f %.3f %.3f]^T\n',opts.Gamma_lambda);
fprintf('sigma_e      : %.4f  [e-modification on theta]\n',opts.sigma_e);
fprintf('checkpoint   : %s\n',file);
fprintf('solver       : untouched\n\n');

checkpoint(file,MC);
total = numel(MC.completed);

% ========================================================================
%  Run paired baseline / PBMRAC campaign
%  ========================================================================
UAV = P_nominal_uav();
if opts.parallel
    MC = runParallel(MC, ref, lam, mass, rc, opts, keepLogs, UAV, file, total);
else
    for c = 1:nc
        for i = 1:n
            if MC.completed(i,c)
                continue
            end
            assertDisk(opts.outdir,opts.min_free_gb);
            t0 = tic;
            try
                in = buildInput(ref, lam, mass, rc, i, MC.configs(c), opts);
                R  = harness.unpack(sim(in));
                [MC, status] = recordRun(MC, R, i, c, lam, opts, UAV, ref);
            catch ME
                MC.failed(i,c) = true;
                MC.errors{i,c} = sprintf('%s: %s',ME.identifier,ME.message);
                status = ['FAILED  ' ME.message];
            end

            % No dense time history survives the iteration.
            clear R in

            try
                Simulink.sdi.clear;
            catch
            end

            MC.completed(i,c) = true;
            MC.run_s(i,c)     = toc(t0);
            MC.elapsed        = MC.elapsed + MC.run_s(i,c);
            MC.meta.elapsed_s = MC.elapsed;

            checkpoint(file,MC);

            done = nnz(MC.completed);
            fprintf('%4d/%-4d %-9s %-58s %5.1fs  eta %s\n', ...
                done,total,MC.tags{c},status,MC.run_s(i,c), ...
                etaString(MC.elapsed,done,total));
        end
    end
end

% ========================================================================
%  Final summary and plots
%  ========================================================================
harness.mc_summary_final(MC);

if any(MC.failed(:))
    fprintf('\n%d of %d simulations failed and are excluded (see MC.errors).\n', ...
        nnz(MC.failed),total);
end

fprintf('\nSaved to %s\n',file);

if opts.plots
    harness.mc_report_final(MC,export=opts.export);
else
    fprintf('Plot later with:\n');
    fprintf('   load(''%s''); harness.mc_report_final(MC);\n',file);
end
end


% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function tag = tagFor(c)
if c
    tag = 'PBMRAC';
else
    tag = 'Baseline';
end
end

function validateRanges(opts)
if opts.mass_range(1) >= opts.mass_range(2)
    error('run_montecarlo_final:badMassRange','mass_range must satisfy lower < upper.');
end
if opts.lambda_range(1) >= opts.lambda_range(2)
    error('run_montecarlo_final:badLambdaRange','lambda_range must satisfy lower < upper.');
end
if any(opts.rc_range(:,1) >= opts.rc_range(:,2))
    error('run_montecarlo_final:badRcRange', ...
          'Each row of rc_range must satisfy lower < upper.');
end
if opts.mass_mu < opts.mass_range(1) || opts.mass_mu > opts.mass_range(2)
    error('run_montecarlo_final:badMassMean','mass_mu must lie inside mass_range.');
end
if opts.lambda_mu < opts.lambda_range(1) || opts.lambda_mu > opts.lambda_range(2)
    error('run_montecarlo_final:badLambdaMean','lambda_mu must lie inside lambda_range.');
end
if any(opts.rc_mu < opts.rc_range(:,1) | opts.rc_mu > opts.rc_range(:,2))
    error('run_montecarlo_final:badRcMean', ...
          'Each rc_mu component must lie inside its rc_range.');
end
end

function MC = ensureMetricFields(MC,n,nc)
fields = {'max_tilt_deg','sat_upper_frac','max_omega_eq', ...
          'min_lambda_hat','proj_region_frac'};
for k = 1:numel(fields)
    if ~isfield(MC,fields{k})
        MC.(fields{k}) = nan(n,nc);
    end
end
end

function checkpoint(file,MC)
% Save to a .pending file, then move it over the checkpoint. A sync client
% (OneDrive) can hold the checkpoint for a moment while it uploads it; the
% move is then retried for up to ~30 s instead of aborting the campaign. If
% the file is still locked, the campaign continues: the .pending file holds
% the latest state and the next checkpoint tries again.
pending = [file '.pending'];
save(pending,'MC','-v7');
for attempt = 1:30
    [ok,msg] = movefile(pending,file,'f');
    if ok, return; end
    pause(1);
end
warning('run_montecarlo_final:checkpoint', ...
        ['Could not install checkpoint (%s). Continuing; the latest state ' ...
         'is in %s.'],msg,pending);
end

function assertDisk(folder,minGb)
for loc = {folder,tempdir}
    free = double(java.io.File(loc{1}).getUsableSpace())/2^30;
    if free < minGb
        error('run_montecarlo_final:diskSpace', ...
              ['Only %.2f GiB free at %s. The checkpoint is intact; ' ...
               'free space and resume.'],free,loc{1});
    end
end
end

function s = etaString(elapsed,done,total)
if done == 0
    s = '--';
    return
end
remaining = elapsed/done*(total-done);
s = string(duration(0,0,remaining,'Format','hh:mm:ss'));
end

function closeQuietly(model)
try
    if bdIsLoaded(model)
        close_system(model,0);
    end
catch
end
end


% ========================================================================
%  One run: input, metrics (shared by the serial and the parallel path)
% ========================================================================

function UAV = P_nominal_uav()
P = harness.config();
UAV = P.UAV;
end

function in = buildInput(ref, lam, mass, rc, i, mrac, opts)
P = harness.config( ...
    package_mass = mass(i), ...
    package_r_c  = rc(i,:).');

% Force the final adaptive tuning explicitly. This removes any
% ambiguity if design.adaptive defaults are edited later.
P.adapt.L_p = opts.L_p;
P.adapt.L_v = opts.L_v;
P.adapt.Gamma_theta  = opts.Gamma_theta;
P.adapt.Gamma_lambda = opts.Gamma_lambda;
P.adapt.sigma_e      = opts.sigma_e;

P.adapt.lambda_proj_center = opts.lambda_proj_center;
P.adapt.lambda_proj_r0     = opts.lambda_proj_r0;
P.adapt.lambda_proj_r1     = opts.lambda_proj_r1;

lamSig = disturb.fault(ref.time, t0 = 0, lambda = lam(i,:).');
robustFlag = opts.robust_adapt * double(mrac ~= 0);

in = harness.input( ...
    ref, ...
    config  = P, ...
    mrac    = mrac, ...
    robust  = robustFlag, ...
    package = opts.package, ...
    lambda  = lamSig, ...
    model   = opts.model);
in = in.setModelParameter('TimeOut',opts.timeout);
end

function [MC, status] = recordRun(MC, R, i, c, lam, opts, UAV, ref)
if isempty(R.t) || R.t(end) < 0.999*ref.time(end)
    error('run_montecarlo_final:incomplete', ...
          'Stopped at %.4g s, expected %.4g s.', R.t(end), ref.time(end));
end

% Standard Task-3.4 metrics
M = harness.metrics(R);
for f = {'rms_p','peak_p','ss_p','rms_v'}
    if ~isfinite(M.(f{1}))
        error('run_montecarlo_final:nonfinite', 'Metric %s is not finite.', f{1});
    end
    MC.(f{1})(i,c) = M.(f{1});
end

% Euler log is ZYX = [yaw pitch roll].
if ~isfield(R,'euler') || size(R.euler,2) < 3
    error('run_montecarlo_final:missingEuler', 'euler log is required for max_tilt_deg.');
end
rp = rad2deg(R.euler(:,2:3));
MC.max_tilt_deg(i,c) = max(abs(rp),[],'all');

% Physical upper actuator-stop occupancy.
if ~isfield(R,'T_out') || size(R.T_out,2) ~= 6
    error('run_montecarlo_final:missingTout', 'T_out log is required for actuator robustness metrics.');
end
lamRow  = reshape(lam(i,:),1,[]);
Tpre    = R.T_out ./ lamRow;
Omega   = sqrt(abs(Tpre)/UAV.k_f);
atUpper = any(Omega >= 0.999*UAV.Omega_max,2);
MC.sat_upper_frac(i,c) = mean(atUpper);
MC.max_omega_eq(i,c)   = max(Omega,[],'all');

% Adaptive-state robustness. For baseline these are kept as NaN.
if MC.configs(c) ~= 0
    if ~isfield(R,'lambda_hat') || isempty(R.lambda_hat)
        error('run_montecarlo_final:missingLambdaHat', ...
              'lambda_hat log is required for PBMRAC robustness metrics.');
    end
    MC.min_lambda_hat(i,c) = min(R.lambda_hat,[],'all');
    inProj = any(R.lambda_hat < (opts.lambda_proj_center - opts.lambda_proj_r0),2);
    MC.proj_region_frac(i,c) = mean(inProj);
end

status = sprintf('RMS %.5f m | tilt %.1f deg | sat %.2f%%', ...
    MC.rms_p(i,c), MC.max_tilt_deg(i,c), 100*MC.sat_upper_frac(i,c));
end

% ========================================================================
%  Parallel path: parsim in batches, checkpoint after every batch
% ========================================================================

function MC = runParallel(MC, ref, lam, mass, rc, opts, keepLogs, UAV, file, total)
pool = gcp('nocreate');
if isempty(pool)
    if opts.workers > 0
        pool = parpool('local', opts.workers);
    else
        pool = parpool('local');
    end
elseif opts.workers > 0 && pool.NumWorkers ~= opts.workers
    fprintf('Using the existing pool (%d workers), not %d.\n', pool.NumWorkers, opts.workers);
end
B = opts.batch;
if B == 0
    B = 4*pool.NumWorkers;
end
fprintf('parallel     : %d workers, %d runs per batch\n\n', pool.NumWorkers, B);

% Each worker loads the model once and applies the same logging guard as
% the serial path (scopes and unused loggers commented out, in memory only).
setup = @() harness.mc_worker_setup(opts.model, keepLogs);

nc = size(MC.completed,2);
for c = 1:nc
    todo = find(~MC.completed(:,c));
    for b0 = 1:B:numel(todo)
        ids = todo(b0:min(numel(todo), b0+B-1));
        assertDisk(opts.outdir,opts.min_free_gb);

        clear ins
        ins(numel(ids)) = Simulink.SimulationInput(opts.model); %#ok<AGROW>
        for k = 1:numel(ids)
            ins(k) = buildInput(ref, lam, mass, rc, ids(k), MC.configs(c), opts);
        end

        t0 = tic;
        outs = parsim(ins, ...
            'ShowProgress', 'off', ...
            'ShowSimulationManager', 'off', ...
            'TransferBaseWorkspaceVariables', 'off', ...
            'SetupFcn', setup, ...
            'StopOnError', 'off');
        wall = toc(t0);

        nFail = 0;
        for k = 1:numel(ids)
            i = ids(k);
            try
                if ~isempty(outs(k).ErrorMessage)
                    error('run_montecarlo_final:sim', '%s', outs(k).ErrorMessage);
                end
                R = harness.unpack(outs(k));
                MC = recordRun(MC, R, i, c, lam, opts, UAV, ref);
            catch ME
                MC.failed(i,c) = true;
                MC.errors{i,c} = sprintf('%s: %s',ME.identifier,ME.message);
                nFail = nFail + 1;
            end
            MC.completed(i,c) = true;
            MC.run_s(i,c)     = wall/numel(ids);
        end
        clear outs ins R

        MC.elapsed        = MC.elapsed + wall;
        MC.meta.elapsed_s = MC.elapsed;
        checkpoint(file,MC);

        done = nnz(MC.completed);
        okIds = ids(~MC.failed(ids,c));
        fprintf('%4d/%-4d %-9s batch of %3d: median RMS %.4f m, max %.3f m, %d failed   %6.1fs  eta %s\n', ...
            done, total, MC.tags{c}, numel(ids), median(MC.rms_p(okIds,c)), max(MC.rms_p(okIds,c)), nFail, wall, ...
            etaString(MC.elapsed, done, total));
    end
end
end
