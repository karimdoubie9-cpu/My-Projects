function S = step_metrics(R, ax, opts)
%STEP_METRICS  Rise time, overshoot and settling time of a setpoint response.
%
%   S = harness.step_metrics(R, 1)                       % x axis
%   S = harness.step_metrics(R, 3, ts_target=15, os_target=35)
%
%   These are the numbers Task 2.2 states its objectives in, so they are what
%   the tuning has to be judged against on the FULL NONLINEAR model - not just
%   on the linear design model in run_task2_frequency.m.
%
%   DEFINITIONS (MATLAB stepinfo conventions)
%     rise time     10% to 90% of the commanded step
%     overshoot     100 * (peak - final) / step, as a percentage of the STEP,
%                   so it is comparable between a 1 m and a 5 m command
%     settling time measured from the step instant to the last time the
%                   response leaves a +-2% band around the commanded value
%     ss error      mean of the last 10% of the run, minus the command
%
%   NORMALISATION. Two corrections are applied before calling stepinfo, and
%   both matter here:
%
%     1. The response is shifted and scaled to go 0 -> 1. The vertical step
%        runs 1 m -> 2 m, and feeding that in raw would have stepinfo measure
%        overshoot against a "final value" of 2 with a start of 1, inflating
%        every percentage.
%     2. Time is re-zeroed at the step instant, and samples before it are
%        dropped. Otherwise the settling time silently includes however long
%        the reference sat at its initial value.
%
%   The step instant is detected from the REFERENCE (p_d), not assumed, so
%   this works for whatever traj.setpoint produced.
%
%   Returns a struct with the metrics plus pass/fail against any targets
%   given, and prints a table unless quiet=true.

arguments
    R              struct
    ax             (1,1) double {mustBeMember(ax,[1 2 3])} = 1
    opts.threshold (1,1) double = 0.02      % settling band, fraction of step
    opts.rise      (1,2) double = [0.1 0.9]
    opts.os_target (1,1) double = NaN       % [%] objective, for the verdict
    opts.ts_target (1,1) double = NaN       % [s] objective
    opts.label     (1,:) char   = ''
    opts.quiet     (1,1) logical = false
end

axname = 'xyz';
t   = R.t(:);
y   = R.p(:,ax);
yd  = R.p_d(:,ax);

%% ---- Locate the step in the reference --------------------------------
y0   = yd(1);
yf   = yd(end);
step = yf - y0;

if abs(step) < 1e-9
    error('harness:step_metrics:noStep', ...
          'The reference does not step on the %s axis (start = end = %g).', ...
          axname(ax), y0);
end

% First sample where the reference has moved more than 1% of the step
k_step = find(abs(yd - y0) > 0.01*abs(step), 1, 'first');
if isempty(k_step), k_step = 1; end
t_step = t(k_step);

%% ---- Normalise and hand to stepinfo ----------------------------------
m  = t >= t_step;
tt = t(m) - t_step;
yy = (y(m) - y0) / step;                 % 0 -> 1 regardless of axis or size

si = stepinfo(yy, tt, 1, ...
              'SettlingTimeThreshold', opts.threshold, ...
              'RiseTimeLimits',        opts.rise);

%% ---- Assemble (back in physical units where meaningful) --------------
n_tail = max(1, round(0.10*numel(yy)));

S.axis          = axname(ax);
S.t_step        = t_step;
S.from          = y0;
S.to            = yf;
S.step          = step;
S.RiseTime      = si.RiseTime;
S.SettlingTime  = si.SettlingTime;
S.Overshoot     = si.Overshoot;              % [%] of the step
S.Undershoot    = si.Undershoot;
S.PeakValue     = y0 + si.Peak*step;         % [m] absolute
S.PeakTime      = si.PeakTime;
S.SSError       = mean(y(end-n_tail+1:end)) - yf;
S.threshold     = opts.threshold;
S.label         = opts.label;

S.os_target = opts.os_target;
S.ts_target = opts.ts_target;
S.os_pass   = isnan(opts.os_target) || S.Overshoot    <= opts.os_target;
S.ts_pass   = isnan(opts.ts_target) || S.SettlingTime <= opts.ts_target;

if opts.quiet, return, end

%% ---- Report ----------------------------------------------------------
name = opts.label;
if isempty(name), name = sprintf('%s-axis step', axname(ax)); end

fprintf('\n=== Step response metrics: %s ===\n', name);
fprintf('  command        : %.3f -> %.3f  (step %.3f m) at t = %.3f s\n', ...
        y0, yf, step, t_step);
fprintf('  rise time      : %7.3f s      (%.0f%% to %.0f%%)\n', ...
        S.RiseTime, 100*opts.rise(1), 100*opts.rise(2));
fprintf('  peak           : %7.3f m      at t = %.3f s\n', S.PeakValue, S.PeakTime);
fprintf('  overshoot      : %7.2f %%      %s\n', S.Overshoot, target(S.os_target, S.os_pass, '%'));
fprintf('  settling time  : %7.3f s      %s   (+-%.0f%% band)\n', ...
        S.SettlingTime, target(S.ts_target, S.ts_pass, 's'), 100*opts.threshold);
fprintf('  steady error   : %7.4f m\n', S.SSError);

if ~isnan(opts.os_target) || ~isnan(opts.ts_target)
    if S.os_pass && S.ts_pass
        fprintf('  VERDICT        : objectives met\n');
    else
        fprintf('  VERDICT        : OBJECTIVES NOT MET -> retune (Task 2.3.1)\n');
    end
end
end

% ---------------------------------------------------------------------
function s = target(val, pass, unit)
%TARGET  Render the objective and whether it was met, or nothing if unset.
if isnan(val)
    s = '';
else
    s = sprintf('(target <= %g %s : %s)', val, unit, verdict(pass));
end
end

function s = verdict(pass)
if pass, s = 'PASS'; else, s = 'FAIL'; end
end
