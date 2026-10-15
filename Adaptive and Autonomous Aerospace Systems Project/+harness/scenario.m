function R = scenario(ref, opts)
%SCENARIO  Run one simulation of Final.slx and return unpacked results.
%
%   R = harness.scenario(ref)
%   R = harness.scenario(ref, mrac=1, package=1, wind=W, lambda=L)
%
%   Thin wrapper: harness.input -> sim -> harness.unpack.  All argument
%   handling lives in harness.input; see there for the option list.
%
%   R contains t, p, p_d, v, v_d (each N-by-3) plus whatever optional signals
%   the model logged, and R.label / R.opts describing the run.
%
%   FAST RESTART is not enabled by a single scenario call -- there is nothing to
%   amortise over one run.  For a BATCH of runs it is both safe and worthwhile:
%   wrap the batch in harness.fastrestart and get roughly a 2.3x speedup.
%
%       guard = harness.fastrestart('Final');
%       for i = 1:N, R(i) = harness.scenario(ref, ...); end
%       clear guard
%
%   This was measured, not assumed: sweeping the package mass, the effectiveness
%   signal, MRAC_on, the reference and the stop time all give bit-identical
%   results with fast restart on. See harness.fastrestart for the details.
%
%   See also harness.input, harness.fastrestart, harness.montecarlo.

arguments
    ref            struct
    opts.mrac      (1,1) double = 0
    opts.package   (1,1) double = 0
    opts.robust    (1,1) double {mustBeMember(opts.robust,[0 1])} = 0
    opts.wind        struct = struct([])
    opts.force       struct = struct([])
    opts.lambda      struct = struct([])
    opts.disturbance struct = struct([])
    opts.config      struct = struct([])
    opts.model     (1,:) char = 'Final'
    opts.stop_time (1,1) double = NaN
    opts.sync_body_velocity (1,1) logical = true
    opts.label     (1,:) char = ''
end

args = namedargs2cell(rmfield(opts, 'label'));
in   = harness.input(ref, args{:});

% Keep the Simulation Data Inspector from stacking up every run of a script:
% a dozen back-to-back simulations otherwise exhaust memory and fail with an
% error that points at logging rather than at the archive. See harness.sdi_reset.
harness.sdi_reset();

R = harness.unpack(sim(in));

if isempty(opts.label)
    R.label = sprintf('%s / %s', ternary(opts.mrac,'PBMRAC','Baseline'), ...
                                 ternary(opts.package,'package','no package'));
else
    R.label = opts.label;
end
R.opts = rmfield(opts, {'config','wind','force','lambda','disturbance'});
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
