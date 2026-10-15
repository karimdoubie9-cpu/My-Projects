function mc_worker_setup(model, keep)
%MC_WORKER_SETUP  parsim SetupFcn for run_montecarlo_final(parallel=true).
%
%   Runs on every pool worker: loads the model and applies the same logging
%   guard as the serial campaign (scopes and unused To-Workspace blocks are
%   commented out in memory; the .slx is never saved). The guard is kept in a
%   persistent variable so it stays active for every simulation the worker
%   runs, and it is applied only once per loaded model.
%
%   See also harness.mc_guard_final, run_montecarlo_final.

persistent guard
if ~isempty(guard) && bdIsLoaded(model)
    return                           % already guarded on this worker
end
guard = [];                          %#ok<NASGU>  release a stale guard first
load_system(model);
guard = harness.mc_guard_final(model, keep = keep, quiet = true);
end
