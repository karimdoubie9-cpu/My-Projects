function guard = fastrestart(model)
%FASTRESTART  Enable fast restart for a batch of runs, restoring it afterwards.
%
%   guard = harness.fastrestart('Final');
%   ... many sim() calls ...
%   clear guard          % or just let it go out of scope
%
%   Fast restart keeps the model compiled between successive sim() calls, which
%   is worth roughly a 2.3x speedup on this model (measured: 1.71 s/run -> 0.75
%   s/run on a 40 s helix).  The saving is pure compile/initialise overhead, so
%   it matters most for many SHORT runs -- exactly the Monte Carlo case.
%
%   MEASURED, NOT ASSUMED: with every variable delivered through a
%   Simulink.SimulationInput, changing the package mass (UAV_pkg.M, .Minv, .S),
%   the rotor effectiveness signal, MRAC_on, package_exist, the reference
%   trajectory AND the stop time between runs all give results bit-identical to
%   running with fast restart off.  There is no stale-parameter hazard here.
%
%   The returned onCleanup guard restores the previous setting even if the batch
%   errors out, so an interrupted sweep cannot leave the model stuck compiled
%   (which would then block edits in the Simulink editor).

arguments
    model (1,:) char = 'Final'
end

if ~bdIsLoaded(model)
    load_system(model);
end

previous = get_param(model, 'FastRestart');
guard    = onCleanup(@() restore(model, previous));
set_param(model, 'FastRestart', 'on');
end

function restore(model, previous)
try
    if bdIsLoaded(model)
        set_param(model, 'FastRestart', previous);
    end
catch
    % Model closed or already torn down: nothing to restore.
end
end
