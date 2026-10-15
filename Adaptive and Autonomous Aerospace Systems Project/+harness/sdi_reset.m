function sdi_reset()
%SDI_RESET  Stop the Simulation Data Inspector accumulating runs.
%
%   harness.sdi_reset()
%
%   WHY THIS IS NEEDED
%   Every sim() archives its logged signals into the Simulation Data Inspector,
%   and by default those runs are RETAINED for the whole session. The task
%   scripts here execute six to twelve simulations back to back, each carrying
%   ~10k samples across a dozen logged signals, so the archive grows until
%   MATLAB either refuses to log ("Data logging exceeded available memory") or
%   dies outright inside the scope buffer code. Neither failure points at the
%   real cause, which is why this is worth doing explicitly rather than leaving
%   to chance.
%
%   WHAT IT DOES
%     * turns OFF auto-archive, so completed runs are not stacked up;
%     * clears the runs already held.
%
%   It deliberately does NOT disable recording. Turning recording off does free
%   the memory, but it also empties the logged signals themselves - lambda_hat
%   comes back 0x0 and the Task 3 effectiveness figure silently has nothing to
%   plot. Clearing keeps the logging intact and only discards the archive.
%
%   Called before each run by harness.scenario, so the most recent run stays
%   inspectable in the SDI while older ones do not accumulate. Everything is
%   guarded: on installations without the SDI these calls simply do nothing.

try
    Simulink.sdi.setAutoArchiveMode(false);
catch
    % Older release, or SDI unavailable - nothing to do.
end

try
    Simulink.sdi.clear;
catch
    % Nothing archived yet, or SDI unavailable.
end
end
