function guard = mc_guard_oldhelix(model, opts)
%MC_GUARD  Strip a model down to the logging a campaign needs, and restore it.
%
%   guard = harness.mc_guard_oldhelix('Final');
%   ... hundreds of sim() calls ...
%   clear guard          % everything is put back
%
%   WHY A LONG CAMPAIGN FAILS WITHOUT THIS
%   Final.slx carries 32 Scope blocks and 18 To-Workspace blocks, and the
%   Simulation Data Inspector defaults to storage mode 'disk'. Every run
%   therefore writes a .dmr file and fills scope buffers, none of which the
%   campaign reads - it keeps four scalars per run. Over several hundred runs
%   that exhausts memory (a crash inside webscopes_foundation, or "Data logging
%   exceeded available memory") and consumes a great deal of disk.
%
%   WHAT IS CHANGED - AND WHY IT CANNOT ALTER THE RESULTS
%     * SDI storage mode -> 'memory', so nothing is spooled to disk;
%     * every Scope block commented out;
%     * every To-Workspace block except the requested ones commented out.
%
%   All three are OUTPUT concerns. Scopes and loggers are sinks: they observe
%   the signals, they do not feed anything back into the model, so commenting
%   them cannot change a single state derivative. The SOLVER IS DELIBERATELY
%   NOT TOUCHED, which is the difference between this and a faster-but-different
%   campaign: results stay bit-comparable with the task scripts.
%
%   NOTHING IS SAVED. The block changes live only in the loaded model, and the
%   returned onCleanup restores every original Commented state - including if
%   the campaign errors part-way. Even an outright MATLAB crash leaves the .slx
%   on disk untouched, because it is never written.

arguments
    model     (1,:) char = 'Final'
    opts.keep cell = {'p','p_d','v','v_d'}   % the signals harness.unpack needs
    opts.quiet (1,1) logical = false
end

if ~bdIsLoaded(model)
    load_system(model);
end

state = struct('block', {}, 'commented', {});

%% ---- Scopes ----------------------------------------------------------
scopes = find_system(model, 'LookUnderMasks','all', 'FollowLinks','on', ...
                            'BlockType','Scope');
for i = 1:numel(scopes)
    state(end+1) = snapshot(scopes{i});             %#ok<AGROW>
    set_param(scopes{i}, 'Commented', 'on');
end

%% ---- To-Workspace blocks --------------------------------------------
logs  = find_system(model, 'LookUnderMasks','all', 'FollowLinks','on', ...
                           'BlockType','ToWorkspace');
kept  = {};
for i = 1:numel(logs)
    name = get_param(logs{i}, 'VariableName');
    state(end+1) = snapshot(logs{i});                %#ok<AGROW>
    if ismember(name, opts.keep)
        set_param(logs{i}, 'Commented', 'off');
        kept{end+1} = name;                          %#ok<AGROW>
    else
        set_param(logs{i}, 'Commented', 'on');
    end
end

missing = setdiff(opts.keep, kept);
if ~isempty(missing)
    % Put the model back before complaining, so a bad 'keep' list does not
    % leave the model in a stripped state.
    restore(state, []);
    error('harness:mc_guard_oldhelix:missingLogger', ...
          'No To-Workspace block logs: %s', strjoin(missing, ', '));
end

%% ---- Simulation Data Inspector --------------------------------------
% Storage mode lives outside the model, so it is captured and restored the
% same way the block states are. An empty value means the SDI is unavailable,
% and restore() then leaves it alone.
prevMode = currentStorageMode();
if ~isempty(prevMode)
    if sdiRunCount() == 0
        % Changing storage mode discards recorded runs, so only switch when
        % there is nothing to discard.
        quietly(@() Simulink.sdi.setStorageMode('memory'));
    elseif ~opts.quiet
        fprintf(['  note: the SDI already holds %d run(s), so storage mode stays ' ...
                 '''%s''.\n        Start a fresh session to get the disk saving.\n'], ...
                sdiRunCount(), prevMode);
    end
end
quietly(@() Simulink.sdi.setAutoArchiveMode(false));
quietly(@() Simulink.sdi.clear);

if ~opts.quiet
    fprintf('mc_guard: %d scopes and %d of %d loggers disabled; keeping %s.\n', ...
            numel(scopes), numel(logs)-numel(kept), numel(logs), strjoin(kept, ', '));
end

guard = onCleanup(@() restore(state, prevMode));
end

% ---------------------------------------------------------------------
function m = currentStorageMode()
%CURRENTSTORAGEMODE  '' when the Simulation Data Inspector is not available.
m = '';
try
    m = Simulink.sdi.getStorageMode();
catch
end
end

function n = sdiRunCount()
n = 0;
try
    n = Simulink.sdi.getRunCount();
catch
end
end

function quietly(fn)
%QUIETLY  Best-effort SDI call: absent on installations without the inspector.
try
    fn();
catch
end
end

function s = snapshot(block)
s = struct('block', block, 'commented', get_param(block, 'Commented'));
end

function restore(state, prevMode)
for i = 1:numel(state)
    try
        set_param(state(i).block, 'Commented', state(i).commented);
    catch
        % Model already closed - nothing to restore into.
    end
end
if ~isempty(prevMode)
    quietly(@() Simulink.sdi.setStorageMode(prevMode));
end
end
