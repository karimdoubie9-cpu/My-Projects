function guard = mc_guard_final(model, opts)
%MC_GUARD_FINAL  Memory-safe logging guard for the final Monte Carlo campaign.
%
%   guard = harness.mc_guard_final('Final')
%   guard = harness.mc_guard_final('Final', ...
%       keep={'p','p_d','v','v_d','T_out','euler','lambda_hat'})
%
% Only output sinks are changed:
%   * Scope blocks are commented out;
%   * all To-Workspace blocks except opts.keep are commented out;
%   * SDI storage is kept in memory and cleared between runs.
%
% The solver and all dynamic/control blocks are untouched. The onCleanup object
% restores the original model state when the campaign exits, including on error.
% The .slx file is never saved by this helper.

arguments
    model (1,:) char = 'Final'
    opts.keep cell = {'p','p_d','v','v_d','T_out','euler','lambda_hat'}
    opts.quiet (1,1) logical = false
end

if ~bdIsLoaded(model)
    load_system(model);
end

state = struct('block',{},'commented',{});

%% Scopes
scopes = find_system(model, ...
    'LookUnderMasks','all', ...
    'FollowLinks','on', ...
    'BlockType','Scope');

for i = 1:numel(scopes)
    state(end+1) = snapshot(scopes{i}); %#ok<AGROW>
    set_param(scopes{i},'Commented','on');
end

%% To-Workspace blocks
logs = find_system(model, ...
    'LookUnderMasks','all', ...
    'FollowLinks','on', ...
    'BlockType','ToWorkspace');

kept = {};

for i = 1:numel(logs)
    name = get_param(logs{i},'VariableName');

    state(end+1) = snapshot(logs{i}); %#ok<AGROW>

    if ismember(name,opts.keep)
        set_param(logs{i},'Commented','off');
        kept{end+1} = name; %#ok<AGROW>
    else
        set_param(logs{i},'Commented','on');
    end
end

missing = setdiff(opts.keep,kept);

if ~isempty(missing)
    restore(state,[]);
    error('harness:mc_guard_final:missingLogger', ...
          'No To-Workspace block logs: %s',strjoin(missing,', '));
end

%% Simulation Data Inspector
prevMode = currentStorageMode();

if ~isempty(prevMode)
    if sdiRunCount() == 0
        quietly(@() Simulink.sdi.setStorageMode('memory'));
    elseif ~opts.quiet
        fprintf(['  note: SDI already holds %d run(s), so storage mode stays ' ...
                 '''%s''.\n'],sdiRunCount(),prevMode);
    end
end

quietly(@() Simulink.sdi.setAutoArchiveMode(false));
quietly(@() Simulink.sdi.clear);

if ~opts.quiet
    fprintf(['mc_guard_final: %d scopes and %d of %d loggers disabled; ' ...
             'keeping %s.\n'], ...
        numel(scopes), ...
        numel(logs)-numel(kept), ...
        numel(logs), ...
        strjoin(kept,', '));
end

guard = onCleanup(@() restore(state,prevMode));
end


function m = currentStorageMode()
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
try
    fn();
catch
end
end

function s = snapshot(block)
s = struct('block',block,'commented',get_param(block,'Commented'));
end

function restore(state,prevMode)
for i = 1:numel(state)
    try
        set_param(state(i).block,'Commented',state(i).commented);
    catch
    end
end

if ~isempty(prevMode)
    quietly(@() Simulink.sdi.setStorageMode(prevMode));
end
end
