function files = export(task, opts)
%EXPORT  Save every open figure as a PNG under Figures/<task>/.
%
%   files = viz.export('Task4_mission')
%   files = viz.export('Task4_mission', dpi=300, clean=false)
%
%   Called once at the end of a task script. Every figure created by that script
%   is written to Figures/<task>/, named from the figure's own Name property, so
%   the file names stay descriptive without maintaining a second list of titles
%   that would drift out of step with the figures themselves.
%
%   NAMING
%     A figure named 'Mission speed profile' created fourth becomes
%
%         Figures/Task4_mission/04_Mission_speed_profile.png
%
%     The prefix is strictly CREATION ORDER, which is the order the script draws
%     them and therefore the order they should be read. Deriving it from the
%     figure name instead would let an unnumbered preamble figure and an
%     explicitly numbered one collide on the same index. Any leading '#N ' left
%     in a figure name is stripped, since the prefix now supersedes it; a
%     non-numeric tag such as '#A ' is kept, because it identifies the scenario
%     rather than the ordering. Unnamed figures fall back to
%     'figure_<handle number>'.
%
%   CLEANING
%     clean=true (the default) removes the PNGs already in that task's folder
%     before writing. Renaming a figure would otherwise leave the old file
%     behind, and a stale PNG that looks plausible is worse than a missing one
%     when the folder is the source for a report. Only *.png inside
%     Figures/<task>/ is touched -- never anything else, and never recursively.
%
%   exportgraphics is used rather than saveas/print: it crops to the axes
%   content and honours the figure background, so the PNGs drop straight into
%   slides without a wide white margin.

arguments
    task       (1,:) char
    opts.dpi   (1,1) double = 200
    opts.clean (1,1) logical = true
    opts.root  (1,:) char = 'Figures'
    opts.quiet (1,1) logical = false
end

folder = fullfile(opts.root, task);
if ~isfolder(folder)
    mkdir(folder);
end

if opts.clean
    stale = dir(fullfile(folder, '*.png'));
    for k = 1:numel(stale)
        delete(fullfile(folder, stale(k).name));
    end
    if ~isempty(stale) && ~opts.quiet
        fprintf('  (removed %d stale PNG(s) from %s)\n', numel(stale), folder);
    end
end

figs = findobj(groot, 'Type', 'figure');
figs = figs(isvalid(figs));            % a figure may have been closed already
if isempty(figs)
    if ~opts.quiet, fprintf('viz.export: no open figures to save.\n'); end
    files = {};
    return
end

% Sort into creation order. findobj returns newest-first, and Number is the
% creation counter -- but a figure made with 'IntegerHandle','off' has an EMPTY
% Number, so [figs.Number] would be shorter than figs and silently misalign the
% sort with the handles. Build the key element by element instead, and give the
% unnumbered ones keys that preserve their relative findobj order at the end.
nums = nan(numel(figs), 1);
for i = 1:numel(figs)
    n = figs(i).Number;
    if isnumeric(n) && isscalar(n), nums(i) = n; end
end
missing = isnan(nums);
if any(missing)
    base = max([nums(~missing); 0]);
    nums(missing) = base + flip(1:nnz(missing)).';   % undo newest-first
end
[~, order] = sort(nums);
figs = figs(order);

% tiledlayout titles and shared labels ignore the groot font defaults, so bring
% them in line with the axes font (viz.defaults) before anything is saved.
fnt = get(groot, 'defaultAxesFontName');
for i = 1:numel(figs)
    for tl = reshape(findall(figs(i), 'Type', 'tiledlayout'), 1, [])
        for part = {'Title', 'Subtitle', 'XLabel', 'YLabel'}
            try, tl.(part{1}).FontName = fnt; catch, end
        end
    end
end

files = cell(numel(figs), 1);

for k = 1:numel(figs)
    f = figs(k);
    % Re-check: plotting code with a timer or callback can close a window while
    % this loop runs, and reading Name off a deleted handle is a hard error.
    if ~isvalid(f)
        if ~opts.quiet
            fprintf('  (figure %d was closed before it could be saved)\n', k);
        end
        continue
    end
    name = f.Name;
    if isempty(name), name = sprintf('figure_%d', k); end

    % Drop a redundant numeric '#N ' tag; the creation-order prefix replaces it.
    % A non-numeric tag like '#A ' is scenario identification, so it survives
    % into the descriptive part of the name.
    name = regexprep(name, '^\s*#\d+\s*', '');

    files{k} = fullfile(folder, sprintf('%02d_%s.png', k, slug(name)));
    try
        exportgraphics(f, files{k}, Resolution = opts.dpi);
    catch ME
        % Name the figure, not just the path: knowing WHICH plot failed is what
        % makes the warning actionable.
        warning('viz:export:failed', 'Could not write "%s" -> %s (%s)', ...
                name, files{k}, ME.message);
        files{k} = '';
    end
end

files = files(~cellfun(@isempty, files));
if ~opts.quiet
    fprintf('\nSaved %d figure(s) to %s\n', numel(files), fullfile(pwd, folder));
end
end

% ---------------------------------------------------------------------
function s = slug(str)
%SLUG  Filename-safe token: keep letters, digits, dash and underscore.
s = char(str);
s = regexprep(s, '[^A-Za-z0-9]+', '_');   % collapse runs of anything else
s = regexprep(s, '_+', '_');
s = regexprep(s, '^_|_$', '');
if isempty(s), s = 'figure'; end
if numel(s) > 60, s = s(1:60); end        % keep paths comfortably short
end
