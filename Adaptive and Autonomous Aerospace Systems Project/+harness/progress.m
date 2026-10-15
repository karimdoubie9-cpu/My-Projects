function step = progress(total, opts)
%PROGRESS  Reusable progress reporter with elapsed time, ETA and rate.
%
%   step = harness.progress(400, label='Monte Carlo');
%   for k = 1:400
%       ... do work ...
%       step(k);
%   end
%   step();          % final line + summary
%
%   Returns a function handle. Call it with the number of items completed; call
%   it with no argument to close out the report.
%
%   MODES
%     'auto'    (default) in-place single line when running in the MATLAB
%               desktop, one line per update otherwise. Backspace overwriting
%               produces unreadable garbage in a diary or a piped log, so it is
%               only used where it actually renders.
%     'inplace' force the single overwriting line
%     'lines'   force one appended line per update (log friendly)
%     'none'    silent
%
%   The update interval adapts: it reports at most every `min_interval` seconds
%   and at least every `every` items, so a fast loop does not drown the console
%   and a slow one still shows signs of life.

arguments
    total             (1,1) double {mustBeInteger, mustBePositive}
    opts.label        (1,:) char = 'progress'
    opts.every        (1,1) double = 0     % 0 = auto (about 40 updates)
    opts.min_interval (1,1) double = 1.0   % [s] minimum time between updates
    opts.mode         (1,:) char {mustBeMember(opts.mode, ...
                        {'auto','inplace','lines','none'})} = 'auto'
end

mode = opts.mode;
if strcmp(mode, 'auto')
    if usejava('desktop'), mode = 'inplace'; else, mode = 'lines'; end
end

every = opts.every;
if every <= 0, every = max(1, round(total/40)); end

t0        = tic;
lastLen   = 0;
lastTime  = -inf;
lastCount = 0;
closed    = false;

step = @update;

% -----------------------------------------------------------------
    function update(k)
        if strcmp(mode,'none') || closed, return, end

        finishing = (nargin == 0);
        if finishing
            k = lastCount;
        else
            lastCount = k;
        end

        elapsed = toc(t0);
        due = finishing || k >= total || ...
              (mod(k, every) == 0 && (elapsed - lastTime) >= opts.min_interval);
        if ~due, return, end
        lastTime = elapsed;

        rate = elapsed / max(k, 1);                 % [s] per item
        eta  = rate * max(total - k, 0);
        txt  = sprintf('  %s  %d/%d  %3.0f%%  |  elapsed %s  eta %s  (%.2f s/run)', ...
                       opts.label, k, total, 100*k/total, ...
                       hms(elapsed), hms(eta), rate);

        switch mode
            case 'inplace'
                fprintf('%s%s', repmat(sprintf('\b'), 1, lastLen), txt);
                lastLen = numel(txt);
            case 'lines'
                fprintf('%s\n', txt);
        end

        if finishing || k >= total
            if strcmp(mode,'inplace'), fprintf('\n'); end
            fprintf('  %s complete: %d runs in %s (%.2f s/run)\n', ...
                    opts.label, k, hms(elapsed), rate);
            closed = true;
        end
    end
end

% ---------------------------------------------------------------------
function s = hms(sec)
%HMS  Compact duration: m:ss below an hour, h:mm:ss above it.
if ~isfinite(sec), s = '--'; return, end
sec = max(0, round(sec));
h = floor(sec/3600); m = floor(mod(sec,3600)/60); s2 = mod(sec,60);
if h > 0, s = sprintf('%d:%02d:%02d', h, m, s2);
else,     s = sprintf('%d:%02d', m, s2);
end
end
