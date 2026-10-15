function T = report(runs, tags, opts)
%REPORT  Print a comparison table of tracking metrics and return it.
%
%   T = harness.report({R1,R2}, {'Baseline','PBMRAC'})
%   T = harness.report(runs, tags, window=[60 80], title='Gust response')
%
%   Prints one row per run and, when exactly two runs are given, the relative
%   improvement of the second over the first -- the number that belongs in the
%   report next to every baseline-vs-adaptive figure.

arguments
    runs        cell
    tags        cell
    opts.window (1,2) double = [NaN NaN]
    opts.title  (1,:) char = 'Tracking metrics'
end

n = numel(runs);
M = cell(1,n);
for k = 1:n
    M{k} = harness.metrics(runs{k}, window=opts.window);
end

hasWin = ~any(isnan(opts.window));

fprintf('\n=== %s ===\n', opts.title);
if hasWin
    fprintf('%-22s %10s %10s %10s %10s %12s\n', ...
            'case','RMS |ep|','peak |ep|','ss |ep|','RMS |ev|','peak in win');
else
    fprintf('%-22s %10s %10s %10s %10s\n', ...
            'case','RMS |ep|','peak |ep|','ss |ep|','RMS |ev|');
end
fprintf('%s\n', repmat('-', 1, 78));

rows = cell(n,1);
for k = 1:n
    if hasWin
        fprintf('%-22s %10.4f %10.4f %10.4f %10.4f %12.4f\n', tags{k}, ...
                M{k}.rms_p, M{k}.peak_p, M{k}.ss_p, M{k}.rms_v, M{k}.window_peak_p);
        rows{k} = {tags{k}, M{k}.rms_p, M{k}.peak_p, M{k}.ss_p, M{k}.rms_v, M{k}.window_peak_p};
    else
        fprintf('%-22s %10.4f %10.4f %10.4f %10.4f\n', tags{k}, ...
                M{k}.rms_p, M{k}.peak_p, M{k}.ss_p, M{k}.rms_v);
        rows{k} = {tags{k}, M{k}.rms_p, M{k}.peak_p, M{k}.ss_p, M{k}.rms_v};
    end
end

if n == 2
    % Neutral wording: these two runs are not always "baseline vs improvement"
    % (nominal vs degraded is a comparison in the opposite direction), so state
    % the reduction of the second relative to the first and let the caller's
    % tags carry the meaning. Positive = the second run tracks better.
    fprintf('%s\n', repmat('-', 1, 78));
    fprintf('%s relative to %s:  RMS %+.1f%%   peak %+.1f%%   steady-state %+.1f%%\n', ...
        tags{2}, tags{1}, ...
        100*(1 - M{2}.rms_p /M{1}.rms_p ), ...
        100*(1 - M{2}.peak_p/M{1}.peak_p), ...
        100*(1 - M{2}.ss_p  /M{1}.ss_p  ));
    fprintf('(positive = %s tracks better)\n', tags{2});
end

vars = {'case','rms_p','peak_p','ss_p','rms_v'};
if hasWin, vars{end+1} = 'window_peak_p'; end
T = cell2table(vertcat(rows{:}), 'VariableNames', vars);
end
