function T = mc_summary_oldhelix(MC)
%MC_SUMMARY  Print the Monte Carlo statistics table and return it.
%
%   T = harness.mc_summary_oldhelix(MC)
%
%   Shared by harness.montecarlo (at the end of a campaign) and by
%   run_task3_montecarlo.m (after reloading a saved .mat), so a reloaded
%   campaign reports exactly the same numbers in exactly the same layout as a
%   freshly run one.
%
%   The p95 column is the number that belongs in the report alongside the mean:
%   the mean says how the controller behaves typically, the 95th percentile says
%   how it behaves in the unlucky corner of the parameter space, which is what a
%   delivery vehicle actually has to survive.

arguments
    MC struct
end

nc   = numel(MC.tags);
rows = cell(nc,1);

fprintf('\n%-10s %8s %10s %10s %10s %10s %8s\n', ...
        'case','runs','mean RMS','median','p95','max','failed');
fprintf('%s\n', repmat('-',1,70));

for c = 1:nc
    x  = MC.rms_p(:,c);
    x  = x(~isnan(x));
    nf = 0;
    if isfield(MC,'failed'), nf = nnz(MC.failed(:,c)); end
    fprintf('%-10s %8d %10.4f %10.4f %10.4f %10.4f %8d\n', MC.tags{c}, ...
            numel(x), mean(x), median(x), harness.percentile_oldhelix(x,95), max(x), nf);
    rows{c} = {MC.tags{c}, numel(x), mean(x), median(x), ...
               harness.percentile_oldhelix(x,95), max(x), nf};
end

if nc == 2
    a  = MC.rms_p(:,1);
    b  = MC.rms_p(:,2);
    ok = ~isnan(a) & ~isnan(b);
    fprintf('\n%s improves mean RMS by %.1f%% and p95 by %.1f%% over the campaign.\n', ...
        MC.tags{2}, ...
        100*(1 - mean(b(ok))/mean(a(ok))), ...
        100*(1 - harness.percentile_oldhelix(b(ok),95)/harness.percentile_oldhelix(a(ok),95)));
end

T = cell2table(vertcat(rows{:}), 'VariableNames', ...
               {'case','runs','mean_rms','median_rms','p95_rms','max_rms','failed'});
end
