function T = mc_summary_final(MC)
%MC_SUMMARY_FINAL  Final Task 3.4 Monte Carlo statistics.
%
% Prints the main RMS table plus paired-improvement and robustness summaries.

arguments
    MC struct
end

nc = numel(MC.tags);
rows = cell(nc,1);

fprintf('\n============================================================\n');
fprintf(' TASK 3.4 - MONTE CARLO TRACKING SUMMARY\n');
fprintf('============================================================\n');

fprintf('\n%-10s %8s %10s %10s %10s %10s %8s\n', ...
    'case','runs','mean RMS','median','p95','max','failed');
fprintf('%s\n',repmat('-',1,70));

for c = 1:nc
    x = MC.rms_p(:,c);
    x = x(~isnan(x));

    nf = 0;
    if isfield(MC,'failed')
        nf = nnz(MC.failed(:,c));
    end

    if isempty(x)
        mn = NaN; med = NaN; p95 = NaN; mx = NaN;
    else
        mn  = mean(x);
        med = median(x);
        p95 = harness.percentile_final(x,95);
        mx  = max(x);
    end

    fprintf('%-10s %8d %10.4f %10.4f %10.4f %10.4f %8d\n', ...
        MC.tags{c},numel(x),mn,med,p95,mx,nf);

    rows{c} = {MC.tags{c},numel(x),mn,med,p95,mx,nf};
end

if nc == 2
    a = MC.rms_p(:,1);
    b = MC.rms_p(:,2);
    ok = isfinite(a) & isfinite(b) & a > eps;

    imp = 100*(1 - b(ok)./a(ok));

    fprintf('\nPaired comparison over %d valid plants:\n',nnz(ok));

    if any(ok)
        fprintf('  PBMRAC improves RMS on        : %6.1f%% of samples\n', ...
            100*mean(b(ok) < a(ok)));
        fprintf('  mean paired RMS improvement  : %6.1f%%\n',mean(imp));
        fprintf('  median paired improvement    : %6.1f%%\n',median(imp));
        fprintf('  5th percentile improvement   : %6.1f%%\n', ...
            harness.percentile_final(imp,5));

        fprintf('\n  campaign mean RMS improvement: %6.1f%%\n', ...
            100*(1-mean(b(ok))/mean(a(ok))));
        fprintf('  campaign p95 RMS improvement : %6.1f%%\n', ...
            100*(1- ...
            harness.percentile_final(b(ok),95)/ ...
            harness.percentile_final(a(ok),95)));
    end
end

%% Robustness summary
if isfield(MC,'max_tilt_deg') && isfield(MC,'sat_upper_frac')
    fprintf('\n============================================================\n');
    fprintf(' ACTUATOR / ATTITUDE ROBUSTNESS\n');
    fprintf('============================================================\n');

    fprintf(['%-10s %8s %12s %12s %12s %12s %12s\n'], ...
        'case','runs','p95 tilt','max tilt','sat runs','mean sat','max sat');
    fprintf('%s\n',repmat('-',1,86));

    for c = 1:nc
        tilt = MC.max_tilt_deg(:,c);
        sat  = MC.sat_upper_frac(:,c);

        ok = isfinite(tilt) & isfinite(sat);

        if any(ok)
            p95tilt = harness.percentile_final(tilt(ok),95);
            maxtilt = max(tilt(ok));
            satRuns = 100*mean(sat(ok) > 1e-6);
            meanSat = 100*mean(sat(ok));
            maxSat  = 100*max(sat(ok));
        else
            p95tilt = NaN; maxtilt = NaN;
            satRuns = NaN; meanSat = NaN; maxSat = NaN;
        end

        fprintf('%-10s %8d %11.2f° %11.2f° %11.1f%% %11.3f%% %11.3f%%\n', ...
            MC.tags{c},nnz(ok),p95tilt,maxtilt,satRuns,meanSat,maxSat);
    end
end

%% PBMRAC Lambda robustness
if isfield(MC,'min_lambda_hat') && size(MC.min_lambda_hat,2) >= 2
    x = MC.min_lambda_hat(:,2);
    p = MC.proj_region_frac(:,2);
    ok = isfinite(x);

    if any(ok)
        fprintf('\nPBMRAC adaptive-effectiveness robustness:\n');
        fprintf('  minimum Lambda_hat over all valid runs : %.4f\n',min(x(ok)));
        fprintf('  p05 of per-run min Lambda_hat          : %.4f\n', ...
            harness.percentile_final(x(ok),5));

        okp = isfinite(p);
        if any(okp)
            fprintf('  runs entering Lambda projection region : %.1f%%\n', ...
                100*mean(p(okp) > 1e-6));
        end
    end
end

T = cell2table(vertcat(rows{:}), ...
    'VariableNames', ...
    {'case','runs','mean_rms','median_rms','p95_rms','max_rms','failed'});
end
