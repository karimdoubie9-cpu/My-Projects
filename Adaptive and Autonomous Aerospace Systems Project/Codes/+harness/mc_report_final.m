function out = mc_report_final(MC, opts)
%MC_REPORT_FINAL  Final Task 3.4 plots and tail diagnostics.
%
% Figure order:
%   1. convergence
%   2. tracking metric distributions
%   3. sampled uncertainty distributions
%   4. sensitivity / spread attribution
%   5. paired RMS improvement
%   6. actuator / attitude / adaptive-state robustness
%
% Then prints the worst baseline and worst PBMRAC samples.

arguments
    MC struct
    opts.metric (1,:) char = 'rms_p'
    opts.tol    (1,1) double = 0.02
    opts.worst  (1,1) double {mustBeInteger,mustBePositive} = 5
    opts.export (1,:) char = ''
    opts.prefix (1,:) char = 'Task 3.4 final Monte Carlo'
end

%% 1. Convergence
[~,out.convergence] = viz.convergence( ...
    MC, ...
    metric=opts.metric, ...
    tol=opts.tol, ...
    name=sprintf('%s - convergence diagnostic',opts.prefix));

%% 2. Tracking metric distributions
viz.histograms( ...
    MC, ...
    metrics={'rms_p','peak_p','ss_p'}, ...
    name=sprintf('%s - tracking metric distributions',opts.prefix));

%% 3. Sampled uncertainties
rc = MC.samples.r_c;

figure( ...
    'Name',sprintf('%s - sampled uncertainties',opts.prefix), ...
    'Color','w', ...
    'Position',[120 120 1050 760]);

subplot(2,2,1)
histogram(MC.samples.mass,30,'EdgeColor','none');
grid on
xlabel('package mass [kg]');
ylabel('count');
title('Sampled payload mass');

subplot(2,2,2)
histogram(MC.samples.lambda(:),30,'EdgeColor','none');
grid on
xlabel('rotor effectiveness \lambda_r [-]');
ylabel('count');
title('Sampled propeller effectiveness');

subplot(2,2,3)
scatter(rc(:,1),rc(:,2),14,'filled');
grid on
axis equal
xlabel('r_{c,x} [m]');
ylabel('r_{c,y} [m]');
title('Sampled package CoM - horizontal plane');

subplot(2,2,4)
histogram(rc(:,3),30,'EdgeColor','none');
grid on
xlabel('r_{c,z} [m]');
ylabel('count');
title('Sampled package CoM - vertical position');

%% 4. Sensitivity
out.sensitivity = harness.sensitivity_final( ...
    MC,metric=opts.metric,quiet=true);

S = out.sensitivity;

if ~isempty(S)
    factorNames = S(1).names;
    nf = numel(factorNames);

    V = nan(numel(S),nf);
    for c = 1:numel(S)
        V(c,:) = 100*S(c).var_frac;
    end

    figure( ...
        'Name',sprintf('%s - variance attribution',opts.prefix), ...
        'Color','w', ...
        'Position',[130 130 950 450]);

    bar(V.','grouped');
    grid on

    set(gca, ...
        'XTick',1:nf, ...
        'XTickLabel',factorNames, ...
        'XTickLabelRotation',25);

    ylabel('share of standardized regression magnitude [%]');
    title(sprintf('Sensitivity of %s',opts.metric),'Interpreter','none');
    legend({S.tag},'Location','best');
end

harness.sensitivity_final(MC,metric=opts.metric);

%% 5. Paired RMS improvement
out.paired = struct();

if size(MC.rms_p,2) >= 2
    base = MC.rms_p(:,1);
    adapt = MC.rms_p(:,2);

    ok = isfinite(base) & isfinite(adapt) & base > eps;
    imp = 100*(1-adapt(ok)./base(ok));

    out.paired.n = nnz(ok);
    out.paired.improvement = imp;
    out.paired.improved_fraction = mean(adapt(ok) < base(ok));

    figure( ...
        'Name',sprintf('%s - paired improvement',opts.prefix), ...
        'Color','w', ...
        'Position',[150 150 900 430]);

    histogram(imp,30,'EdgeColor','none');
    grid on
    xline(0,'--','no improvement');
    xlabel('paired RMS improvement of PBMRAC over baseline [%]');
    ylabel('count');
    title('Paired performance across identical sampled plants');

    fprintf('\n=== Paired RMS improvement ===\n');
    fprintf('valid pairs             : %d\n',nnz(ok));
    fprintf('PBMRAC better than base : %.1f%% of samples\n', ...
        100*out.paired.improved_fraction);
    fprintf('mean improvement        : %.1f%%\n',mean(imp));
    fprintf('median improvement      : %.1f%%\n',median(imp));
    fprintf('5th percentile          : %.1f%%\n', ...
        harness.percentile_final(imp,5));
end

%% 6. Robustness figure
hasRobust = isfield(MC,'max_tilt_deg') && ...
            isfield(MC,'sat_upper_frac') && ...
            isfield(MC,'max_omega_eq');

if hasRobust
    figure( ...
        'Name',sprintf('%s - actuator and attitude robustness',opts.prefix), ...
        'Color','w', ...
        'Position',[170 80 1100 800]);

    tags = MC.tags;

    subplot(2,2,1)
    hold on
    for c = 1:numel(tags)
        x = MC.max_omega_eq(:,c);
        histogram(x(isfinite(x)),30,'DisplayStyle','stairs','LineWidth',1.5);
    end
    if isfield(MC,'meta') && isfield(MC.meta,'uav')
        xline(MC.meta.uav.Omega_max,'--','\Omega_{max}');
    end
    grid on
    xlabel('maximum equivalent rotor speed [rad/s]');
    ylabel('count');
    title('Actuator headroom');
    legend(tags,'Location','best');

    subplot(2,2,2)
    hold on
    for c = 1:numel(tags)
        x = 100*MC.sat_upper_frac(:,c);
        histogram(x(isfinite(x)),30,'DisplayStyle','stairs','LineWidth',1.5);
    end
    grid on
    xlabel('upper-limit occupancy [% of run]');
    ylabel('count');
    title('Physical actuator saturation');
    legend(tags,'Location','best');

    subplot(2,2,3)
    hold on
    for c = 1:numel(tags)
        x = MC.max_tilt_deg(:,c);
        histogram(x(isfinite(x)),30,'DisplayStyle','stairs','LineWidth',1.5);
    end
    grid on
    xlabel('max(|roll|,|pitch|) [deg]');
    ylabel('count');
    title('Maximum attitude excursion');
    legend(tags,'Location','best');

    subplot(2,2,4)
    if isfield(MC,'min_lambda_hat') && size(MC.min_lambda_hat,2) >= 2
        x = MC.min_lambda_hat(:,2);
        x = x(isfinite(x));
        histogram(x,30,'EdgeColor','none');
        hold on

        if isfield(MC,'meta') && isfield(MC.meta,'controller')
            C = MC.meta.controller;
            innerLow = C.lambda_proj_center-C.lambda_proj_r0;
            outerLow = C.lambda_proj_center-C.lambda_proj_r1;
            xline(innerLow,'--','projection starts');
            xline(outerLow,':','outer wall');
        end

        grid on
        xlabel('minimum \Lambda_hat in each PBMRAC run');
        ylabel('count');
        title('Adaptive-effectiveness margin');
    else
        axis off
        text(0.5,0.5,'Lambda_hat metric unavailable', ...
            'HorizontalAlignment','center');
    end
end

%% Console summary
harness.mc_summary_final(MC);

%% Worst tails
out.worst_baseline = [];
out.worst_pbmrac   = [];

if size(MC.rms_p,2) >= 1
    [~,ord] = sort(MC.rms_p(:,1),'descend','MissingPlacement','last');
    ord = ord(isfinite(MC.rms_p(ord,1)));
    out.worst_baseline = reshape(ord(1:min(opts.worst,numel(ord))),1,[]);

    printWorst( ...
        MC, ...
        out.worst_baseline, ...
        'Worst baseline RMS runs');
end

if size(MC.rms_p,2) >= 2
    [~,ord] = sort(MC.rms_p(:,2),'descend','MissingPlacement','last');
    ord = ord(isfinite(MC.rms_p(ord,2)));
    out.worst_pbmrac = reshape(ord(1:min(opts.worst,numel(ord))),1,[]);

    printWorst( ...
        MC, ...
        out.worst_pbmrac, ...
        'Worst PBMRAC RMS runs');
end

if isfield(MC,'failed') && any(MC.failed(:))
    fprintf('\n%d simulation(s) failed and were excluded.\n',nnz(MC.failed));
end

if ~isempty(opts.export)
    viz.export(opts.export);
end

if isfield(MC,'meta') && isfield(MC.meta,'file')
    fprintf('\nResults file: %s\n',MC.meta.file);
end
end


function printWorst(MC,ord,titleText)
fprintf('\n=== %s ===\n',titleText);

hasAdapt = size(MC.rms_p,2) >= 2;
hasRobust = isfield(MC,'sat_upper_frac') && isfield(MC,'max_tilt_deg');

if hasAdapt && hasRobust
    fprintf(['%5s %7s %-22s %-28s %9s %9s ' ...
             '%8s %8s %8s\n'], ...
        'run','mass','r_c [m]','lambda (6 rotors)', ...
        'RMS base','RMS adapt','sat A','tilt A','minLam');
    fprintf('%s\n',repmat('-',1,115));
else
    fprintf('%5s %7s %-22s %-28s %9s\n', ...
        'run','mass','r_c [m]','lambda (6 rotors)','RMS');
    fprintf('%s\n',repmat('-',1,80));
end

for i = ord
    rcStr  = sprintf('%+.3f %+.3f %+.3f',MC.samples.r_c(i,:));
    lamStr = strtrim(sprintf('%.2f ',MC.samples.lambda(i,:)));

    if hasAdapt && hasRobust
        minLam = NaN;
        if isfield(MC,'min_lambda_hat') && size(MC.min_lambda_hat,2) >= 2
            minLam = MC.min_lambda_hat(i,2);
        end

        fprintf(['%5d %7.3f [%-20s] [%-26s] %9.4f %9.4f ' ...
                 '%7.2f%% %7.2f %8.4f\n'], ...
            i, ...
            MC.samples.mass(i), ...
            rcStr, ...
            lamStr, ...
            MC.rms_p(i,1), ...
            MC.rms_p(i,2), ...
            100*MC.sat_upper_frac(i,2), ...
            MC.max_tilt_deg(i,2), ...
            minLam);
    else
        fprintf('%5d %7.3f [%-20s] [%-26s] %9.4f\n', ...
            i,MC.samples.mass(i),rcStr,lamStr,MC.rms_p(i,1));
    end
end
end
