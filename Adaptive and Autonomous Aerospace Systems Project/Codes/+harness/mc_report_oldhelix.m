function out = mc_report_oldhelix(MC, opts)
%MC_REPORT_OLDHELIX  Task 3.4 plots/tables for the isolated old-helix campaign.
%
%   harness.mc_report_oldhelix(MC)
%   harness.mc_report_oldhelix(MC, export='Task3_montecarlo_oldhelix')
%
%   This helper is deliberately given a unique name so it can coexist with the
%   project's original harness.mc_report.  It understands the extended sample
%   set [package mass, package CoM, six rotor-effectiveness factors].

arguments
    MC            struct
    opts.metric   (1,:) char = 'rms_p'
    opts.tol      (1,1) double = 0.02
    opts.worst    (1,1) double {mustBeInteger,mustBePositive} = 5
    opts.export   (1,:) char = ''
    opts.prefix   (1,:) char = 'Task 3.4 old helix'
end

%% 1. Convergence
[~, out.convergence] = viz.convergence(MC, metric=opts.metric, tol=opts.tol, ...
    name=sprintf('%s - convergence diagnostic', opts.prefix));

%% 2. Metric distributions
viz.histograms(MC, metrics={'rms_p','peak_p','ss_p'}, ...
    name=sprintf('%s - metric distributions', opts.prefix));

%% 3. Sampled uncertainties
hasRc = isfield(MC.samples,'r_c');
if hasRc
    rc = MC.samples.r_c;
    figure('Name',sprintf('%s - sampled uncertainties',opts.prefix), ...
        'Color','w','Position',[120 120 1000 760]);
    subplot(2,2,1)
    histogram(MC.samples.mass,30,'EdgeColor','none'); grid on
    xlabel('package mass [kg]'); ylabel('count'); title('Sampled payload mass');
    subplot(2,2,2)
    histogram(MC.samples.lambda(:),30,'EdgeColor','none'); grid on
    xlabel('rotor effectiveness \lambda_r [-]'); ylabel('count');
    title('Sampled propeller effectiveness');
    subplot(2,2,3)
    scatter(rc(:,1),rc(:,2),14,'filled'); grid on; axis equal
    xlabel('r_{c,x} [m]'); ylabel('r_{c,y} [m]');
    title('Sampled package CoM: horizontal plane');
    subplot(2,2,4)
    histogram(rc(:,3),30,'EdgeColor','none'); grid on
    xlabel('r_{c,z} [m]'); ylabel('count');
    title('Sampled package CoM: vertical position');
else
    figure('Name',sprintf('%s - sampled uncertainties',opts.prefix), ...
        'Color','w','Position',[120 120 900 380]);
    subplot(1,2,1)
    histogram(MC.samples.mass,30,'EdgeColor','none'); grid on
    xlabel('package mass [kg]'); ylabel('count'); title('Sampled payload mass');
    subplot(1,2,2)
    histogram(MC.samples.lambda(:),30,'EdgeColor','none'); grid on
    xlabel('rotor effectiveness \lambda_r [-]'); ylabel('count');
    title('Sampled propeller effectiveness');
end

%% 4. Sensitivity / spread attribution
out.sensitivity = harness.sensitivity_oldhelix(MC, metric=opts.metric, quiet=true);
S = out.sensitivity;
if ~isempty(S)
    factorNames = S(1).names;
    nf = numel(factorNames);
    V = nan(numel(S),nf);
    for c = 1:numel(S), V(c,:) = 100*S(c).var_frac; end
    figure('Name',sprintf('%s - variance attribution',opts.prefix), ...
        'Color','w','Position',[130 130 900 430]);
    bar(V.','grouped'); grid on
    set(gca,'XTick',1:nf,'XTickLabel',factorNames,'XTickLabelRotation',25);
    ylabel('share of standardized regression magnitude [%]');
    title(sprintf('Sensitivity of %s',opts.metric),'Interpreter','none');
    legend({S.tag},'Location','best');
end
harness.sensitivity_oldhelix(MC, metric=opts.metric);

%% 5. Worst baseline cases
hasAdapt = size(MC.rms_p,2) >= 2;
[~,ord] = sort(MC.rms_p(:,1),'descend','MissingPlacement','last');
ord = reshape(ord(1:min(opts.worst,numel(ord))),1,[]);
fprintf('\n=== Worst %d baseline runs ===\n',numel(ord));
if hasRc
    if hasAdapt
        fprintf('%5s %7s  %-24s  %-29s %10s %10s %8s\n', ...
            'run','mass','r_c [m]','lambda (6 rotors)','RMS base','RMS adapt','change');
    else
        fprintf('%5s %7s  %-24s  %-29s %10s\n', ...
            'run','mass','r_c [m]','lambda (6 rotors)','RMS base');
    end
else
    if hasAdapt
        fprintf('%5s %8s  %-31s %11s %11s %9s\n', ...
            'run','mass','lambda (6 rotors)','RMS base','RMS adapt','change');
    else
        fprintf('%5s %8s  %-31s %11s\n', ...
            'run','mass','lambda (6 rotors)','RMS base');
    end
end

for i = ord
    lam_str = strtrim(sprintf('%.2f ',MC.samples.lambda(i,:)));
    if hasRc
        rc_str = sprintf('%+.3f %+.3f %+.3f',MC.samples.r_c(i,:));
        if hasAdapt
            fprintf('%5d %7.3f  [%-22s]  [%-27s] %10.4f %10.4f %7.1f%%\n', ...
                i,MC.samples.mass(i),rc_str,lam_str,MC.rms_p(i,1),MC.rms_p(i,2), ...
                100*(1-MC.rms_p(i,2)/MC.rms_p(i,1)));
        else
            fprintf('%5d %7.3f  [%-22s]  [%-27s] %10.4f\n', ...
                i,MC.samples.mass(i),rc_str,lam_str,MC.rms_p(i,1));
        end
    else
        if hasAdapt
            fprintf('%5d %8.3f  [%-29s] %11.4f %11.4f %8.1f%%\n', ...
                i,MC.samples.mass(i),lam_str,MC.rms_p(i,1),MC.rms_p(i,2), ...
                100*(1-MC.rms_p(i,2)/MC.rms_p(i,1)));
        else
            fprintf('%5d %8.3f  [%-29s] %11.4f\n', ...
                i,MC.samples.mass(i),lam_str,MC.rms_p(i,1));
        end
    end
end
out.worst = ord;

if isfield(MC,'failed') && any(MC.failed(:))
    fprintf('\n%d run(s) failed to produce usable output and were excluded.\n',nnz(MC.failed));
end

if ~isempty(opts.export), viz.export(opts.export); end
if isfield(MC,'meta') && isfield(MC.meta,'file')
    fprintf('\nResults file: %s\n',MC.meta.file);
end
end
