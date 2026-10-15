function out = mc_report(MC, opts)
%MC_REPORT  Every Task 3.4 figure and table, from a campaign struct.
%
%   harness.mc_report(MC)
%   harness.mc_report(MC, export='Task3_montecarlo')
%   out = harness.mc_report(MC, worst=10)
%
%   Works on ANY campaign struct - one just computed, or one reloaded from a
%   .mat - because it only reads the reduced metrics and the stored samples.
%   Shared by run_task3_montecarlo and run_montecarlo so the two runners cannot
%   drift into producing different figures from the same data.
%
%   READ THE OUTPUT IN THIS ORDER
%     1. convergence  - is N large enough for any of the rest to mean anything?
%                       Deliberately first: if the campaign has not converged,
%                       the distributions and percentiles below are not
%                       trustworthy and nothing else should be quoted.
%     2. distributions- the metric histograms themselves
%     3. samples      - what was actually drawn, for the appendix
%     4. sensitivity  - which uncertainty is responsible for the spread
%     5. worst cases  - the tail, named and reproducible from the stored seed
%
%   Returns a struct with the convergence and sensitivity diagnostics, so a
%   caller can assert on them rather than reading the printout.

arguments
    MC            struct
    opts.metric   (1,:) char = 'rms_p'
    opts.tol      (1,1) double = 0.02
    opts.worst    (1,1) double {mustBeInteger,mustBePositive} = 5
    opts.export   (1,:) char = ''      % non-empty -> viz.export to that folder
    opts.prefix   (1,:) char = 'Task 3.4'
end

%% ---- 1. Convergence --------------------------------------------------
[~, out.convergence] = viz.convergence(MC, metric = opts.metric, tol = opts.tol, ...
                        name = sprintf('%s - convergence diagnostic', opts.prefix));

%% ---- 2. Metric distributions ----------------------------------------
viz.histograms(MC, metrics = {'rms_p','peak_p','ss_p'}, ...
               name = sprintf('%s - metric distributions', opts.prefix));

%% ---- 3. What was drawn ----------------------------------------------
figure('Name', sprintf('%s - sampled uncertainties', opts.prefix), ...
       'Color','w', 'Position',[120 120 900 380]);
subplot(1,2,1)
histogram(MC.samples.mass, 30, 'FaceColor',[0.2 0.5 0.8], 'EdgeColor','none'); grid on
xlabel('package mass  [kg]'); ylabel('count'); title('Sampled payload mass');
subplot(1,2,2)
histogram(MC.samples.lambda(:), 30, 'FaceColor',[0.85 0.4 0.1], 'EdgeColor','none'); grid on
xlabel('rotor effectiveness  \lambda_r  [-]'); ylabel('count');
title('Sampled propeller efficiency (all rotors)');

%% ---- 4. Attribution of the spread -----------------------------------
% The campaign varies payload mass and rotor efficiency TOGETHER in every run,
% which is the realistic deployment case. The price is that a raw scatter is
% confounded: when mass happens to be high, lambda is whatever it happened to
% be. The regression recovers each factor's contribution from those same joint
% samples, holding the others fixed, so no second campaign is needed.
[~, out.sensitivity] = viz.sensitivity(MC, metric = opts.metric, ...
                        name = sprintf('%s - variance attribution', opts.prefix));
harness.sensitivity(MC, metric = opts.metric);

%% ---- 5. The tail ----------------------------------------------------
% The unlucky corner of the parameter space is what a delivery vehicle actually
% has to survive, so list the draws that hurt the baseline most - and what the
% augmentation managed to do with those same draws.
hasAdapt = size(MC.rms_p, 2) >= 2;
[~, ord]  = sort(MC.rms_p(:,1), 'descend', 'MissingPlacement','last');
ord       = reshape(ord(1:min(opts.worst, numel(ord))), 1, []);

fprintf('\n=== Worst %d baseline runs ===\n', numel(ord));
if hasAdapt
    fprintf('%5s %8s  %-31s %11s %11s %9s\n', ...
            'run','mass','lambda (6 rotors)','RMS base','RMS adapt','change');
else
    fprintf('%5s %8s  %-31s %11s\n', 'run','mass','lambda (6 rotors)','RMS base');
end

for i = ord
    % sprintf with a repeated specifier cycles over the vector and returns a
    % plain char row. compose/join return string or cell depending on the input
    % shape, and fprintf will not accept a cell.
    lam_str = strtrim(sprintf('%.2f ', MC.samples.lambda(i,:)));
    if hasAdapt
        fprintf('%5d %8.3f  [%-29s] %11.4f %11.4f %8.1f%%\n', ...
                i, MC.samples.mass(i), lam_str, MC.rms_p(i,1), MC.rms_p(i,2), ...
                100*(1 - MC.rms_p(i,2)/MC.rms_p(i,1)));
    else
        fprintf('%5d %8.3f  [%-29s] %11.4f\n', ...
                i, MC.samples.mass(i), lam_str, MC.rms_p(i,1));
    end
end
out.worst = ord;

if isfield(MC,'failed') && any(MC.failed(:))
    fprintf('\n%d run(s) failed to produce usable output and were excluded.\n', ...
            nnz(MC.failed));
end

%% ---- Save ------------------------------------------------------------
if ~isempty(opts.export)
    viz.export(opts.export);
end
if isfield(MC,'meta') && isfield(MC.meta,'file')
    fprintf('\nResults file: %s\n', MC.meta.file);
end
end
