function fig = effectiveness(runs, tags, opts)
%EFFECTIVENESS  The adaptive control-effectiveness estimate Lambda_hat.
%
%   viz.effectiveness({R_ad}, {'PBMRAC'})
%   viz.effectiveness(runs, tags, lambda=lam, t_fault=60)
%
%   TWO DIFFERENT SPACES - the point of the lower panel.
%   lambda_sig is a per-ROTOR effectiveness: six numbers, one per arm, and it
%   is what the plant is actually given. Lambda_hat is a per-AXIS estimate:
%   three numbers, one per translational force channel, and it is what the
%   augmentation infers. They are not the same quantity and cannot be compared
%   entry by entry. A rotor fault reaches Lambda_hat only through the
%   allocation, which spreads the loss of one arm across whichever axes that
%   arm was contributing to.
%
%   WHAT DRIVES Lambda_hat AWAY FROM 1 (both at once, in this scenario)
%     * added payload mass - the augmentation's model divides by the NOMINAL
%       mass, so a heavier vehicle produces less acceleration per unit
%       commanded force and the estimator reads that as reduced effectiveness;
%     * genuine rotor degradation - less thrust delivered per unit command.
%   The estimate cannot separate them, and it does not need to: the augmented
%   input inverts whatever total it finds.
%
%   The live law has no sigma-modification leak (robust_on uses projection):
%   Lambda_hat is not pulled back towards 1, and without persistent excitation
%   its final value is not an identified physical parameter.

arguments
    runs         cell
    tags         cell
    opts.lambda  struct = struct([])   % the commanded per-rotor lambda_sig
    opts.t_fault (1,1) double = NaN
    opts.name    (1,:) char = 'Estimated control effectiveness'
end

hasLam   = ~isempty(fieldnames(opts.lambda));
nPanels  = 1 + double(hasLam);

fig = figure('Name', opts.name, 'Color','w', 'Position',[110 110 900 300*nPanels+120]);

%% ---- Panel 1: the estimate ------------------------------------------
subplot(nPanels,1,1); hold on; grid on
axlbl  = {'\Lambda_{f,x}','\Lambda_{f,y}','\Lambda_{f,z}'};
styles = {'-','--',':','-.'};
plotted = false;

for k = 1:numel(runs)
    R = runs{k};
    if ~isfield(R,'lambda_hat') || isempty(R.lambda_hat), continue, end
    set(gca,'ColorOrderIndex',1);
    plot(R.lambda_hat_t, R.lambda_hat, styles{1+mod(k-1,numel(styles))}, ...
         'LineWidth', 2);
    plotted = true;
end

if ~plotted
    close(fig); fig = gobjects(1);
    warning('viz:effectiveness:notLogged', ...
            'No run carries lambda_hat - enable that To-Workspace block.');
    return
end

yline(1, ':', 'nominal', 'Color',[0.4 0.4 0.4], 'HandleVisibility','off');
if ~isnan(opts.t_fault)
    xline(opts.t_fault, '--', 'fault', 'Color',[0.8 0.2 0.2], ...
          'LineWidth', 2, 'HandleVisibility','off');
end
ylabel('$\hat{\Lambda}_f$  [-]', 'Interpreter','latex');
legend(axlbl, 'Location','best');
if numel(runs) > 1
    title(sprintf('%s  (solid: %s, dashed: %s)', opts.name, tags{1}, tags{2}));
else
    title(opts.name);
end
if nPanels == 1, xlabel('t  [s]'); end

%% ---- Panel 2: the per-rotor truth that caused it ---------------------
if hasLam
    subplot(nPanels,1,2); hold on; grid on
    stairs(opts.lambda.time, opts.lambda.signals.values, 'LineWidth', 1.8);
    if ~isnan(opts.t_fault)
        xline(opts.t_fault, '--', 'Color',[0.8 0.2 0.2], 'HandleVisibility','off');
    end
    ylim([0 1.15]);
    xlabel('t  [s]'); ylabel('\lambda_r  [-]');
    title('Commanded per-rotor effectiveness (the applied fault)');
    legend(arrayfun(@(r) sprintf('rotor %d', r), 1:size(opts.lambda.signals.values,2), ...
           'UniformOutput', false), 'Location','eastoutside');
end

%% ---- Numbers ---------------------------------------------------------
fprintf('\n=== Estimated control effectiveness ===\n');
fprintf('%-16s %10s %10s %10s %10s\n', 'case','axis','min','post-fault','final');
fprintf('%s\n', repmat('-',1,62));
for k = 1:numel(runs)
    R = runs{k};
    if ~isfield(R,'lambda_hat') || isempty(R.lambda_hat), continue, end
    L = R.lambda_hat;  tL = R.lambda_hat_t;
    for j = 1:size(L,2)
        post = NaN;
        if ~isnan(opts.t_fault)
            m = tL >= opts.t_fault;
            if any(m), post = min(L(m,j)); end
        end
        name = ''; if j == 1, name = tags{k}; end
        axn = 'xyz';    % a char literal cannot be indexed inline in MATLAB
        fprintf('%-16s %10s %10.4f %10.4f %10.4f\n', name, axn(j), ...
                min(L(:,j)), post, L(end,j));
    end
end
end
