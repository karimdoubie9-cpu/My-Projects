%% plot_montecarlo_final_polished.m
% Post-process an EXISTING Task 3.4 Monte Carlo .mat file.
%
% No simulations are launched.
%
% Main outputs:
%   01_sampled_uncertainties
%   02_paired_performance
%   03_tracking_distributions
%   04_challenging_tail
%   05_convergence_robust
%
% The "typical" / "challenging-tail" split is ONLY for visualization.
% ALL samples are retained in the numerical statistics and CSV summaries.

%clear; close all; clc;

%% ------------------------------------------------------------------------
% USER SETTINGS
% -------------------------------------------------------------------------

MAT_FILE = 'results\montecarlo_final_glamz10_20260925_160627.mat';   % unified gains, hedging, new helix
                    % set to '' -> file chooser

OUT_DIR = fullfile('Figures','Task3_montecarlo_polished');

% Visualization-only definition of a challenging/tail realization.
% A plant is labelled "tail" if EITHER controller has:
%   - RMS position error > 1 m, OR
%   - upper actuator saturation for > 5% of the run, OR
%   - max roll/pitch excursion > 30 deg.
%
% This does NOT remove a sample from statistics.
TAIL_RMS_M      = 1.0;
TAIL_SAT_FRAC   = 0.05;
TAIL_TILT_DEG   = 30;

N_WORST = 10;       % number of challenging cases listed / highlighted
N_BINS  = 24;

% Report styling
FONT_NAME = 'Arial';
FONT_SIZE = 11;
LINE_W    = 1.8;
MARKER_SZ = 34;

% Consistent controller colors
C_BASE = [0.35 0.35 0.35];
C_ADAP = [0.80 0.00 0.70];
C_TAIL = [0.85 0.20 0.15];
C_TYP  = [0.20 0.45 0.75];

%% ------------------------------------------------------------------------
% LOAD EXISTING MONTE CARLO
% -------------------------------------------------------------------------

if isempty(MAT_FILE)
    [f,p] = uigetfile('*.mat','Select Task 3.4 Monte Carlo results file');
    if isequal(f,0)
        error('No .mat file selected.');
    end
    MAT_FILE = fullfile(p,f);
end

S = load(MAT_FILE);

if isfield(S,'MC')
    MC = S.MC;
else
    names = fieldnames(S);
    idx = find(structfun(@isstruct,S),1,'first');
    if isempty(idx)
        error('No MC structure found in %s.', MAT_FILE);
    end
    MC = S.(names{idx});
    warning('Variable MC not found; using structure "%s".', names{idx});
end

required = {'rms_p','peak_p','ss_p','samples','tags'};
for k = 1:numel(required)
    assert(isfield(MC,required{k}), 'Missing MC.%s', required{k});
end
assert(isfield(MC.samples,'mass'),   'Missing MC.samples.mass');
assert(isfield(MC.samples,'lambda'), 'Missing MC.samples.lambda');

tags = string(MC.tags(:));
nc = numel(tags);
assert(nc >= 2,'Expected at least Baseline and PBMRAC configurations.');

iB = find(contains(lower(tags),'base'),1);
iA = find(contains(lower(tags),'pbmrac') | contains(lower(tags),'adapt'),1);
if isempty(iB), iB = 1; end
if isempty(iA), iA = min(2,nc); end

tagB = char(tags(iB));
tagA = char(tags(iA));

rB = MC.rms_p(:,iB);
rA = MC.rms_p(:,iA);
pB = MC.peak_p(:,iB);
pA = MC.peak_p(:,iA);
sB = MC.ss_p(:,iB);
sA = MC.ss_p(:,iA);

n = numel(rB);
mass = MC.samples.mass(:);
lam  = MC.samples.lambda;

assert(numel(mass)==n, 'Mass sample count does not match MC.rms_p.');
assert(size(lam,1)==n,  'Lambda sample count does not match MC.rms_p.');

valid = isfinite(rB) & isfinite(rA);

satB  = metric_col(MC,'sat_upper_frac',iB,n);
satA  = metric_col(MC,'sat_upper_frac',iA,n);
tiltB = metric_col(MC,'max_tilt_deg',iB,n);
tiltA = metric_col(MC,'max_tilt_deg',iA,n);
minHatA = metric_col(MC,'min_lambda_hat',iA,n);

% If a robustness metric was not stored, do not let NaNs classify the plant.
satForTail  = [zero_if_nan(satB),  zero_if_nan(satA)];
tiltForTail = [zero_if_nan(tiltB), zero_if_nan(tiltA)];

tail = valid & ( ...
       max([rB rA],[],2) > TAIL_RMS_M ...
    | max(satForTail,[],2) > TAIL_SAT_FRAC ...
    | max(tiltForTail,[],2) > TAIL_TILT_DEG );

typical = valid & ~tail;

fprintf('\n============================================================\n');
fprintf(' POLISHED TASK 3.4 POST-PROCESSING\n');
fprintf('============================================================\n');
fprintf('File                       : %s\n', MAT_FILE);
fprintf('Valid paired plants        : %d / %d\n', nnz(valid),n);
fprintf('Typical-region plants      : %d\n',nnz(typical));
fprintf('Challenging/tail plants    : %d\n',nnz(tail));
fprintf('Tail visualization rule    : RMS > %.2f m OR sat > %.1f%% OR tilt > %.1f deg\n',...
    TAIL_RMS_M,100*TAIL_SAT_FRAC,TAIL_TILT_DEG);
fprintf('IMPORTANT: all %d valid plants remain in the reported statistics.\n\n',nnz(valid));

if ~exist(OUT_DIR,'dir')
    mkdir(OUT_DIR);
end

%% ------------------------------------------------------------------------
% ALL-SAMPLE NUMERICAL SUMMARY
% -------------------------------------------------------------------------

stats = table( ...
    string({tagB;tagA}), ...
    [mean(rB(valid)); mean(rA(valid))], ...
    [median(rB(valid)); median(rA(valid))], ...
    [pct(rB(valid),95); pct(rA(valid),95)], ...
    [pct(rB(valid),99); pct(rA(valid),99)], ...
    [max(rB(valid)); max(rA(valid))], ...
    'VariableNames',{'Controller','Mean_RMS_m','Median_RMS_m','P95_RMS_m','P99_RMS_m','Max_RMS_m'});

disp(stats);
writetable(stats,fullfile(OUT_DIR,'all_sample_tracking_summary.csv'));

imp = nan(n,1);
imp(valid) = 100*(1 - rA(valid)./rB(valid));

fprintf('PBMRAC improves paired RMS on : %.1f %% of valid plants\n', ...
    100*mean(imp(valid)>0));
fprintf('Median paired improvement    : %.1f %%\n', median(imp(valid)));
fprintf('5th percentile improvement   : %.1f %%\n', pct(imp(valid),5));

%% ------------------------------------------------------------------------
% FIGURE 1 — SAMPLED UNCERTAINTIES
% -------------------------------------------------------------------------

fig = newfig('Task 3.4 - sampled uncertainties');
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

nexttile;
histogram(mass,N_BINS,'FaceColor',C_TYP,'EdgeColor','none');
xlabel('package mass [kg]'); ylabel('count');
title('Payload mass');
grid on; box on;

nexttile;
histogram(lam(:),N_BINS,'FaceColor',C_TYP,'EdgeColor','none');
xlabel('rotor effectiveness \lambda_r [-]'); ylabel('count');
title(sprintf('All six rotor-effectiveness draws (%d total)',numel(lam)));
grid on; box on;

lamMin = min(lam,[],2);
nexttile;
histogram(lamMin,N_BINS,'FaceColor',C_TYP,'EdgeColor','none');
xlabel('min_i(\lambda_{r,i}) per plant [-]'); ylabel('count');
title(sprintf('Weakest rotor per plant (%d plants)',n));
grid on; box on;

lamSpread = max(lam,[],2)-min(lam,[],2);
nexttile;
histogram(lamSpread,N_BINS,'FaceColor',C_TYP,'EdgeColor','none');
xlabel('max(\lambda_r)-min(\lambda_r) [-]'); ylabel('count');
title(sprintf('Rotor-effectiveness asymmetry (%d plants)',n));
grid on; box on;

title(tl,'Task 3.4 sampled uncertainties');
style_all(fig,FONT_NAME,FONT_SIZE);
save_both(fig,OUT_DIR,'01_sampled_uncertainties');

%% ------------------------------------------------------------------------
% FIGURE 2 — PAIRED PERFORMANCE
% -------------------------------------------------------------------------

fig = newfig('Task 3.4 - paired performance');
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

% ---- A: paired log-log scatter, ALL samples retained -------------------
ax = nexttile;
hold(ax,'on');

scatter(ax,rB(typical),rA(typical),MARKER_SZ,C_TYP,'filled', ...
    'MarkerFaceAlpha',0.55,'DisplayName','Typical region');
scatter(ax,rB(tail),rA(tail),MARKER_SZ*1.25,C_TAIL,'filled', ...
    'MarkerEdgeColor','k','LineWidth',0.5,'DisplayName','Challenging tail');

positiveVals = [rB(valid);rA(valid)];
positiveVals = positiveVals(positiveVals>0 & isfinite(positiveVals));
lo = 10^floor(log10(min(positiveVals)));
hi = 10^ceil(log10(max(positiveVals)));

plot(ax,[lo hi],[lo hi],'k--','LineWidth',1.3,'DisplayName','equal performance');

set(ax,'XScale','log','YScale','log');
xlim(ax,[lo hi]); ylim(ax,[lo hi]);
axis(ax,'square');
xlabel(ax,'Baseline RMS |e_p| [m]');
ylabel(ax,'PBMRAC RMS |e_p| [m]');
title(ax,'Paired performance — all sampled plants');
grid(ax,'on'); box(ax,'on');
legend(ax,'Location','northwest');

text(ax,0.04,0.04, ...
    sprintf('PBMRAC better: %.1f%%',100*mean(imp(valid)>0)), ...
    'Units','normalized','FontWeight','bold','BackgroundColor','w','Margin',4);

% ---- B: paired improvement in typical region only ----------------------
ax = nexttile;
impTyp = imp(typical);
if isempty(impTyp)
    text(ax,0.5,0.5,'No samples satisfy the typical-region mask.', ...
        'HorizontalAlignment','center');
else
    loImp = pct(impTyp,1);
    hiImp = pct(impTyp,99);
    pad = max(5,0.08*(hiImp-loImp));
    edges = linspace(loImp-pad,hiImp+pad,N_BINS+1);

    histogram(ax,impTyp,edges,'FaceColor',C_ADAP,'EdgeColor','none');
    hold(ax,'on');
    xline(ax,0,'k--','LineWidth',1.4,'Label','no improvement');
    xline(ax,median(impTyp),'Color',C_ADAP,'LineStyle','-', ...
        'LineWidth',1.8,'Label',sprintf('median %.1f%%',median(impTyp)));
    xlabel(ax,'paired RMS improvement of PBMRAC [%]');
    ylabel(ax,'count');
    title(ax,{'Typical operating region','challenging tail shown separately'});
    grid(ax,'on'); box(ax,'on');
end

title(tl,'Paired controller performance');
style_all(fig,FONT_NAME,FONT_SIZE);
save_both(fig,OUT_DIR,'02_paired_performance');

%% ------------------------------------------------------------------------
% FIGURE 3 — TRACKING DISTRIBUTIONS
% -------------------------------------------------------------------------

fig = newfig('Task 3.4 - tracking distributions');
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% Zoomed distributions: typical region only for readability.
% ALL samples remain in the statistics and in the full-range CDF panel.

overlay_hist(nexttile, rB(typical),rA(typical),N_BINS,C_BASE,C_ADAP, ...
    'RMS |e_p| [m]','Typical-region RMS');

overlay_hist(nexttile, pB(typical),pA(typical),N_BINS,C_BASE,C_ADAP, ...
    'peak |e_p| [m]','Typical-region peak error');

overlay_hist(nexttile, sB(typical),sA(typical),N_BINS,C_BASE,C_ADAP, ...
    'steady-state |e_p| [m]','Typical-region steady-state error');

% Full empirical CDF of RMS on log-x: ALL valid samples.
ax = nexttile;
hold(ax,'on');
[xb,yb] = empirical_cdf(rB(valid));
[xa,ya] = empirical_cdf(rA(valid));
plot(ax,xb,yb,'Color',C_BASE,'LineWidth',LINE_W,'DisplayName',tagB);
plot(ax,xa,ya,'Color',C_ADAP,'LineWidth',LINE_W,'DisplayName',tagA);
set(ax,'XScale','log');

positiveCDF = [xb(:); xa(:)];
positiveCDF = positiveCDF(isfinite(positiveCDF) & positiveCDF > 0);
if ~isempty(positiveCDF)
    xlim(ax,[min(positiveCDF) max(positiveCDF)]);
end

% Make the CDF interpretation explicit:
% y = fraction of Monte Carlo runs whose RMS error is <= the x-axis value.
yline(ax,0.50,':','Color',[0.55 0.55 0.55], ...
    'LineWidth',1.0,'HandleVisibility','off');
yline(ax,0.95,':','Color',[0.55 0.55 0.55], ...
    'LineWidth',1.0,'HandleVisibility','off');

medB = median(rB(valid));
medA = median(rA(valid));
p95B = pct(rB(valid),95);
p95A = pct(rA(valid),95);

plot(ax,medB,0.50,'o','Color',C_BASE,'MarkerFaceColor',C_BASE, ...
    'MarkerSize',6,'HandleVisibility','off');
plot(ax,medA,0.50,'o','Color',C_ADAP,'MarkerFaceColor',C_ADAP, ...
    'MarkerSize',6,'HandleVisibility','off');
plot(ax,p95B,0.95,'s','Color',C_BASE,'MarkerFaceColor',C_BASE, ...
    'MarkerSize',6,'HandleVisibility','off');
plot(ax,p95A,0.95,'s','Color',C_ADAP,'MarkerFaceColor',C_ADAP, ...
    'MarkerSize',6,'HandleVisibility','off');

xlabel(ax,'RMS |e_p| [m]');
ylabel(ax,'fraction of runs with RMS \leq x');
title(ax,'Empirical CDF — all 500 runs');
ylim(ax,[0 1]);
grid(ax,'on'); box(ax,'on');
legend(ax,'Location','southeast');

text(ax,0.03,0.08,'curve farther left = lower tracking error', ...
    'Units','normalized','FontSize',9,'BackgroundColor','w','Margin',3);

title(tl,{ ...
    'Tracking-error distributions', ...
    'Zoomed panels show typical region; full CDF retains every sample'});
style_all(fig,FONT_NAME,FONT_SIZE);
save_both(fig,OUT_DIR,'03_tracking_distributions');

%% ------------------------------------------------------------------------
% FIGURE 4 — CHALLENGING TAIL
% -------------------------------------------------------------------------

% Rank by the larger RMS of the two controllers.
severity = max([rB rA],[],2);
severity(~valid) = -Inf;

tailCandidates = find(tail);
[~,ordLocal] = sort(severity(tailCandidates),'descend');
tailCandidates = tailCandidates(ordLocal);

nShow = min(N_WORST,numel(tailCandidates));

if nShow > 0
    worst = tailCandidates(1:nShow);
else
    % Fallback: if no case crosses the explicit tail rule, show worst overall.
    [~,ord] = sort(severity,'descend');
    worst = ord(1:min(N_WORST,nnz(valid)));
    warning('No explicit tail cases found; Figure 4 shows worst valid runs instead.');
end

runLabels = compose('%d',worst);

fig = newfig('Task 3.4 - challenging-tail cases');
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% ---- A: RMS in worst cases ---------------------------------------------
ax = nexttile;
Y = [rB(worst),rA(worst)];
bar(ax,Y,'grouped');
set(ax,'YScale','log');
xticks(ax,1:numel(worst)); xticklabels(ax,runLabels);
xlabel(ax,'Monte Carlo run index'); ylabel(ax,'RMS |e_p| [m] — log scale');
title(ax,'Worst challenging cases');
legend(ax,{tagB,tagA},'Location','northwest');
grid(ax,'on'); box(ax,'on');

% ---- B: actuator saturation in same cases ------------------------------
ax = nexttile;
Y = 100*[satB(worst),satA(worst)];
bar(ax,Y,'grouped');
xticks(ax,1:numel(worst)); xticklabels(ax,runLabels);
xlabel(ax,'Monte Carlo run index'); ylabel(ax,'upper-limit occupancy [% of run]');
title(ax,'Actuator saturation');
legend(ax,{tagB,tagA},'Location','northwest');
ymax = max(Y(isfinite(Y)));
if isempty(ymax), ymax = 5; end
ylim(ax,[0 max(5,1.05*ymax)]);
grid(ax,'on'); box(ax,'on');

% ---- C: attitude in same cases -----------------------------------------
ax = nexttile;
Y = [tiltB(worst),tiltA(worst)];
bar(ax,Y,'grouped');
xticks(ax,1:numel(worst)); xticklabels(ax,runLabels);
xlabel(ax,'Monte Carlo run index'); ylabel(ax,'max(|roll|,|pitch|) [deg]');
title(ax,'Attitude excursion');
legend(ax,{tagB,tagA},'Location','northwest');
grid(ax,'on'); box(ax,'on');

% ---- D: where tail appears in uncertainty space ------------------------
ax = nexttile;
hold(ax,'on');
scatter(ax,mass(typical),lamMin(typical),24,[0.70 0.70 0.70],'filled', ...
    'MarkerFaceAlpha',0.35,'DisplayName','Typical');
scatter(ax,mass(tail),lamMin(tail),42,C_TAIL,'filled', ...
    'MarkerEdgeColor','k','LineWidth',0.4,'DisplayName','Tail');

for j = 1:numel(worst)
    text(ax,mass(worst(j)),lamMin(worst(j)),sprintf('  %d',worst(j)), ...
        'FontSize',8,'Color',[0.15 0.15 0.15]);
end

xlabel(ax,'package mass [kg]');
ylabel(ax,'minimum rotor effectiveness in plant [-]');
title(ax,'Location of challenging cases');
grid(ax,'on'); box(ax,'on');
legend(ax,'Location','southwest');

title(tl,{ ...
    'Challenging-tail cases', ...
    'These cases are retained in all campaign statistics'});
style_all(fig,FONT_NAME,FONT_SIZE);
save_both(fig,OUT_DIR,'04_challenging_tail');

% Save explicit tail table.
meanLam = mean(lam,2);
spreadLam = max(lam,[],2)-min(lam,[],2);

% Classify the dominant mechanism for each displayed tail case.
TailClass = strings(numel(worst),1);
for j = 1:numel(worst)
    i = worst(j);
    sat_i = max([satB(i), satA(i)],[],'omitnan');
    if isfinite(sat_i) && sat_i > TAIL_SAT_FRAC
        TailClass(j) = "actuator saturation";
    elseif isfinite(minHatA(i)) && minHatA(i) <= 0.45
        TailClass(j) = "adaptive projection region";
    else
        TailClass(j) = "tracking/attitude outlier";
    end
end

TailTable = table( ...
    worst(:), ...
    mass(worst), ...
    meanLam(worst), ...
    lamMin(worst), ...
    spreadLam(worst), ...
    rB(worst), ...
    rA(worst), ...
    100*satB(worst), ...
    100*satA(worst), ...
    tiltB(worst), ...
    tiltA(worst), ...
    minHatA(worst), ...
    TailClass, ...
    'VariableNames',{ ...
    'Run','Mass_kg','MeanRotorEffectiveness','MinRotorEffectiveness', ...
    'RotorEffectivenessSpread','RMS_Baseline_m','RMS_PBMRAC_m', ...
    'Sat_Baseline_pct','Sat_PBMRAC_pct', ...
    'Tilt_Baseline_deg','Tilt_PBMRAC_deg','MinLambdaHat_PBMRAC', ...
    'TailClass'});

fprintf('\n=== Challenging-tail table ===\n');
disp(TailTable);
writetable(TailTable,fullfile(OUT_DIR,'challenging_tail_cases.csv'));

%% ------------------------------------------------------------------------
% FIGURE 5 — ROBUST CONVERGENCE DIAGNOSTICS
% -------------------------------------------------------------------------
% Replaces the mean +/- CI / re-ordering-spaghetti figure.
%
% These three quantities remain interpretable for a heavy-tailed campaign:
%   1) running median RMS
%   2) running p95 RMS
%   3) running fraction of paired plants improved by PBMRAC

k0 = min(10,n);
kk = (k0:n).';

runMedB = nan(size(kk));
runMedA = nan(size(kk));
runP95B = nan(size(kk));
runP95A = nan(size(kk));
runBetter = nan(size(kk));

for q = 1:numel(kk)
    k = kk(q);
    vb = rB(1:k); va = rA(1:k);
    ok = isfinite(vb) & isfinite(va);
    vb = vb(ok); va = va(ok);

    runMedB(q) = median(vb);
    runMedA(q) = median(va);
    runP95B(q) = pct(vb,95);
    runP95A(q) = pct(va,95);
    runBetter(q) = 100*mean(va < vb);
end

fig = newfig('Task 3.4 - robust convergence');
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

ax = nexttile;
hold(ax,'on');
plot(ax,kk,runMedB,'Color',C_BASE,'LineWidth',LINE_W,'DisplayName',tagB);
plot(ax,kk,runMedA,'Color',C_ADAP,'LineWidth',LINE_W,'DisplayName',tagA);
% Retrospective reference: final statistic computed using all 500 samples.
% It is NOT a theoretical target or an independent truth value.
yline(ax,median(rB(valid)),'--','Color',C_BASE,'HandleVisibility','off');
yline(ax,median(rA(valid)),'--','Color',C_ADAP,'HandleVisibility','off');
xlabel(ax,'samples used'); ylabel(ax,'running median RMS [m]');
title(ax,'Running median RMS');
grid(ax,'on'); box(ax,'on');
legend(ax,'Location','best');

ax = nexttile;
hold(ax,'on');
plot(ax,kk,runP95B,'Color',C_BASE,'LineWidth',LINE_W,'DisplayName',tagB);
plot(ax,kk,runP95A,'Color',C_ADAP,'LineWidth',LINE_W,'DisplayName',tagA);
% Retrospective reference: final p95 computed using all 500 samples.
yline(ax,pct(rB(valid),95),'--','Color',C_BASE,'HandleVisibility','off');
yline(ax,pct(rA(valid),95),'--','Color',C_ADAP,'HandleVisibility','off');
xlabel(ax,'samples used'); ylabel(ax,'running p95 RMS [m]');
title(ax,'Running p95 RMS');
grid(ax,'on'); box(ax,'on');

ax = nexttile;
plot(ax,kk,runBetter,'Color',C_ADAP,'LineWidth',LINE_W);
hold(ax,'on');
yline(ax,100*mean(rA(valid)<rB(valid)),'--k', ...
    sprintf('500-run value %.1f%%',100*mean(rA(valid)<rB(valid))));
xlabel(ax,'samples used');
ylabel(ax,'paired plants improved by PBMRAC [%]');
title(ax,'Running paired-improvement fraction');
ylim(ax,[0 100]);
grid(ax,'on'); box(ax,'on');

title(tl,{ ...
    'Monte Carlo running-statistic stability', ...
    'Dashed horizontal lines = values obtained from the complete 500-run campaign'});
style_all(fig,FONT_NAME,FONT_SIZE);
save_both(fig,OUT_DIR,'05_convergence_robust');

%% ------------------------------------------------------------------------
% OPTIONAL: PRINT SHORT REPORT-READY SUMMARY
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' REPORT-READY SUMMARY\n');
fprintf('============================================================\n');
fprintf('Paired valid plants        : %d\n',nnz(valid));
fprintf('PBMRAC lower RMS           : %.1f %%\n',100*mean(rA(valid)<rB(valid)));
fprintf('Median RMS   Baseline      : %.4f m\n',median(rB(valid)));
fprintf('Median RMS   PBMRAC        : %.4f m\n',median(rA(valid)));
fprintf('p95 RMS      Baseline      : %.4f m\n',pct(rB(valid),95));
fprintf('p95 RMS      PBMRAC        : %.4f m\n',pct(rA(valid),95));
fprintf('Median paired improvement  : %.1f %%\n',median(imp(valid)));
fprintf('5th pct paired improvement : %.1f %%\n',pct(imp(valid),5));
fprintf('Challenging-tail plants    : %d (%.1f %%)\n', ...
    nnz(tail),100*nnz(tail)/nnz(valid));
fprintf('Figures + CSVs written to  : %s\n',OUT_DIR);

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function x = metric_col(MC,name,col,n)
    x = nan(n,1);
    if ~isfield(MC,name)
        return;
    end
    A = MC.(name);
    if isvector(A) && numel(A)==n
        x = A(:);
    elseif size(A,1)==n && size(A,2)>=col
        x = A(:,col);
    end
end

function x = zero_if_nan(x)
    x(~isfinite(x)) = 0;
end

function q = pct(x,p)
%PCT  Same percentile convention as +harness/percentile_final.m.
%
% This is intentionally duplicated locally so this post-processing script
% reproduces the Monte Carlo percentile values without depending on MATLAB
% package/path resolution.

    x = x(:);
    x = sort(x(~isnan(x)));

    if isempty(x)
        q = nan(size(p));
        return
    end

    if isscalar(x)
        q = repmat(x,size(p));
        return
    end

    p = max(0,min(100,p));

    % Same plotting positions used by harness.percentile_final:
    %   pos_i = 100*(i-0.5)/N
    pos = 100*((1:numel(x))-0.5)/numel(x);
    q   = interp1(pos,x,p,'linear','extrap');
    q   = min(max(q,x(1)),x(end));
end

function [x,y] = empirical_cdf(v)
% Toolbox-free empirical CDF.
    x = sort(v(isfinite(v) & v>0));
    y = (1:numel(x)).'/numel(x);
end

function overlay_hist(ax,a,b,nBins,c1,c2,xlab,ttl)
    a = a(isfinite(a));
    b = b(isfinite(b));
    vals = [a(:);b(:)];
    if isempty(vals)
        text(ax,0.5,0.5,'No data','HorizontalAlignment','center');
        return;
    end
    lo = min(vals); hi = max(vals);
    if hi<=lo
        hi = lo + max(1e-6,0.05*abs(lo)+1e-6);
    end
    edges = linspace(lo,hi,nBins+1);
    histogram(ax,a,edges,'DisplayStyle','stairs','LineWidth',1.8, ...
        'EdgeColor',c1,'DisplayName','Baseline');
    hold(ax,'on');
    histogram(ax,b,edges,'DisplayStyle','stairs','LineWidth',1.8, ...
        'EdgeColor',c2,'DisplayName','PBMRAC');
    xline(ax,median(a),'--','Color',c1,'LineWidth',1.2,'HandleVisibility','off');
    xline(ax,median(b),'--','Color',c2,'LineWidth',1.2,'HandleVisibility','off');
    xlabel(ax,xlab); ylabel(ax,'count');
    title(ax,ttl);
    grid(ax,'on'); box(ax,'on');
    legend(ax,'Location','best');
end

function fig = newfig(name)
    fig = figure('Name',name,'Color','w', ...
        'Units','pixels','Position',[100 80 1450 780]);
    try
        fig.ToolBar = 'none';
    catch
    end
end

function style_all(fig,fontName,fontSize)
    axs = findall(fig,'Type','axes');
    for k = 1:numel(axs)
        set(axs(k),'FontName',fontName,'FontSize',fontSize, ...
            'LineWidth',0.9,'TickDir','out');
    end
end

function save_both(fig,outDir,stem)
    pngFile = fullfile(outDir,[stem '.png']);
    pdfFile = fullfile(outDir,[stem '.pdf']);
    axs = findall(fig,'Type','axes');
    for k = 1:numel(axs)
        try
            axs(k).Toolbar.Visible = 'off';
        catch
        end
    end
    drawnow;
    exportgraphics(fig,pngFile,'Resolution',300);
    exportgraphics(fig,pdfFile,'ContentType','vector');
end
%%
%% ========================================================================
% MONTE CARLO SAMPLE-SIZE JUSTIFICATION
% ========================================================================
%
% The principal Monte Carlo question is treated as a paired Bernoulli
% experiment:
%
%       X_i = 1   if RMS_PBMRAC(i) < RMS_Baseline(i)
%             0   otherwise
%
% where Baseline and PBMRAC are evaluated on the SAME sampled plant
% (same payload mass and same six rotor-effectiveness values).
%
% Hence
%
%       p_hat = (1/N) * sum_i X_i
%
% estimates the probability that PBMRAC improves RMS tracking for a plant
% drawn from the prescribed uncertainty distributions.
%
%
% -------------------------------------------------------------------------
% 1. WHY N = 500 IS A REASONABLE MONTE CARLO SIZE
% -------------------------------------------------------------------------
%
% For a binomial proportion, the approximate sampling uncertainty scales as
%
%       SE(p_hat) ~ sqrt( p*(1-p) / N )
%
% and therefore decreases only as
%
%       1/sqrt(N).
%
% Before knowing p, the most conservative possible case is
%
%       p = 0.5,
%
% because p(1-p) is maximised at 0.25.
%
% Using a 95% confidence level (z = 1.96), the usual conservative
% sample-size estimate for a desired absolute precision epsilon is
%
%       N >= z^2 * 0.25 / epsilon^2.
%
% For epsilon = 0.05 (±5 percentage points):
%
%       N >= 1.96^2 * 0.25 / 0.05^2
%         ~= 384.2
%
% so at least 385 independent plants are required in the worst-case
% binomial setting.
%
% The present campaign uses
%
%       N = 500 > 385,
%
% and therefore exceeds this conservative requirement.
%
% This gives a theoretical justification for using 500 plants without
% relying on the observed controller performance itself.
%
%
% -------------------------------------------------------------------------
% 2. ACHIEVED PRECISION: WILSON SCORE INTERVAL
% -------------------------------------------------------------------------
%
% After the Monte Carlo has been run, the achieved uncertainty in the
% improvement probability can be quantified using a Wilson score interval.
%
% The Wilson interval is preferred to the simple Wald interval
%
%       p_hat ± z*sqrt(p_hat*(1-p_hat)/N)
%
% because the Wald approximation becomes unreliable when p_hat is close to
% 0 or 1.  Here p_hat is close to 1, so the Wilson interval is more
% appropriate.
%
% For confidence level 1-alpha, with z = 1.96 for 95% confidence:
%
%       centre =
%           ( p_hat + z^2/(2N) ) / ( 1 + z^2/N )
%
%       half_width =
%           z/(1 + z^2/N) *
%           sqrt( p_hat*(1-p_hat)/N + z^2/(4*N^2) )
%
%       CI_95 = [centre - half_width, centre + half_width].
%
% In the present campaign:
%
%       N          = 500
%       improved   = 489
%       p_hat      = 489/500 = 0.978 = 97.8%
%
% giving approximately
%
%       95% Wilson CI = [0.961, 0.988]
%
% or
%
%       96.1% <= p <= 98.8%.
%
% Therefore the statement that PBMRAC improves RMS tracking for the vast
% majority of sampled plants is not based only on the raw 97.8% count:
% the corresponding binomial confidence interval is also narrow.
%
%
% -------------------------------------------------------------------------
% 3. WHY WE DO NOT USE THE OLD "2% OF THE MEAN" RULE
% -------------------------------------------------------------------------
%
% An earlier convergence criterion considered the relative confidence
% interval of the arithmetic mean RMS:
%
%       1.96 * s / (sqrt(N)*mean(RMS)).
%
% This is not used as the primary sample-size criterion here because the
% RMS distribution is strongly right-skewed / heavy-tailed.
%
% A small number of actuator-limited plants produce RMS errors of several
% metres or tens of metres, while the majority of runs remain below
% approximately a few tenths of a metre.
%
% Consequently:
%
%       - the arithmetic mean is strongly affected by rare extreme cases;
%       - the sample standard deviation becomes very large;
%       - requiring the mean to be known to ±2% would imply an
%         unrealistically large number of simulations;
%       - such a requirement is not directly related to the principal
%         controller-comparison question.
%
% The Monte Carlo is therefore assessed primarily through:
%
%       * paired probability of PBMRAC improvement;
%       * median RMS;
%       * p95 RMS;
%       * challenging-tail / actuator-saturation cases.
%
%
% -------------------------------------------------------------------------
% 4. WHY 1000 RUNS ARE NOT AUTOMATICALLY NECESSARY
% -------------------------------------------------------------------------
%
% Statistical uncertainty decreases approximately as 1/sqrt(N).
%
% Therefore doubling the campaign from
%
%       N = 500  ->  N = 1000
%
% only reduces standard statistical uncertainty by
%
%       sqrt(500/1000) = 1/sqrt(2) ~= 0.707,
%
% i.e. about a 29% reduction.
%
% Thus 1000 runs would improve precision, but at approximately twice the
% simulation cost and without changing the expected qualitative conclusion.
%
% A larger campaign would mainly be beneficial if a more detailed
% characterisation of rare challenging-tail events were required.
%
%
% -------------------------------------------------------------------------
% 5. INTERPRETATION
% -------------------------------------------------------------------------
%
% The final justification is therefore:
%
%   - N = 500 exceeds the conservative N ~= 385 requirement for estimating
%     an unknown binomial proportion to ±5 percentage points at 95%
%     confidence.
%
%   - The observed paired result is 489/500 = 97.8% PBMRAC improvements.
%
%   - The corresponding 95% Wilson score interval is approximately
%     [96.1%, 98.8%].
%
%   - Running median, p95 and paired-improvement plots are additionally
%     inspected to verify that the main Monte Carlo conclusions are stable
%     as the sample count increases.
%
% Hence N = 500 is considered sufficient for the scope of Task 3.4.
% ========================================================================


% Optional: reproduce the numerical Wilson-interval calculation

N        = nnz(valid);                    % from the loaded campaign
nImprove = nnz(rA(valid) < rB(valid));
p_hat    = nImprove / N;

z = 1.96;                  % 95% confidence

wilsonCentre = ...
    (p_hat + z^2/(2*N)) / (1 + z^2/N);

wilsonHalfWidth = ...
    z/(1 + z^2/N) * ...
    sqrt(p_hat*(1-p_hat)/N + z^2/(4*N^2));

wilsonCI = [ ...
    wilsonCentre - wilsonHalfWidth, ...
    wilsonCentre + wilsonHalfWidth ];

% Conservative pre-run sample-size requirement:
epsilon = 0.05;            % desired +/- 5 percentage-point precision
N_required_conservative = ceil(z^2 * 0.25 / epsilon^2);

fprintf('\n=== Monte Carlo sample-size justification ===\n');
fprintf('Samples                              : %d\n', N);
fprintf('PBMRAC improvements                  : %d / %d\n', ...
    nImprove, N);
fprintf('Observed improvement probability     : %.2f %%\n', ...
    100*p_hat);
fprintf('95%% Wilson interval                  : [%.2f, %.2f] %%\n', ...
    100*wilsonCI(1), 100*wilsonCI(2));
fprintf('Worst-case N for +/-5%% at 95%% CI    : %d\n', ...
    N_required_conservative);