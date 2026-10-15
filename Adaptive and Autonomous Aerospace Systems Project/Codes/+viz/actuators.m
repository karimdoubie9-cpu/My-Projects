function fig = actuators(runs, tags, UAV, opts)
%ACTUATORS  Rotor commands against the physical envelope, and control effort.
%
%   viz.actuators({R}, {'Baseline'}, UAV)
%   viz.actuators(runs, tags, UAV, events=E)
%
%   Task 2.2 asks for the gains to be tuned "taking into account constraints
%   from actuator dynamics (bandwidth, saturation)". This is the figure that
%   closes that loop: it shows what the tuned controller actually DEMANDS of
%   the rotors, against what they can deliver.
%
%   PANEL 1 - rotor speed envelope. Each rotor obeys T_r = k_f*Omega_r^2 and
%   Omega_r is bounded by [Omega_min, Omega_max]. Plotting all six traces for
%   two controllers is unreadable, so each run is drawn as a BAND spanning
%   min..max across its six rotors: the band's width is how asymmetric the
%   allocation is, and any excursion past a red limit line is saturation - at
%   which point the wrench the vehicle receives is no longer the wrench the
%   controller asked for, and every linear conclusion is void.
%
%   PANEL 2 - collective thrust T_c against the vehicle weight. In steady
%   level flight these coincide; the gap is what is being spent on
%   accelerating and on fighting the disturbance.
%
%   PANEL 3 - body moments from the applied wrench w_a(4:6). Moments come from
%   DIFFERENTIAL thrust, so this is the quantity that actually consumes the
%   spread between the rotors seen in panel 1.

arguments
    runs        cell
    tags        cell
    UAV         struct
    opts.events double = zeros(0,2)
    opts.name   (1,:) char = 'Actuator usage'
end

s   = viz.style();
fig = figure('Name', opts.name, 'Color','w', 'Position',[100 100 900 760]);

have = @(R,f) isfield(R,f) && ~isempty(R.(f));

%% ---- Panel 1: rotor speed envelope ---------------------------------
subplot(3,1,1); hold on; grid on
for k = 1:numel(runs)
    R = runs{k};
    if ~have(R,'omega_c'), continue, end
    col = s.colors{1+mod(k-1,numel(s.colors))};
    lo  = min(R.omega_c, [], 2);
    hi  = max(R.omega_c, [], 2);
    t   = R.omega_c_t;
    fill([t; flipud(t)], [hi; flipud(lo)], col, 'FaceAlpha', 0.25, ...
         'EdgeColor','none', 'DisplayName', sprintf('%s (rotor spread)', tags{k}));
    plot(t, hi, 'LineWidth', 1.7, 'Color', col, 'HandleVisibility','off');
    plot(t, lo, 'LineWidth', 1.7, 'Color', col, 'HandleVisibility','off');
end
yline(UAV.Omega_max, '--', '\Omega_{max}', 'Color',[0.8 0 0], 'LineWidth',1.8, ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');
yline(UAV.Omega_min, '--', '\Omega_{min}', 'Color',[0.8 0 0], 'LineWidth',1.8, ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');
markEvents(opts.events);
ylabel('\Omega_r  [rad/s]'); title('Commanded rotor speeds vs actuator envelope');
legend('Location','best');

%% ---- Panel 2: collective thrust -------------------------------------
subplot(3,1,2); hold on; grid on
for k = 1:numel(runs)
    R = runs{k};
    if ~have(R,'T_c'), continue, end
    plot(R.T_c_t, R.T_c, 'LineWidth', s.lw, ...
         'Color', s.colors{1+mod(k-1,numel(s.colors))}, 'DisplayName', tags{k});
end
yline(UAV.m*9.81, ':', 'bare weight', 'Color',[0.3 0.3 0.3], 'HandleVisibility','off');
markEvents(opts.events);
ylabel('T_c  [N]'); title('Collective thrust command'); legend('Location','best');

%% ---- Panel 3: body moments -----------------------------------------
subplot(3,1,3); hold on; grid on
lbl = {'\tau_x','\tau_y','\tau_z'};
for k = 1:numel(runs)
    R = runs{k};
    if ~have(R,'w_a'), continue, end
    st = ternary(k == 1, '-', '--');
    set(gca,'ColorOrderIndex',1);
    plot(R.w_a_t, R.w_a(:,4:6), st, 'LineWidth', 1.7);
end
markEvents(opts.events);
xlabel('t  [s]'); ylabel('\tau  [N m]');
title('Applied body moments (solid: first case, dashed: second)');
legend(lbl, 'Location','best');

%% ---- Envelope utilisation report ------------------------------------
fprintf('\n=== Actuator envelope (limits %.0f - %.0f rad/s) ===\n', ...
        UAV.Omega_min, UAV.Omega_max);
fprintf('%-22s %10s %10s %10s %12s\n', ...
        'case','min \Omega','max \Omega','max spread','saturated');
fprintf('%s\n', repmat('-',1,70));
for k = 1:numel(runs)
    R = runs{k};
    if ~have(R,'omega_c'), continue, end
    W   = R.omega_c;
    sat = nnz(any(W > UAV.Omega_max | W < UAV.Omega_min, 2));
    fprintf('%-22s %10.1f %10.1f %10.1f %11.1f%%\n', tags{k}, ...
            min(W(:)), max(W(:)), max(max(W,[],2) - min(W,[],2)), ...
            100*sat/size(W,1));
    if sat > 0
        fprintf(['    WARNING: the rotors saturate for %.1f%% of the run. The\n' ...
                 '    delivered wrench is not the commanded one there, so the\n' ...
                 '    linear analysis does not describe this case.\n'], ...
                 100*sat/size(W,1));
    end
end
end

% ---------------------------------------------------------------------
function markEvents(events)
for e = 1:size(events,1)
    xline(events(e,1), '--', 'Color',[0.8 0.2 0.2], 'HandleVisibility','off');
end
end

function out = ternary(c,a,b)
if c, out = a; else, out = b; end
end
