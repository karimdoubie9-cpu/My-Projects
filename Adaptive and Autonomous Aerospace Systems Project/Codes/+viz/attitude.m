function fig = attitude(runs, tags, opts)
%ATTITUDE  Inner-loop response: Euler angles, heading tracking, body rates.
%
%   viz.attitude({R}, {'Baseline'})
%   viz.attitude(runs, tags, events=E)
%
%   The position controller can only produce horizontal force by TILTING, so
%   roll and pitch are not incidental - they are the actuation channel for
%   horizontal motion. Reading them alongside the position error is what
%   distinguishes the two failure modes that look identical from outside:
%
%     * large position error WITH large tilt   -> the outer loop is asking for
%       force and the vehicle is delivering it; the error is thrust or
%       bandwidth limited.
%     * large position error WITHOUT tilt      -> the outer loop is not
%       commanding it; the fault is in the position loop or its reference.
%
%   PANEL 1  roll and pitch. Small in normal flight; the excursion during a
%            gust is the horizontal force being generated to resist it.
%   PANEL 2  yaw against its command yaw_d. Heading is commanded INDEPENDENTLY
%            of position, so this is the one attitude channel with a genuine
%            reference to track rather than a value implied by the outer loop.
%   PANEL 3  body angular rates. These should settle far faster than the
%            position error - that visible time-scale gap IS the cascade
%            assumption the whole architecture rests on.

arguments
    runs        cell
    tags        cell
    opts.events double = zeros(0,2)
    opts.name   (1,:) char = 'Attitude response'
end

s   = viz.style();
fig = figure('Name', opts.name, 'Color','w', 'Position',[110 110 900 760]);
have = @(R,f) isfield(R,f) && ~isempty(R.(f));

% quat2eul(...,'ZYX') orders the output [yaw pitch roll]
IYAW = 1; IPITCH = 2; IROLL = 3;

%% ---- Panel 1: roll and pitch ---------------------------------------
subplot(3,1,1); hold on; grid on
for k = 1:numel(runs)
    R = runs{k};
    if ~have(R,'euler'), continue, end
    col = s.colors{1+mod(k-1,numel(s.colors))};
    plot(R.euler_t, rad2deg(R.euler(:,IROLL)),  '-',  'LineWidth',1.8, ...
         'Color', col, 'DisplayName', sprintf('%s roll', tags{k}));
    plot(R.euler_t, rad2deg(R.euler(:,IPITCH)), '--', 'LineWidth',1.8, ...
         'Color', col, 'DisplayName', sprintf('%s pitch', tags{k}));
end
markEvents(opts.events);
ylabel('angle  [deg]'); title('Roll (solid) and pitch (dashed) - the horizontal actuation channel');
legend('Location','best');

%% ---- Panel 2: heading ----------------------------------------------
subplot(3,1,2); hold on; grid on
drewRef = false;
for k = 1:numel(runs)
    R = runs{k};
    if ~have(R,'euler'), continue, end
    if ~drewRef && have(R,'yaw_d')
        plot(R.yaw_d_t, rad2deg(unwrap(R.yaw_d)), s.ref_style, ...
             'Color', s.ref_color, 'LineWidth', s.lw_ref, 'DisplayName','\psi_d');
        drewRef = true;
    end
    plot(R.euler_t, rad2deg(unwrap(R.euler(:,IYAW))), 'LineWidth', s.lw, ...
         'Color', s.colors{1+mod(k-1,numel(s.colors))}, 'DisplayName', tags{k});
end
markEvents(opts.events);
ylabel('\psi  [deg]'); title('Heading tracking'); legend('Location','best');

%% ---- Panel 3: body rates -------------------------------------------
subplot(3,1,3); hold on; grid on
for k = 1:numel(runs)
    R = runs{k};
    if ~have(R,'omegab'), continue, end
    st = ternary(k == 1, '-', '--');
    set(gca,'ColorOrderIndex',1);
    plot(R.omegab_t, R.omegab, st, 'LineWidth', 1.6);
end
markEvents(opts.events);
xlabel('t  [s]'); ylabel('\omega_b  [rad/s]');
title('Body angular rates (solid: first case, dashed: second)');
legend({'p','q','r'}, 'Location','best');

%% ---- Numbers --------------------------------------------------------
fprintf('\n=== Attitude excursions ===\n');
fprintf('%-22s %13s %13s %11s %13s\n', ...
        'case','max |roll|','max |pitch|','max |omega|','yaw RMS err');
fprintf('%s\n', repmat('-',1,72));
for k = 1:numel(runs)
    R = runs{k};
    if ~have(R,'euler'), continue, end
    yerr = NaN;
    if have(R,'yaw_d')
        yd   = interp1(R.yaw_d_t, unwrap(R.yaw_d), R.euler_t, 'linear','extrap');
        yerr = rad2deg(sqrt(mean((unwrap(R.euler(:,IYAW)) - yd).^2)));
    end
    wmax = NaN;
    if have(R,'omegab'), wmax = max(vecnorm(R.omegab,2,2)); end
    % Plain 'deg' rather than char(176): the degree sign does not survive
    % every console encoding, and a mangled glyph in a results table is
    % worse than a slightly longer label.
    fprintf('%-22s %9.2f deg %9.2f deg %11.3f %9.2f deg\n', tags{k}, ...
            rad2deg(max(abs(R.euler(:,IROLL)))), ...
            rad2deg(max(abs(R.euler(:,IPITCH)))), ...
            wmax, yerr);
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
